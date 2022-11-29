module slbm_2d_simulation

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_const,    only : CS2
    use slbm_2d_const,    only : CS4

    use slbm_2d_const,    only : LATTICE_Q
    use slbm_2d_const,    only : LATTICE_E
    use slbm_2d_const,    only : LATTICE_W

    use slbm_2d_const,    only : NUM_BNDRY
    use slbm_2d_const,    only : BNDRY_SIDE_LEFT
    use slbm_2d_const,    only : BNDRY_SIDE_RIGHT
    use slbm_2d_const,    only : BNDRY_SIDE_TOP
    use slbm_2d_const,    only : BNDRY_SIDE_BOTTOM
    use slbm_2d_const,    only : BNDRY_SIDE_UNKNOWN
    use slbm_2d_const,    only : BNDRY_SIDE_NAMES
    use slbm_2d_const,    only : BNDRY_SIDE_IDS

    use slbm_2d_const,    only : BNDRY_COND_INFLOW
    use slbm_2d_const,    only : BNDRY_COND_MOVING
    use slbm_2d_const,    only : BNDRY_COND_OUTFLOW
    use slbm_2d_const,    only : BNDRY_COND_NOSLIP
    use slbm_2d_const,    only : BNDRY_COND_SLIP

    use slbm_2d_const,    only : STOP_COND_TIME
    use slbm_2d_const,    only : STOP_COND_STEADY

    use slbm_2d_state,    only : state_t
    use slbm_2d_bndry,    only : bndry_t
    use slbm_2d_bndry,    only : side_id_to_name
    use slbm_2d_config,   only : config_t
    use slbm_2d_vector,   only : vector_t
    use slbm_2d_vector,   only : dot
    use slbm_2d_vector,   only : mag 

    use stdlib_io_npy,    only : save_npy

    implicit none
    private

    integer(ip), parameter :: LAST_STATE = 1
    integer(ip), parameter :: PRED_STATE = 2
    integer(ip), parameter :: CURR_STATE = 3

    type, public :: simulation_t
        type(config_t)  :: config  ! configuration data
        type(state_t)   :: last    ! state from last step
        type(state_t)   :: pred    ! state from predictor step
        type(state_t)   :: curr    ! state from current step
    contains
        private
        procedure, public  :: run
        procedure          :: set_bndry_cond
        procedure          :: incr_time
        procedure          :: predictor
        procedure          :: corrector
        procedure          :: check_stop_cond
        procedure          :: equilib_func
        procedure          :: steady_conv_test
        procedure          :: check_save_dir
        procedure          :: save_data
        procedure          :: print_info
    end type simulation_t

    interface simulation_t
        procedure simulation_constructor
    end interface simulation_t
    
contains

    function simulation_constructor(config) result(simulation)
        type(config_t), intent(in) :: config
        type(simulation_t)         :: simulation
        simulation % config = config

        simulation % last = state_t(config)
        call simulation % set_bndry_cond(LAST_STATE)

        simulation % pred = state_t(config)
        call simulation % set_bndry_cond(PRED_STATE)

        simulation % curr = state_t(config)
        call simulation % set_bndry_cond(CURR_STATE) 
    end function simulation_constructor


    subroutine run(this)
        class(simulation_t), intent(inout) :: this
        logical     :: done  = .false. ! flag indicating completion  
        integer(ip) :: iter  = 0_ip    ! current iteration number
        real(wp)    :: time  = 0.0_wp  ! current time in the simulation
        real(wp)    :: serr  = 0.0_wp  ! steady state error  

        call this % check_save_dir()

        do while ( .not. done )
            call this % incr_time(iter, time)
            call this % predictor()
            call this % set_bndry_cond(PRED_STATE)
            call this % corrector()
            call this % set_bndry_cond(CURR_STATE)
            call this % check_stop_cond(time, done, serr)
            call this % save_data(iter, time)
            call this % print_info(iter, time, serr)
        end do
    end subroutine run


    subroutine predictor(this)
        class(simulation_t), intent(inout), target :: this

        type(state_t), pointer :: last     ! alias for last state
        type(state_t), pointer :: pred     ! alias for pred state
        type(state_t), pointer :: curr     ! alias for curr state

        real(wp)               :: feq      ! equilibrium distribution
        integer(ip)            :: num_x    ! size of mesh in x dimension
        integer(ip)            :: num_y    ! size of mesh in y dimension
        integer(ip)            :: i, j, k  ! loop indices

        ! Get aliases for last, pred and curr states (less verbose)
        last => this % last
        pred => this % pred
        curr => this % curr

        ! Get mesh size along x and y dimensions
        num_x = this % config % num_x
        num_y = this % config % num_y

        ! Perform predictor update
        pred = 0.0_wp

        do j = 2, num_y-1_ip
            do i = 2, num_x-1_ip
                last % rho(i,j) = curr % rho(i,j)
                last % u(i,j)   = curr % u(i,j)
                do k = 1, LATTICE_Q 
                    feq = this % equilib_func(k,i,j,CURR_STATE)
                    pred % rho(i,j) = pred % rho(i,j) + feq
                    pred % u(i,j)   = pred % u(i,j)   + feq*LATTICE_E(k)
                end do
                pred % u(i,j) = pred % u(i,j) / pred % rho(i,j)
            end do
        end do
    end subroutine predictor


    subroutine set_bndry_cond(this, state_id)
        class(simulation_t), intent(inout), target :: this
        integer(ip), intent(in)                    :: state_id

        type(state_t), pointer   :: s           ! alias for state
        integer(ip)              :: num_x       ! mesh size in x dim
        integer(ip)              :: num_y       ! mesh size in x dim
        integer(ip)              :: i1          ! 1st x index
        integer(ip)              :: i2          ! 2nd x index
        integer(ip)              :: ik          ! x offset direction index
        integer(ip)              :: j1          ! 1st y index
        integer(ip)              :: j2          ! 2nd y index
        integer(ip)              :: jk          ! y offset direction index
        integer(ip)              :: side_id     ! boundry side id
        integer(ip)              :: bndry_index ! boundry index
        type(bndry_t), pointer   :: bndry       ! alias for boundary
        real(wp), pointer        :: rho1(:,:)   ! alias, 1st offset for density 
        real(wp), pointer        :: rho2(:,:)   ! alias, 2nd offset for density 

        select case(state_id)
        case (LAST_STATE)
            s => this % last
        case (PRED_STATE)
            s => this % pred
        case (CURR_STATE)
            s => this % curr
        case default
            print *, 'set_bndry_cond: unknown state_id'
            stop
        end select

        num_x = this % config % num_x
        num_y = this % config % num_y

        do bndry_index = 1, NUM_BNDRY 
            ! Get boundary data
            side_id = BNDRY_SIDE_IDS(bndry_index)
            call this % config % get_bndry(side_id, bndry)
            call bndry % get_indices(num_x, num_y, i1, i2, j1, j2)
            call bndry % get_offsets(ik, jk)

            ! Set velocity boundary conditions
            select case (bndry % cond_id)
            case (BNDRY_COND_INFLOW, BNDRY_COND_MOVING, BNDRY_COND_NOSLIP)
                s % u(i1:i2, j1:j2) = bndry % velocity
            case (BNDRY_COND_OUTFLOW)
                s % u(i1:i2, j1:j2) = s % u(i1+ik:i2+ik, j1+jk:j2+jk)
            case (BNDRY_COND_SLIP)
                select case( bndry % side_id)
                case (BNDRY_SIDE_LEFT, BNDRY_SIDE_RIGHT)
                    s % u(i1:i2, j1:j2) % x = 0.0_wp
                    s % u(i1:i2, j1:j2) % y = s % u(i1+ik:i2+ik, j1+jk:j2+jk) % y 
                case (BNDRY_SIDE_TOP, BNDRY_SIDE_BOTTOM)
                    s % u(i1:i2, j1:j2) % x = s % u(i1+ik:i2+ik, j1+jk:j2+jk) % x 
                    s % u(i1:i2, j1:j2) % y = 0.0_wp
                end select
            end select

            ! Set density boundary conditions
            rho1 => s % rho( i1+1*ik : i2+1*ik, j1+1*jk : j2+1*jk)
            rho2 => s % rho( i1+2*ik : i2+2*ik, j1+2*jk : j2+2*jk) 
            s % rho(i1:i2, j1:j2) = (4.0_wp*rho1 - rho2)/3.0_wp
        end do
    end subroutine set_bndry_cond


    subroutine corrector(this)
        class(simulation_t), intent(inout), target :: this

        type(state_t), pointer :: last       ! alias for last state 
        type(state_t), pointer :: pred       ! alias for pred state
        type(state_t), pointer :: curr       ! alias for curr state
        real(wp)               :: feq        ! equilibrium distribution
        real(wp)               :: tau        ! Relaxation parameter
        integer(ip)            :: num_x      ! size of mesh in x dimension
        integer(ip)            :: num_y      ! size of mesh in y dimension
        integer(ip)            :: i, j, k    ! loop indices

        ! Get aliases for last, pred and curr states (less verbose)
        last => this % last
        pred => this % pred
        curr => this % curr

        ! Get mesh size along x and y dimensions
        tau   = this % config % tau
        num_x = this % config % num_x
        num_y = this % config % num_y

        do j = 2, num_y-1
            do i = 2, num_x-1
                curr % u(i,j) = pred % rho(i,j) * pred % u(i,j)
                do k = 1, LATTICE_Q
                    feq = this % equilib_func(k,i,j,PRED_STATE)
                    curr % u(i,j) = curr % u(i,j) + (tau-1.0_wp) * LATTICE_E(k)*feq
                end do
                curr % rho(i,j) = pred % rho(i,j)
                curr % u(i,j)   = curr % u(i,j) - (tau-1.0_wp) * last % rho(i,j) * last % u(i,j)
                curr % u(i,j)   = curr % u(i,j) / curr % rho(i,j)
            end do
        end do

    end subroutine corrector


    subroutine check_stop_cond(this, time, done, serr) 
        class(simulation_t), intent(inout) :: this
        real(wp), intent(in)               :: time
        logical, intent(inout)             :: done
        real(wp), intent(out)              :: serr 

        select case (this % config % stop_cond)
        case (STOP_COND_TIME)
            if (time >= this % config % stop_time) then
                done = .true.
            end if
        case (STOP_COND_STEADY)
            serr = 0.0_wp
            call this % steady_conv_test(serr)
            if (serr <= this % config % stop_etol) then
                done = .true.
            end if
        case default
            print *, 'unknown stop condition'
            stop
        end select
    end subroutine check_stop_cond


    subroutine incr_time(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(inout)      :: iter
        real(wp), intent(inout)         :: time
        iter = iter + 1_ip
        time = time + (this % config % dt)
    end subroutine incr_time


    function equilib_func(this, k, i, j, state_id) result(equilib)
        class(simulation_t), intent(in), target :: this
        integer(ip), intent(in)                 :: k
        integer(ip), intent(in)                 :: i
        integer(ip), intent(in)                 :: j
        integer(ip), intent(in)                 :: state_id
        real(wp)                                :: equilib 

        ! Constant used in equilibrium calculation
        real(wp), parameter  :: k1 = 1.0_wp/CS2
        real(wp), parameter  :: k2 = 1.0_wp/(2.0_wp*CS4)
        real(wp), parameter  :: k3 = 1.0_wp/(2.0_wp*CS2)

        type(state_t), pointer   :: state ! alias for state
        integer(ip)              :: ie    ! index offset by lattice vector x-comp
        integer(ip)              :: je    ! index offest by lattice vector y-comp
        type(vector_t)           :: u     ! velocity vector at ie, je
        real(wp)                 :: wt    ! k-th lattice weight
        real(wp)                 :: rho
        real(wp)                 :: uu    ! squared magnitude of velocity
        real(wp)                 :: eu    ! lattice velocity (ex,ey) to velocity (ux,uy)
        real(wp)                 :: eu2   ! square of eu

        select case (state_id)
        case (LAST_STATE)
            state => this % last
        case (PRED_STATE)
            state => this % pred
        case (CURR_STATE)
            state => this % curr
        end select

        ie  = i - nint(LATTICE_E(k) % x)
        je  = j - nint(LATTICE_E(k) % y)
        rho = state % rho(ie,je)
        u   = state % u(ie,je)

        wt  = LATTICE_W(k)
        uu  = dot(u,u)
        eu  = dot(LATTICE_E(k),u)
        eu2 = eu**2 
        equilib = rho*wt*(1.0_wp + k1*eu + k2*eu2 - k3*uu)

    end function equilib_func


    subroutine steady_conv_test(this, max_err)
        class(simulation_t), intent(in), target :: this
        real(wp), intent(out)                   :: max_err

        type(state_t), pointer :: last
        type(state_t), pointer :: curr
        integer(ip)            :: num_x
        integer(ip)            :: num_y
        integer(ip)            :: i,j
        real(wp)               :: mag_last
        real(wp)               :: mag_curr
        real(wp)               :: err

        last => this % last
        curr => this % curr

        num_x = this % config % num_x
        num_y = this % config % num_y

        do j = 2, num_y-1
            do i = 2, num_x-1
                mag_last = mag(last % u(i,j))
                mag_curr = mag(curr % u(i,j))
                err = abs(mag_curr - mag_last)/mag_last
                if (err > max_err) then
                    max_err = err
                end if
            end do
        end do
    end subroutine steady_conv_test


    subroutine check_save_dir(this)
        class(simulation_t), intent(in) :: this
        character(:), allocatable       :: cmd
        cmd = 'mkdir -p ' // this % config % save_dir
        call execute_command_line(cmd)
        cmd = 'rm ' // this % config % save_dir // "/*.npy"
        call execute_command_line(cmd)
    end subroutine check_save_dir


    subroutine save_data(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(in)         :: iter

        integer(ip), parameter          :: num_item = 3 
        integer(ip), save               :: save_cnt = 0_ip
        integer(ip)                     :: num_x
        integer(ip)                     :: num_y
        real(wp), intent(in)            :: time
        real(wp), allocatable           :: out_data(:,:,:)
        character(:), allocatable       :: save_cnt_str
        character(:), allocatable       :: out_filename

        if ( modulo(iter, this % config % save_nstep) == 0) then
            save_cnt = save_cnt + 1
            num_x = this % config % num_x
            num_y = this % config % num_y
            allocate(out_data(num_item, num_x, num_y))
            out_data(1,:,:) = this % curr % rho * CS2 ! pressure
            out_data(2,:,:) = this % curr % u % x     ! velocity x-comp
            out_data(3,:,:) = this % curr % u % y     ! velocity y-comp
            save_cnt_str = get_save_cnt_str(save_cnt)
            out_filename = this % config % save_dir // '/data_' // save_cnt_str // '.npy'
            call save_npy(out_filename, out_data)
            print *, 'saving ' // out_filename
        end if

    end subroutine save_data


    subroutine print_info(this, iter, time, serr)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(in)         :: iter
        real(wp), intent(in)            :: time
        real(wp), intent(in)            :: serr
        if (this % config % stop_cond == STOP_COND_TIME) then
            print *, iter, time
        else
            print *, iter, time, serr
        end if
    end subroutine print_info
    

    function get_save_cnt_str(cnt) result(res)
        integer,intent(in)        :: cnt
        character(:),allocatable  :: res
        character(range(cnt)+2)   :: tmp
        write(tmp,'(i6.6)') cnt 
        res = trim(tmp)
    end function


    function itoa(i) result(res)
        integer,intent(in)        :: i
        character(:),allocatable  :: res
        character(range(i)+2)     :: tmp
        write(tmp,'(i0)') i
        res = trim(tmp)
    end function


end module slbm_2d_simulation
