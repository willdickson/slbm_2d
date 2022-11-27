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

    use slbm_2d_bndry,    only : bndry_t
    use slbm_2d_bndry,    only : side_id_to_name
    use slbm_2d_config,   only : config_t
    use slbm_2d_vector,   only : vector_t
    use slbm_2d_domain,   only : domain_t
    use slbm_2d_density,  only : density_t
    use slbm_2d_velocity, only : velocity_t

    implicit none
    private

    integer(ip), parameter :: PRED_STEP = 1
    integer(ip), parameter :: CORR_STEP = 2

    type, public :: simulation_t
        type(config_t)  :: config  ! configuration data
        type(domain_t)  :: domain  ! simulation domain
        real(wp)        :: err
    contains
        private
        procedure, public  :: run
        procedure          :: predictor
        procedure          :: corrector
        procedure          :: set_bndry_cond
        procedure          :: check_stop_cond
        procedure          :: incr_iter_time
        procedure          :: equilib_func
        procedure          :: steady_conv_test
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
        simulation % domain = domain_t(config)
    end function simulation_constructor


    subroutine run(this)
        class(simulation_t), intent(inout) :: this
        real(wp)    :: time = 0.0_wp  ! current time in the simulation
        integer(ip) :: iter = 0_ip    ! current iteration number
        logical     :: done = .false. ! flag indicating completion  

        do while ( .not. done )
            call this % incr_iter_time(iter, time)
            call this % predictor(time)
            call this % set_bndry_cond(time, PRED_STEP)
            call this % corrector(time)
            call this % set_bndry_cond(time, CORR_STEP)
            call this % check_stop_cond(time, done)
            call this % print_info(iter, time)
        end do
    end subroutine run


    subroutine predictor(this, time)
        class(simulation_t), intent(inout), target :: this
        real(wp), intent(in)                       :: time

        real(wp), parameter  :: ex(LATTICE_Q) = LATTICE_E % x
        real(wp), parameter  :: ey(LATTICE_Q) = LATTICE_E % y
        real(wp), pointer    :: ro_l(:,:)  ! alias for last density
        real(wp), pointer    :: ro_p(:,:)  ! alias for pred density 
        real(wp), pointer    :: ro_c(:,:)  ! alias for curr density 
        real(wp), pointer    :: ux_l(:,:)  ! alias for last velocity x-component
        real(wp), pointer    :: uy_l(:,:)  ! alias for last velocity y-component
        real(wp), pointer    :: ux_p(:,:)  ! alias for pred velocity x-component
        real(wp), pointer    :: uy_p(:,:)  ! alias for pred velocity y-component
        real(wp), pointer    :: ux_c(:,:)  ! alias for curr velocity x-component
        real(wp), pointer    :: uy_c(:,:)  ! alias for curr velocity y-component
        real(wp)             :: feq        ! equilibrium distribution
        integer(ip)          :: num_x      ! size of mesh in x dimension
        integer(ip)          :: num_y      ! size of mesh in y dimension
        integer(ip)          :: i, j, k    ! loop indices
        integer(ip)          :: ie, je     ! loop indices

        ! Get mesh size along x and y dimensions
        num_x = this % config % num_x
        num_y = this % config % num_y

        ! Get aliases for equilibrium and predictor values of density and 
        ! the x and y components of the velocity field 
        ro_l => this % domain % density % last
        ro_p => this % domain % density % pred
        ro_c => this % domain % density % curr 

        ux_l => this % domain % velocity % last % x
        ux_p => this % domain % velocity % pred % x
        ux_c => this % domain % velocity % curr % x

        uy_l => this % domain % velocity % last % y
        uy_c => this % domain % velocity % curr % y
        uy_p => this % domain % velocity % pred % y

        ro_p= 0.0_wp
        ux_p= 0.0_wp
        uy_p= 0.0_wp

        ! Perform predictor update
        do j = 2, num_y-1_ip
            do i = 2, num_x-1_ip
                ro_l(i,j) = ro_c(i,j)
                ux_l(i,j) = ux_c(i,j)
                uy_l(i,j) = uy_c(i,j)
                do k = 1, LATTICE_Q 
                    ie = i - nint(ex(k))
                    je = j - nint(ey(k))
                    feq = this % equilib_func(k,ie,je,PRED_STEP) 
                    ro_p(i,j) = ro_p(i,j) + feq
                    ux_p(i,j) = ux_p(i,j) + feq*ex(k)
                    uy_p(i,j) = uy_p(i,j) + feq*ey(k)
                end do
                ux_p(i,j) = ux_p(i,j)/ro_p(i,j)
                uy_p(i,j) = uy_p(i,j)/ro_p(i,j)
            end do
        end do

    end subroutine predictor


    subroutine set_bndry_cond(this, time, step)
        class(simulation_t), intent(inout), target :: this
        real(wp), intent(in)                       :: time
        integer(ip), intent(in)                    :: step

        integer(ip)              :: num_x       ! mesh size in x dim
        integer(ip)              :: num_y       ! mesh size in x dim
        integer(ip)              :: i1          ! 1st x index
        integer(ip)              :: i2          ! 2nd x index
        integer(ip)              :: ik          ! x offset direction index
        integer(ip)              :: j1          ! 1st y index
        integer(ip)              :: j2          ! 2nd y index
        integer(ip)              :: jk          ! y offset direction index
        integer(ip)              :: side_id     ! boundry side id
        integer(ip)              :: bindx       ! boundry index
        type(bndry_t), pointer   :: bndry       ! alias for boundary
        type(vector_t), pointer  :: u(:,:)      ! alias for velocty
        real(wp), pointer        :: ro(:,:)    ! alias for density 
        real(wp), pointer        :: ro1(:,:)   ! alias, 1st offset for density 
        real(wp), pointer        :: ro2(:,:)   ! alias, 2nd offset for density 

        num_x = this % config % num_x
        num_y = this % config % num_y

        select case (step)
        case (PRED_STEP)
            u  => this % domain % velocity % pred
            ro => this % domain % density  % pred
        case (CORR_STEP)
            u  => this % domain % velocity % curr
            ro => this % domain % density  % curr 
        end select

        do bindx = 1, NUM_BNDRY 
            ! Get boundary data
            side_id = BNDRY_SIDE_IDS(bindx)
            call this % config % get_bndry(side_id, bndry)
            call bndry % get_indices(num_x, num_y, i1, i2, j1, j2)
            call bndry % get_offsets(ik, jk)

            ! Set velocity boundary conditions
            select case (bndry % cond_id)
            case (BNDRY_COND_INFLOW, BNDRY_COND_MOVING, BNDRY_COND_NOSLIP)
                u(i1:i2, j1:j2) = bndry % velocity
            case (BNDRY_COND_OUTFLOW)
                u(i1:i2, j1:j2) = u(i1+ik:i2+ik, j1+jk:j2+jk)
            case (BNDRY_COND_SLIP)
                select case( bndry % side_id)
                case (BNDRY_SIDE_LEFT, BNDRY_SIDE_RIGHT)
                    u(i1:i2, j1:j2) % x = 0.0_wp
                    u(i1:i2, j1:j2) % y = u(i1+ik:i2+ik, j1+jk:j2+jk) % y 
                case (BNDRY_SIDE_TOP, BNDRY_SIDE_BOTTOM)
                    u(i1:i2, j1:j2) % x = u(i1+ik:i2+ik, j1+jk:j2+jk) % x 
                    u(i1:i2, j1:j2) % y = 0.0_wp
                end select
            end select

            ! Set density boundary conditions
            ro1 => ro( i1+1*ik : i2+1*ik, j1+1*jk : j2+1*jk)
            ro2 => ro( i1+2*ik : i2+2*ik, j1+2*jk : j2+2*jk) 
            ro(i1:i2, j1:j2) = (4.0_wp*ro1 - ro2)/3.0_wp
        end do
    end subroutine set_bndry_cond


    subroutine corrector(this, time)
        class(simulation_t), intent(inout), target :: this
        real(wp), intent(in)                       :: time

        real(wp), parameter  :: ex(LATTICE_Q) = LATTICE_E % x
        real(wp), parameter  :: ey(LATTICE_Q) = LATTICE_E % y
        real(wp), pointer    :: ro_l(:,:)  ! alias for last step density
        real(wp), pointer    :: ro_p(:,:)  ! alias for pred step density   
        real(wp), pointer    :: ro_c(:,:)  ! alias for curr step density 
        real(wp), pointer    :: ux_l(:,:)  ! alias for last step velocity x-component
        real(wp), pointer    :: uy_l(:,:)  ! alias for last step velocity y-component
        real(wp), pointer    :: ux_p(:,:)  ! alias for pred step velocity x-component
        real(wp), pointer    :: uy_p(:,:)  ! alias for pred step velocity y-component
        real(wp), pointer    :: ux_c(:,:)  ! alias for curr step velocity x-component
        real(wp), pointer    :: uy_c(:,:)  ! alias for curr step velocity y-component
        real(wp)             :: tau        ! Relaxation parameter
        real(wp)             :: feq        ! equilibrium distribution 
        integer(ip)          :: num_x      ! size of mesh in x dimension
        integer(ip)          :: num_y      ! size of mesh in y dimension
        integer(ip)          :: i, j, k    ! loop indices
        integer(ip)          :: ie, je     

        ! Get mesh size along x and y dimensions
        tau   = this % config % tau
        num_x = this % config % num_x
        num_y = this % config % num_y

        ! Get aliases for equilibrium and predictor values of density and 
        ! the x and y components of the velocity field 
        ro_l => this % domain % density % last 
        ro_p => this % domain % density % pred
        ro_c => this % domain % density % curr 

        ux_l => this % domain % velocity % last % x
        ux_p => this % domain % velocity % pred % x
        ux_c => this % domain % velocity % curr % x

        uy_l => this % domain % velocity % last % y
        uy_p => this % domain % velocity % pred % y
        uy_c => this % domain % velocity % curr % y

        do j = 2, num_y-1
            do i = 2, num_x-1
                ux_c(i,j) = ro_p(i,j)*ux_p(i,j)
                uy_c(i,j) = ro_p(i,j)*uy_p(i,j)
                do k = 1, LATTICE_Q
                    ie = i - nint(ex(k))
                    je = j - nint(ey(k))
                    feq = this % equilib_func(k,ie,je,CORR_STEP)
                    ux_c(i,j) = ux_c(i,j) + (tau - 1.0_wp)*ex(k)*feq
                    uy_c(i,j) = uy_c(i,j) + (tau - 1.0_wp)*ey(k)*feq
                end do
                ux_c(i,j) = (ux_c(i,j) - (tau - 1.0_wp)*ro_l(i,j)*ux_l(i,j))/ro_p(i,j) 
                uy_c(i,j) = (uy_c(i,j) - (tau - 1.0_wp)*ro_l(i,j)*uy_l(i,j))/ro_p(i,j) 
                ro_c(i,j) = ro_p(i,j)
            end do
        end do

    end subroutine corrector


    subroutine check_stop_cond(this, time, done) 
        class(simulation_t), intent(inout) :: this
        real(wp), intent(in)               :: time
        logical, intent(inout)             :: done
        real(wp)                           :: err
        select case (this % config % stop_cond)
        case (STOP_COND_TIME)
            if (time >= this % config % stop_time) then
                done = .true.
            end if
        case (STOP_COND_STEADY)
            call this % steady_conv_test(err)
            if (err <= this % config % stop_etol) then
                done = .true.
            end if
            this % err = err
        case default
            print *, 'unknown stop condition'
            stop
        end select
    end subroutine check_stop_cond


    subroutine incr_iter_time(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(inout)      :: iter
        real(wp), intent(inout)         :: time
        iter = iter + 1_ip
        time = time + (this % config % dt)
    end subroutine incr_iter_time


    pure function equilib_func(this, k, i, j, step) result(equilib)
        class(simulation_t), intent(in), target :: this
        integer(ip), intent(in)                 :: k
        integer(ip), intent(in)                 :: i
        integer(ip), intent(in)                 :: j
        integer(ip), intent(in)                 :: step
        real(wp)                                :: equilib 

        ! Constant used in equilibrium calculation
        real(wp), parameter  :: k1 = 1.0_wp/CS2
        real(wp), parameter  :: k2 = 1.0_wp/(2.0_wp*CS4)
        real(wp), parameter  :: k3 = 1.0_wp/(2.0_wp*CS2)

        real(wp)             :: ro  ! density at i,j
        real(wp)             :: ux  ! velocity x-component at i,j
        real(wp)             :: uy  ! velocity y-component at i,j
        real(wp)             :: ex  ! k-th lattice velocity x-component at i,j
        real(wp)             :: ey  ! k-th lattice velocity y-component at i,j
        real(wp)             :: wt  ! k-th lattice weight
        real(wp)             :: uu  ! squared magnitude of velocity
        real(wp)             :: eu  ! lattice velocity (ex,ey) to velocity (ux,uy)
        real(wp)             :: eu2 ! square of eu

        ! Get denisty and velocity vector components
        select case (step)
        case (PRED_STEP)
            ro = this % domain % density  % curr(i,j)
            ux = this % domain % velocity % curr(i,j) % x
            uy = this % domain % velocity % curr(i,j) % y
        case (CORR_STEP)
            ro = this % domain % density  % pred(i,j)
            ux = this % domain % velocity % pred(i,j) % x
            uy = this % domain % velocity % pred(i,j) % y
        end select

        ! Get lattice velocity components and weight
        ex = LATTICE_E(k) % x
        ey = LATTICE_E(k) % y
        wt = LATTICE_W(k)

        ! Compute equilibrium value
        uu = ux**2 + uy**2
        eu = ex*ux + ey*uy
        eu2 = eu**2
        equilib = ro*wt*(1.0_wp + k1*eu + k2*eu2 - k3*uu)
    end function equilib_func

    subroutine steady_conv_test(this, max_err)
        class(simulation_t), intent(in), target :: this
        real(wp), intent(out)                   :: max_err

        real(wp), pointer  :: ux_l(:,:)
        real(wp), pointer  :: uy_l(:,:)
        real(wp), pointer  :: ux_c(:,:)
        real(wp), pointer  :: uy_c(:,:)
        real(wp)           :: umag_l
        real(wp)           :: umag_c
        real(wp)           :: err
        integer(ip)        :: num_x
        integer(ip)        :: num_y
        integer(ip)        :: i, j

        num_x = this % config % num_x
        num_y = this % config % num_y

        ux_l => this % domain % velocity % last % x
        ux_c => this % domain % velocity % curr % x

        uy_l => this % domain % velocity % last % y
        uy_c => this % domain % velocity % curr % y

        do j = 2, num_y-1
            do i = 2, num_x-1
                umag_l = sqrt(ux_l(i,j)**2 + uy_l(i,j)**2)
                umag_c = sqrt(ux_c(i,j)**2 + uy_c(i,j)**2)
                err = abs(umag_c - umag_l)/umag_l
                if (err > max_err) then
                    max_err = err
                end if
            end do
        end do

    end subroutine steady_conv_test


    subroutine print_info(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(in)         :: iter
        real(wp), intent(in)            :: time
        if (this % config % stop_cond == STOP_COND_TIME) then
            print *, iter, time
        else
            print *, iter, time, this % err
        end if
    end subroutine print_info
    

    
end module slbm_2d_simulation
