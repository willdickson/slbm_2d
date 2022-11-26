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
    use slbm_2d_distrib,  only : distrib_t
    use slbm_2d_density,  only : density_t
    use slbm_2d_velocity, only : velocity_t

    implicit none
    private

    type, public :: simulation_t
        type(config_t)  :: config  ! configuration data
        type(domain_t)  :: domain  ! simulation domain
    contains
        private
        procedure, public  :: run
        procedure          :: predictor
        procedure          :: corrector
        procedure          :: set_bndry_cond
        procedure          :: check_stop_cond
        procedure          :: incr_iter_time
        procedure          :: equilib_func
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
            call this % set_bndry_cond(time)
            call this % corrector(time)
            call this % check_stop_cond(time, done)
            call this % print_info(iter, time)
        end do
    end subroutine run


    subroutine predictor(this, time)
        class(simulation_t), intent(inout), target :: this
        real(wp), intent(in)                       :: time

        real(wp), parameter  :: ex(LATTICE_Q) = LATTICE_E % x
        real(wp), parameter  :: ey(LATTICE_Q) = LATTICE_E % y
        real(wp), pointer    :: feq(:,:,:) ! alias for equilibrium distribution
        real(wp), pointer    :: rho(:,:)   ! alias for density 
        real(wp), pointer    :: ux(:,:)    ! alias for velocity x-component
        real(wp), pointer    :: uy(:,:)    ! alias for velocity y-component
        integer(ip)          :: num_x      ! size of mesh in x dimension
        integer(ip)          :: num_y      ! size of mesh in y dimension
        integer(ip)          :: i, j, k    ! loop indices

        num_x = this % config % num_x
        num_y = this % config % num_y

        ! Get aliases for equilibrium and predictor values of density and 
        ! the x and y components of the velocity field 
        feq => this % domain % equilib % val
        rho => this % domain % density % pred
        ux  => this % domain % velocity % pred % x
        uy  => this % domain % velocity % pred % y

        ! Perform predictor update
        do j = 2, num_y-1_ip
            do i = 2, num_x-1_ip
                do k = 1, LATTICE_Q 
                    feq(k,i,j) = this % equilib_func(k,i,j) 
                    rho(i,j) = rho(i,j) + feq(k,i,j)
                    ux(i,j)  = ux(i,j)  + feq(k,i,j)*ex(k)
                    uy(i,j)  = uy(i,j)  + feq(k,i,j)*ey(k)
                end do
                ux(i,j) = ux(i,j)/rho(i,j)
                uy(i,j) = uy(i,j)/rho(i,j)
            end do
        end do

    end subroutine predictor


    subroutine set_bndry_cond(this, time)
        class(simulation_t), intent(inout), target :: this
        real(wp), intent(in)                       :: time

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
        real(wp), pointer        :: feq(:,:,:)  ! alias for equilib distribution
        real(wp), pointer        :: rho(:,:)    ! alias for density 
        real(wp), pointer        :: rho1(:,:)   ! alias, 1st offset for density 
        real(wp), pointer        :: rho2(:,:)   ! alias, 2nd offset for density 

        num_x = this % config % num_x
        num_y = this % config % num_y

        u   => this % domain % velocity % pred
        feq => this % domain % equilib % val
        rho => this % domain % density % pred

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

            ! Set equilibrium distribution boundary conditions
            feq(:,i1:i2, j1:j2) = feq(:,i1+ik:i2+ik, j1+jk:j2+jk)

            ! Set density boundary conditions
            rho1 => rho( i1+1*ik : i2+1*ik, j1+1*jk : j2+1*jk)
            rho2 => rho( i1+2*ik : i2+2*ik, j1+2*jk : j2+2*jk) 
            rho(i1:i2, j1:j2) = (4.0_wp*rho1 - rho2)/3.0_wp

            ! Debug
            ! ---------------------------------------------------
            !print *, side_id_to_name(side_id) 
            !print *, 'side_id     = ', side_id 
            !print *, 'i1, i2, ik  = ', i1, i2, ik 
            !print *, 'j1, j2, jk  = ', j1, j2, jk
            !print *, ''
            ! ---------------------------------------------------
        end do


    end subroutine set_bndry_cond


    subroutine corrector(this, time)
        class(simulation_t), intent(inout) :: this
        real(wp), intent(in)               :: time
        !print *, 'corrector' 
    end subroutine corrector


    subroutine check_stop_cond(this, time, done) 
        class(simulation_t), intent(in) :: this
        real(wp), intent(in)            :: time
        logical, intent(inout)          :: done
        select case (this % config % stop_cond)
        case (STOP_COND_TIME)
            if (time >= this % config % stop_time) then
                done = .true.
            end if
        case (STOP_COND_STEADY)
            print *, 'steady state stop condition not implemented yet'
            stop
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


    pure function equilib_func(this, k, i, j) result(equilib)
        class(simulation_t), intent(in), target :: this
        integer(ip), intent(in)                 :: k
        integer(ip), intent(in)                 :: i
        integer(ip), intent(in)                 :: j
        real(wp)                                :: equilib 

        ! Constant used in equilibrium calculation
        real(wp), parameter  :: k1 = 1.0_wp/CS2
        real(wp), parameter  :: k2 = 1.0_wp/(2.0_wp*CS4)
        real(wp), parameter  :: k3 = 1.0_wp/(2.0_wp*CS2)

        real(wp)             :: rho ! density at i,j
        real(wp)             :: ux  ! velocity x-component at i,j
        real(wp)             :: uy  ! velocity y-component at i,j
        real(wp)             :: ex  ! k-th lattice velocity x-component at i,j
        real(wp)             :: ey  ! k-th lattice velocity y-component at i,j
        real(wp)             :: wt  ! k-th lattice weight
        real(wp)             :: uu  ! squared magnitude of velocity
        real(wp)             :: eu  ! lattice velocity (ex,ey) to velocity (ux,uy)
        real(wp)             :: eu2 ! square of eu

        ! Get, velocity vector components
        ux = this % domain % velocity % curr(i,j) % x
        uy = this % domain % velocity % curr(i,j) % y

        ! Get lattice velocity components and weight
        ex = LATTICE_E(k) % x
        ey = LATTICE_E(k) % y
        wt = LATTICE_W(k)

        ! Compute equilibrium value
        uu = ux*ux + uy*uy
        eu = ex*ux + ey*uy
        eu2 = eu**2
        equilib = rho*wt*(1.0_wp + k1*eu + k2*eu2 - k3*uu)
    end function equilib_func


    subroutine print_info(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(in)         :: iter
        real(wp), intent(in)            :: time
        print *, iter, time
    end subroutine print_info

    
end module slbm_2d_simulation
