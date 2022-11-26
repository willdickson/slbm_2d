module slbm_2d_simulation

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_const,    only : CS2
    use slbm_2d_const,    only : CS4

    use slbm_2d_const,    only : LATTICE_Q
    use slbm_2d_const,    only : LATTICE_E
    use slbm_2d_const,    only : LATTICE_W

    use slbm_2d_const,    only : BNDRY_SIDE_LEFT
    use slbm_2d_const,    only : BNDRY_SIDE_RIGHT
    use slbm_2d_const,    only : BNDRY_SIDE_TOP
    use slbm_2d_const,    only : BNDRY_SIDE_BOTTOM

    use slbm_2d_const,    only : STOP_COND_TIME
    use slbm_2d_const,    only : STOP_COND_STEADY

    use slbm_2d_bndry,    only : bndry_t
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
        integer(ip)          :: i, j, k    ! loop indices

        ! Get aliases for equilibrium and predictor values of density and 
        ! the x and y components of the velocity field 
        feq => this % domain % equilib % val
        rho => this % domain % density % pred
        ux  => this % domain % velocity % pred % x
        uy  => this % domain % velocity % pred % y

        ! Perform predictor update
        do j = 2, (this % config % num_y - 1)
            do i = 2, (this % config % num_x - 1)
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
        class(simulation_t), intent(inout) :: this
        real(wp), intent(in)               :: time
        type(bndry_t), pointer             :: bndry

        ! left boundary
        call this % config % get_bndry(BNDRY_SIDE_LEFT, bndry)

        ! right boundary
        call this % config % get_bndry(BNDRY_SIDE_RIGHT, bndry)

        ! top boundary
        call this % config % get_bndry(BNDRY_SIDE_TOP, bndry)

        ! bottom boundary
        call this % config % get_bndry(BNDRY_SIDE_BOTTOM, bndry)

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


    function equilib_func(this, k, i, j) result(equilib)
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
