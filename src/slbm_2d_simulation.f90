module slbm_2d_simulation

    use slbm_2d_kinds,   only : wp, ip
    use slbm_2d_const,   only : STOP_COND_TIME
    use slbm_2d_const,   only : STOP_COND_STEADY
    use slbm_2d_config,  only : config_t
    use slbm_2d_domain,  only : domain_t

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
            !call this % print_info(iter, time)
        end do
    end subroutine run


    subroutine print_info(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(in)         :: iter
        real(wp), intent(in)            :: time
        print *, iter, time
    end subroutine print_info


    subroutine predictor(this, time)
        class(simulation_t), intent(inout) :: this
        real(wp), intent(in)               :: time
        !print *, 'predictor'
    end subroutine predictor


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


    subroutine set_bndry_cond(this, time)
        class(simulation_t), intent(inout) :: this
        real(wp), intent(in)               :: time

        ! left boundary

        ! right boundary

        ! top boundary

        ! bottom boundary

    end subroutine set_bndry_cond


    subroutine incr_iter_time(this, iter, time)
        class(simulation_t), intent(in) :: this
        integer(ip), intent(inout)      :: iter
        real(ip), intent(inout)         :: time
        iter = iter + 1_ip
        time = time + (this % config % dt)
    end subroutine incr_iter_time

    
end module slbm_2d_simulation
