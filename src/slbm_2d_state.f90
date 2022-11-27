module slbm_2d_state
    use slbm_2d_kinds,      only : wp, ip
    use slbm_2d_config,     only : config_t
    use slbm_2d_vector,     only : vector_t
    use slbm_2d_init,       only : init_file_t
    use slbm_2d_init,       only : init_const_t
    implicit none
    private

    type, public :: state_t
        real(wp), allocatable       :: rho(:,:) ! fluid density 
        type(vector_t), allocatable :: u(:,:)   ! fluid velocity 
    contains
        private
        procedure :: set_init_cond
        procedure :: assign_scalar
        procedure :: assign_state
        generic, public :: assignment(=) => assign_scalar, &
                                            assign_state
    end type

    interface state_t
        procedure :: state_constructor
    end interface state_t

contains

    function state_constructor(config) result(state)
        type(config_t), intent(in) :: config
        type(state_t)              :: state
        allocate(state % rho(config % num_x, config % num_y))
        allocate(state % u(config % num_x, config % num_y))
        call state % set_init_cond(config)
    end function state_constructor

    subroutine set_init_cond(this, config)
        class(state_t), intent(inout) :: this
        type(config_t), intent(in)    :: config
        this % rho = config % density
        select type(init => config % init)
        class is (init_const_t)
            this % u = init % velocity
        class is (init_file_t)
            print *, "setting initial cond. from file net implemented yet"
            stop
        class default
            print *, "error unknown method for settign initial cond."
            stop
        end select
    end subroutine set_init_cond

    pure subroutine assign_scalar(this,scalar)
        class(state_t), intent(inout) :: this
        real(wp), intent(in)          :: scalar
        this % rho = scalar
        this % u   = vector_t(scalar, scalar)
    end subroutine assign_scalar

    pure subroutine assign_state(this, state)
        class(state_t), intent(inout) :: this
        type(state_t), intent(in)     :: state
        this % rho = state % rho
        this % u   = state % u
    end subroutine assign_state

end module slbm_2d_state
