module slbm_2d_velocity
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_init,   only : init_file_t
    use slbm_2d_init,   only : init_const_t
    use slbm_2d_vector, only : vector_t
    use slbm_2d_config, only : config_t
    implicit none
    private

    type, public :: velocity_t
        type(vector_t), allocatable :: curr(:,:)
        type(vector_t), allocatable :: pred(:,:)
    contains
        private
        procedure          :: set_initial_cond
        procedure, public  :: deallocate => velocity_deallocate
        final              :: velocity_destructor
    end type velocity_t

    interface velocity_t
        procedure velocity_constructor
    end interface velocity_t

contains

    function velocity_constructor(config) result(velocity)
        type(config_t), intent(in) :: config
        type(velocity_t)        :: velocity
        allocate(velocity % curr(config % num_x, config % num_y))
        allocate(velocity % pred(config % num_x, config % num_y))
        call velocity % set_initial_cond(config)
    end function velocity_constructor

    subroutine set_initial_cond(this, config)
        class(velocity_t), intent(inout) :: this
        type(config_t), intent(in)       :: config
        select type(init => config % init)
        class is (init_const_t)
            this % curr = init % velocity
            this % pred = vector_t(0.0_wp, 0.0_wp)
        class is (init_file_t)
            print *, "setting initial cond. from file net implemented yet"
            stop
        class default
            print *, "error unknown method for settign initial cond."
            stop
        end select
    end subroutine set_initial_cond

    subroutine velocity_deallocate(this)
        class(velocity_t), intent(inout) :: this
        if ( allocated(this % curr) ) then
            deallocate(this % curr)
        end if
        if ( allocated(this % pred) ) then
            deallocate(this % pred)
        end if
    end subroutine velocity_deallocate

    subroutine velocity_destructor(this)
        type(velocity_t), intent(inout) :: this
        call this % deallocate()
    end subroutine velocity_destructor

end module slbm_2d_velocity
