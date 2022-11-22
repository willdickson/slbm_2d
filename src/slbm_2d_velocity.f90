module slbm_2d_velocity
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_vector, only : vector_t
    implicit none
    private

    type, public :: velocity_t
        integer(ip)                 :: num_x = 0_ip
        integer(ip)                 :: num_y = 0_ip
        type(vector_t), allocatable :: curr(:,:)
        type(vector_t), allocatable :: pred(:,:)
    contains
        private
        procedure, public  :: deallocate => velocity_deallocate
        final              :: velocity_destructor
    end type velocity_t

    interface velocity_t
        procedure velocity_constructor
    end interface velocity_t

contains

    function velocity_constructor(num_x, num_y) result(velocity)
        integer(ip), intent(in) :: num_x
        integer(ip), intent(in) :: num_y
        type(velocity_t)        :: velocity
        velocity % num_x = num_x
        velocity % num_y = num_y
        allocate(velocity % curr(num_x, num_y))
        allocate(velocity % pred(num_x, num_y))
    end function velocity_constructor

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
