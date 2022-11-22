module slbm_2d_density
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_vector, only : vector_t
    implicit none
    private

    type, public :: density_t
        integer(ip)                 :: num_x = 0_ip
        integer(ip)                 :: num_y = 0_ip
        type(vector_t), allocatable :: curr(:,:)
        type(vector_t), allocatable :: pred(:,:)

    contains
        private
        procedure, public :: deallocate => density_deallocate
        final    :: density_destructor 
    end type density_t

    interface density_t
        procedure density_constructor
    end interface density_t


contains

    function density_constructor(num_x, num_y) result(density)
        integer(ip), intent(in) :: num_x
        integer(ip), intent(in) :: num_y
        type(density_t)         :: density
        allocate(density % curr(num_x, num_y))
        allocate(density % pred(num_x, num_y))
    end function density_constructor

    subroutine density_deallocate(this)
        class(density_t), intent(inout) :: this
        if ( allocated(this % curr) ) then
            deallocate(this % curr)
        end if
        if ( allocated(this % pred) ) then
            deallocate(this % pred)
        end if
    end subroutine density_deallocate

    subroutine density_destructor(this)
        type(density_t), intent(inout) :: this
        call this % deallocate()
    end subroutine density_destructor

end module slbm_2d_density
