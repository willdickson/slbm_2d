module slbm_2d_celltype
    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    integer(ip), parameter :: CELLTYPE_FLUID = 0_ip
    integer(ip), parameter :: CELLTYPE_BNDRY = 1_ip
    integer(ip), parameter :: CELLTYPE_EMPTY = 2_ip

    type, public :: celltype_t
        integer(ip)              :: num_x
        integer(ip)              :: num_y
        integer(ip), allocatable :: id(:,:)
    end type celltype_t

    interface celltype
        procedure celltype_constructor
    end interface celltype

contains

    function celltype_constructor(num_x, num_y) result(celltype)
        integer(ip), intent(in) :: num_x
        integer(ip), intent(in) :: num_y
        type(celltype_t)        :: celltype
        allocate( celltype % id(num_x, num_y) )
        celltype % id = CELLTYPE_FLUID 
    end function celltype_constructor

end module slbm_2d_celltype
