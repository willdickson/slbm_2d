module slbm_2d_spmat
    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    type, public :: spmat_t
        integer(ip)              :: nnz = 0_ip
        integer(ip), allocatable :: ix(:)
        integer(ip), allocatable :: iy(:)
        real(wp),    allocatable :: val(:) 
    end type spmat_t

    interface spmat_t
        procedure :: spmat_constructor
    end interface spmat_t

contains

    function spmat_constructor(n) result(spmat)
        integer, intent(in) :: n 
        type(spmat_t)       :: spmat
        spmat % nnz = 0_ip
        allocate(spmat % ix(n),  source=0_ip)
        allocate(spmat % iy(n),  source=0_ip)
        allocate(spmat % val(n), source=0.0_wp)
    end function spmat_constructor

end module slbm_2d_spmat
