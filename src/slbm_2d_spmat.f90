module slbm_2d_spmat
    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    type, public :: spmat_t
        integer(ip)              :: n = 0_ip    ! Matrix size  (nxn)
        integer(ip)              :: nnz = 0_ip  ! Number nonzero elements
        integer(ip), allocatable :: ix(:)       ! Nonzero element row indices
        integer(ip), allocatable :: jy(:)       ! Nonzero element col indices
        real(wp),    allocatable :: val(:)      ! Nonzero element values
    contains
        private
        procedure, public :: density      => spmat_density
        procedure, public :: zero         => spmat_zero
        procedure, public :: as_dense_mat => spmat_as_dense_mat
    end type spmat_t

    interface spmat_t
        procedure :: spmat_constructor
    end interface spmat_t

contains

    function spmat_constructor(n) result(spmat)
        integer, intent(in) :: n 
        type(spmat_t)       :: spmat
        spmat % n = n 
        spmat % nnz = 0_ip
        allocate(spmat % ix(n**2),  source=0_ip)
        allocate(spmat % jy(n**2),  source=0_ip)
        allocate(spmat % val(n**2), source=0.0_wp)
    end function spmat_constructor


    function spmat_density(this) result(density)
        class(spmat_t), intent(in) :: this
        real(wp)                   :: density
        if (this % n == 0_ip) then
            density = 0.0_wp
        else
            density = real(this % nnz, wp)/real((this % n)**2, wp)
        end if
    end function spmat_density


    subroutine spmat_zero(this)
        class(spmat_t), intent(inout) :: this
        this % nnz = 0_ip
        this % ix  = 0_ip
        this % jy  = 0_ip
        this % val = 0.0_wp
    end subroutine spmat_zero


    function spmat_as_dense_mat(this) result(a)
        class(spmat_t), intent(in) :: this
        real(wp), allocatable      :: a(:,:)
        integer(ip)                :: i, j, k
        allocate(a(this % n, this % n), source=0.0_wp)
        do k = 1, this % nnz
            i = this % ix(k)
            j = this % jy(k)
            a(i,j) = this % val(k)
        end do
    end function spmat_as_dense_mat

end module slbm_2d_spmat
