module slbm_2d_ibsol

    use slbm_2d_kinds, only : wp, ip
    use slbm_2d_body,  only : body_t
    use slbm_2d_spmat, only : spmat_t
    implicit none
    private

    type, public :: ibsol_t
        type(body_t), allocatable :: body(:)   ! Array of immersed bodies
        type(spmat_t)             :: a         ! A matrix for velocity correction 
        real(wp), allocatable     :: bx(:)     ! b vector for velocity corr. x comp.
        real(wp), allocatable     :: by(:)     ! b vector for velocity corr. y comp. 
        real(wp), allocatable     :: vx(:)     ! velocity correction x component 
        real(wp), allocatable     :: vy(:)     ! velocity correction y component
    contains
        private
        procedure, public :: num_body => ibsol_num_body
    end type ibsol_t

    interface ibsol_t
        procedure :: ibsol_constructor
    end interface ibsol_t

contains

    function ibsol_constructor(body) result(ibsol)
        type(body_t), intent(in) :: body(:)
        type(ibsol_t)            :: ibsol
        integer(ip)              :: num_pos
        integer(ip)              :: i
        allocate(ibsol % body(size(body)))
        ibsol % body = body
        num_pos = 0_ip
        do i = 1, size(ibsol % body) 
            num_pos = num_pos + ibsol % body(i) % num_pos()
        end do
        ibsol % a = spmat_t(num_pos)
        allocate(ibsol % bx(num_pos), source=0.0_wp)
        allocate(ibsol % by(num_pos), source=0.0_wp)
        allocate(ibsol % vx(num_pos), source=0.0_wp)
        allocate(ibsol % vy(num_pos), source=0.0_wp)
    end function ibsol_constructor


    function ibsol_num_body(this) result(num)
        class(ibsol_t), intent(in) :: this
        integer(ip)                :: num
        num = size(this % body)
    end function ibsol_num_body

end module slbm_2d_ibsol
