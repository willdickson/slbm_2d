module slbm_2d_nbrs
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_const,  only : VECTOR_ZERO
    use slbm_2d_const,  only : BODY_PT_MAX_NBRS
    use slbm_2d_vector, only : vector_t

    implicit none
    private

    type, public :: nbrs_t
        integer(ip)    :: num = 0_ip            ! Number of mesh neighbors
        integer(ip)    :: ix(BODY_PT_MAX_NBRS)  ! x indices of mesh neighbors
        integer(ip)    :: jy(BODY_PT_MAX_NBRS)  ! y indices of mesh neighbors
        type(vector_t) :: pos(BODY_PT_MAX_NBRS) ! position of mesh neighbors
        type(vector_t) :: u(BODY_PT_MAX_NBRS)   ! fluid velocity of mesh neighbors
    contains
        private
        procedure, public :: set_to_zero
    end type

contains

    subroutine set_to_zero(this)
        class(nbrs_t), intent(inout) :: this
        this % num = 0_ip
        this % ix  = 0_ip
        this % jy  = 0_ip
        this % pos = VECTOR_ZERO 
        this % u   = VECTOR_ZERO 
    end subroutine set_to_zero

end module slbm_2d_nbrs
