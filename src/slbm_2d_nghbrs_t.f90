module slbm_2d_nghbrs
    use slbm_2d_kinds, only : wp, ip
    use slbm_2d_const, only : BODY_PT_MAX_NGHBRS

    implicit none
    private

    type, public :: nghbrs_t
        integer(ip) :: num = 0_ip              ! Number of mesh neighbors
        integer(ip) :: ix(BODY_PT_MAX_NGHBRS)  ! x indices of mesh neighbors
        integer(ip) :: iy(BODY_PT_MAX_NGHBRS)  ! y indices of mesh neihbors
    end type

end module slbm_2d_nghbrs
