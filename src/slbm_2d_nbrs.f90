module slbm_2d_nbrs
    use slbm_2d_kinds, only : wp, ip
    use slbm_2d_const, only : BODY_PT_MAX_NBRS

    implicit none
    private

    type, public :: nbrs_t
        integer(ip) :: num = 0_ip            ! Number of mesh neighbors
        integer(ip) :: ix(BODY_PT_MAX_NBRS)  ! x indices of mesh neighbors
        integer(ip) :: iy(BODY_PT_MAX_NBRS)  ! y indices of mesh neihbors
    end type

end module slbm_2d_nbrs
