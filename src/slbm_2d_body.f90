module slbm_2d_body
    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    type, public :: body_t
        integer(ip)           :: type_id
        real(wp), allocatable :: x(:)
        real(ip), allocatable :: y(:)
    contains
    end type body_t

    interface body_t
        procedure :: body_constructor
    end interface body_t

contains

    function body_constructor() result(body)
        type(body_t) :: body
    end function body_constructor

end module slbm_2d_body
