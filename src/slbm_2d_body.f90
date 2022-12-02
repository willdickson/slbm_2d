module slbm_2d_body
    use slbm_2d_kinds,   only : wp, ip
    use slbm_2d_const,   only : BODY_TYPE_OPEN
    use slbm_2d_const,   only : BODY_TYPE_CLOSED
    use slbm_2d_const,   only : BODY_TYPE_UNKNOWN
    use slbm_2d_vector,  only : vector_t
    implicit none
    private

    type, public :: body_t
        integer(ip)                 :: type_id = BODY_TYPE_UNKNOWN
        type(vector_t), allocatable :: points(:)
        logical                     :: moving = .false.
    contains
    end type body_t

    interface body_t
        procedure :: body_constructor
    end interface body_t

contains

    function body_constructor(type_id, pts, moving) result(body)
        integer(ip), intent(in)    :: type_id
        type(vector_t), intent(in) :: points(:)
        logical, intent(in)        :: moving

        body % type_id = type_id
        body % points  = points
        body % moving  = moving

        type(body_t) :: body
    end function body_constructor

end module slbm_2d_body
