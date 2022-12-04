module slbm_2d_body
    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_const,    only : BODY_TYPE_OPEN
    use slbm_2d_const,    only : BODY_TYPE_CLOSED
    use slbm_2d_const,    only : BODY_TYPE_UNKNOWN
    use slbm_2d_vector,   only : vector_t
    use slbm_2d_line_seg, only : line_seg_t
    use slbm_2d_line_seg, only : intersect
    use slbm_2d_line_seg, only : is_chain
    implicit none
    private

    type, public :: body_t
        integer(ip)                 :: type_id = BODY_TYPE_UNKNOWN
        type(vector_t), allocatable :: points(:)
    contains
        private
        procedure  :: check_points
    end type body_t

    interface body_t
        procedure :: body_constructor
    end interface body_t

contains

    function body_constructor(type_id, points) result(body)
        integer(ip), intent(in)    :: type_id
        type(vector_t), intent(in) :: points(:)
        type(body_t)               :: body
        body % type_id = type_id
        body % points  = points
    end function body_constructor


    subroutine check_points(this)
        class(body_t), intent(in) :: this
        type(line_seg_t)          :: seg1
        type(line_seg_t)          :: seg2
        integer(ip)               :: i
        integer(ip)               :: j
        logical                   :: th_test
        logical                   :: on_test
        logical                   :: ok
        ok = .true.
        do i = 1, size(this % points)-2
            do j = i+1, size(this % points)-1
                seg1 = line_seg_t(this % points(i), this % points(i+1))
                seg2 = line_seg_t(this % points(j), this % points(j+1))
                if (intersect(seg1, seg2)) then
                    if (j /= (i+1)) then
                        ok = .false.
                    else
                        if (.not. is_chain(seg1, seg2)) then
                            ok = .false.
                        end if
                    end if
                end if
            end do
        end do
    end subroutine check_points

end module slbm_2d_body
