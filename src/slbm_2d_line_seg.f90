module slbm_2d_line_seg
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_vector, only : vector_t
    implicit none
    private

    type, public :: line_seg_t
        type(vector_t) :: p
        type(vector_t) :: q
    end type line_seg_t

    public :: intersect



contains

    function intersect(seg1, seg2) result(res)
        type(line_seg_t), intent(in) :: seg1
        type(line_seg_t), intent(in) :: seg2
        logical                      :: test1
        logical                      :: test2
        logical                      :: res
        test1 = ccw(seg1 % p, seg2 % p, seg2 % q) .neqv. ccw(seg1 % q, seg2 % p, seg2 % q) 
        test2 = ccw(seg1 % p, seg1 % q, seg2 % p) .neqv. ccw(seg1 % p, seg1 % q, seg2 % q) 
        res = test1 .and. test2
    end function intersect


    function ccw(a,b,c) result(res)
        type(vector_t), intent(in) :: a
        type(vector_t), intent(in) :: b
        type(vector_t), intent(in) :: c
        logical                    :: res
        res = (c % y - a % y) * (b % x - a % x) > (b % y-a % y) * (c % x - a % x)
    end function ccw

end module slbm_2d_line_seg
