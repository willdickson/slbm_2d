module slbm_2d_line_seg
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_vector, only : vector_t
    implicit none
    private

    type, public :: line_seg_t
        type(vector_t) :: p
        type(vector_t) :: q
    contains
        private
        procedure, public :: reverse
        procedure, public :: equals
        generic, public   :: operator(==) => equals 
    end type line_seg_t

    public :: intersect 
    public :: ccw



contains

    function reverse(this) result(seg)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t)              :: seg
        seg = line_seg_t(this % q, this % p)
    end function reverse


    function equals(this, other) result(res)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t), intent(in)  :: other
        logical                       :: res
        res = (this % p == other % p) .and. (this % q == other % q)
    end function equals


    function intersect(seg1, seg2) result(res)
        type(line_seg_t), intent(in) :: seg1
        type(line_seg_t), intent(in) :: seg2
        logical                      :: test1
        logical                      :: test2
        logical                      :: test3
        logical                      :: test4
        logical                      :: test5
        logical                      :: res
        ! Fails if two segments are identical
        test1 = ccw(seg1 % p, seg2 % p, seg2 % q) 
        test2 = ccw(seg1 % q, seg2 % p, seg2 % q) 
        test3 = ccw(seg1 % p, seg1 % q, seg2 % p) 
        test4 = ccw(seg1 % p, seg1 % q, seg2 % q) 
        test5 = (seg1 == seg2) .or. (seg1 == seg2 % reverse())
        res = ((test2 .neqv. test1) .and. (test3 .neqv. test4)) .or. test5
    end function intersect


    function ccw(a,b,c) result(res)
        type(vector_t), intent(in) :: a
        type(vector_t), intent(in) :: b
        type(vector_t), intent(in) :: c
        real(wp)                   :: val1
        real(wp)                   :: val2
        logical                    :: res
        val1 = (c % y - a % y) * (b % x - a % x) 
        val2 = (b % y - a % y) * (c % x - a % x)
        res = val1 > val2
    end function ccw

end module slbm_2d_line_seg
