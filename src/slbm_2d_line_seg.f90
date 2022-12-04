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
        procedure, public :: tail_equal_head
        procedure, public :: head_equal_tail
        procedure, public :: on_segment
        procedure, public :: equal
        procedure, public :: not_equal
        generic, public   :: operator(==) => equal 
        generic, public   :: operator(/=) => not_equal
    end type line_seg_t

    public :: collinear
    public :: intersect 
    public :: orientation
    public :: is_chain

    interface collinear
        procedure :: collinear_seg_and_pt
        procedure :: collinear_pt_and_pt
    end interface collinear

    integer(ip), parameter :: ORIENT_CL  = 0 ! collinear
    integer(ip), parameter :: ORIENT_CW  = 1 ! clockwise
    integer(ip), parameter :: ORIENT_CCW = 2 ! counter clockwise

contains

    elemental function reverse(this) result(seg)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t)              :: seg
        seg = line_seg_t(this % q, this % p)
    end function reverse


    elemental function tail_equal_head(this, other) result(res)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t), intent(in)  :: other
        logical                       :: res
        res = (this % q == other % p)
    end function tail_equal_head


    elemental function on_segment(this, r) result(res)
        class(line_seg_t), intent(in) :: this
        type(vector_t), intent(in)   :: r
        logical                      :: res
        logical                      :: test1
        logical                      :: test2
        logical                      :: test3
        logical                      :: test4
        if (collinear(this,r)) then
            test1 = r % x <= max(this % p % x, this % q % x)
            test2 = r % x >= min(this % p % x, this % q % x)
            test3 = r % y <= max(this % p % y, this % q % y)
            test4 = r % y >= min(this % p % y, this % q % y)
            res = test1 .and. test2 .and. test3 .and. test4
        else
            res = .false.
        end if
    end function on_segment


    elemental function head_equal_tail(this, other) result(res)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t), intent(in)  :: other
        logical                       :: res
        res = (this % p == other % q)
    end function head_equal_tail


    elemental function equal(this, other) result(res)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t), intent(in)  :: other
        logical                       :: res
        res = (this % p == other % p) .and. (this % q == other % q)
    end function equal


    elemental function not_equal(this, other) result(res)
        class(line_seg_t), intent(in) :: this
        type(line_seg_t), intent(in)  :: other
        logical                       :: res
        res = .not. (this == other)
    end function not_equal


    elemental function intersect(seg1, seg2) result(res)
        type(line_seg_t), intent(in) :: seg1
        type(line_seg_t), intent(in) :: seg2
        integer(ip)                  :: ori1
        integer(ip)                  :: ori2
        integer(ip)                  :: ori3
        integer(ip)                  :: ori4
        logical                      :: res 

        res = .false.
        ori1 = orientation(seg1 % p, seg1 % q, seg2 % p) 
        ori2 = orientation(seg1 % p, seg1 % q, seg2 % q) 
        ori3 = orientation(seg2 % p, seg2 % q, seg1 % p) 
        ori4 = orientation(seg2 % p, seg2 % q, seg1 % q) 
        if ((ori1 /= ori2) .and. (ori3 /= ori4)) then
            res = .true.
        else
            if ((ori1 == 0) .and.  seg1 % on_segment(seg2 % p))then
                if (seg2 % p /= seg1 % q) then
                    res = .true.
                end if
            else if ((ori2 == 0) .and. seg1 % on_segment(seg2 % q))  then
                res = .true.
            else if ((ori3 == 0) .and. seg2 % on_segment(seg1 % p))  then
                res = .true.
            else if ((ori4 == 0) .and. seg2 % on_segment(seg1 % q))  then
                res = .true.
            end if
        end if
    end function intersect



    elemental function collinear_seg_and_pt(seg, r) result(res)
        type(line_seg_t), intent(in) :: seg
        type(vector_t), intent(in)   :: r
        real(wp)                     :: val
        logical                      :: res
        select case(orientation(seg % p, seg %q, r))
        case (ORIENT_CL)
            res = .true.
        case default
            res = .false.
        end select
    end function collinear_seg_and_pt


    elemental function collinear_pt_and_pt(a,b,c) result(res)
        type(vector_t), intent(in) :: a
        type(vector_t), intent(in) :: b
        type(vector_t), intent(in) :: c
        logical                :: res
        select case(orientation(a,b,c))
        case (ORIENT_CL)
            res = .true.
        case default
            res = .false.
        end select
    end function collinear_pt_and_pt


    elemental function orientation(a,b,c) result(res)
        type(vector_t), intent(in) :: a
        type(vector_t), intent(in) :: b
        type(vector_t), intent(in) :: c
        real(wp)                   :: val
        integer(ip)                :: res
        val = (b % y - a % y) * (c % x - b % x) - (b % x - a % x) * (c % y - b % y) 
        if (val == 0.0_wp) then
            res = ORIENT_CL 
        else if (val > 0.0_wp) then
            res = ORIENT_CW
        else
            res = ORIENT_CCW 
        end if
    end function orientation 

    elemental function is_chain(seg1,seg2) result(res)
        type(line_seg_t), intent(in) :: seg1
        type(line_seg_t), intent(in) :: seg2
        logical                      :: res
        logical                      :: test1
        logical                      :: test2
        test1 = seg1 % tail_equal_head(seg2) 
        test2 = seg1 % on_segment(seg2 % q)
        res = test1 .and. (.not. test2)
    end function is_chain


end module slbm_2d_line_seg
