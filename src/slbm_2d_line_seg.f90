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

    integer(ip), parameter :: CLOCKWISE      =  1_ip
    integer(ip), parameter :: CNTR_CLOCKWISE = -1_ip
    integer(ip), parameter :: COLLINEAR      =  0_ip


contains

    function intersect(seg1, seg2) result(test)
        type(line_seg_t), intent(in) :: seg1
        type(line_seg_t), intent(in) :: seg2
        logical                      :: test
        test = .false.

    end function intersect

    function orientation(p, q, r) result(ori_val)
        type(vector_t), intent(in) :: p
        type(vector_t), intent(in) :: q
        type(vector_t), intent(in) :: r
        real(wp)                   :: ori_val
        integer(ip)                :: ori_int
        ori_val = (q % y - p % y) * (r % x - q % x) - (q % x - p % x) * (r % y - q % y)
        if (ori_val > 0) then
            ori_int = CLOCKWISE 
        else if (ori_val < 0) then 
            ori_int = CNTR_CLOCKWISE 
        else
            ori_int = COLLINEAR 
        end if
    end function orientation

end module slbm_2d_line_seg
