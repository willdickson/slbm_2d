module slbm_2d_body
    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_const,    only : BODY_TYPE_OPEN
    use slbm_2d_const,    only : BODY_TYPE_CLOSED
    use slbm_2d_const,    only : BODY_TYPE_UNKNOWN
    use slbm_2d_vector,   only : vector_t
    use slbm_2d_nghbrs,   only : nghbrs_t
    use slbm_2d_line_seg, only : line_seg_t
    use slbm_2d_line_seg, only : intersect
    use slbm_2d_line_seg, only : is_chain
    implicit none
    private

    type, public :: body_t
        integer(ip)                 :: type_id = BODY_TYPE_UNKNOWN
        type(vector_t), allocatable :: points(:)  ! points defining body
        type(nghbrs_t), allocatable :: nghbrs(:)  ! neighboring mesh points
    contains
        private
        procedure, public  :: num_pts
        procedure, public  :: update_nghbrs
        procedure, public  :: check_points
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
        allocate(body % nghbrs(size(points)))
        body % nghbrs = nghbrs_t(0_ip, 0_ip, 0_ip)
        call body % check_points()
    end function body_constructor


    elemental function num_pts(this) result(num)
        class(body_t), intent(in) :: this
        integer(ip)               :: num
        if (.not. allocated(this % points)) then
            num = 0_ip
        else
            num = size(this % points)
        end if
    end function num_pts


    subroutine update_nghbrs(this, num_x, num_y, ds)
        class(body_t), intent(inout) :: this
        integer(ip), intent(in)      :: num_x
        integer(ip), intent(in)      :: num_y
        real(wp), intent(in)         :: ds

        integer(ip)                  :: num
        integer(ip)                  :: i_min
        integer(ip)                  :: i_max
        integer(ip)                  :: j_min
        integer(ip)                  :: j_max
        integer(ip)                  :: i, j, k

        this % nghbrs  = nghbrs_t(0_ip, 0_ip, 0_ip)
        do k = 1, size(this % points)

            i_min = ceiling((this % points(k) % x - 2.0_wp*ds)/ds)
            j_min = ceiling((this % points(k) % y - 2.0_wp*ds)/ds)
            i_max = floor((this % points(k) % x + 2.0_wp*ds)/ds)
            j_max = floor((this % points(k) % y + 2.0_wp*ds)/ds)

            i_min = min(max(i_min,1), num_x)
            j_min = min(max(j_min,1), num_y)
            i_max = min(max(i_max,1), num_x)
            j_max = min(max(j_max,1), num_y)

            do i = i_min, i_max
                do j = j_min, j_max
                    num = this % nghbrs(k) % num + 1
                    this % nghbrs(k) % num  = num
                    this % nghbrs(k) % ix(num) = i
                    this % nghbrs(k) % iy(num) = j
                end do
            end do 
        end do

    end subroutine update_nghbrs


    subroutine check_points(this)
        class(body_t), intent(in) :: this
        type(line_seg_t)          :: seg1
        type(line_seg_t)          :: seg2
        integer(ip)               :: i
        integer(ip)               :: j
        logical                   :: th_test
        logical                   :: on_test
        logical                   :: ok

        ! Make sure that body is simple curve.
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

        ! May want to check spacing??

    end subroutine check_points

end module slbm_2d_body
