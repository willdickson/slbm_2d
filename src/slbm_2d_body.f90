module slbm_2d_body

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_const,    only : BODY_TYPE_OPEN
    use slbm_2d_const,    only : BODY_TYPE_CLOSED
    use slbm_2d_const,    only : BODY_TYPE_UNKNOWN
    use slbm_2d_vector,   only : vector_t
    use slbm_2d_nbrs,     only : nbrs_t
    use slbm_2d_line_seg, only : line_seg_t
    use slbm_2d_line_seg, only : intersect
    use slbm_2d_line_seg, only : is_chain

    implicit none
    private


    type, public :: body_t
        integer(ip)                 :: type_id = BODY_TYPE_UNKNOWN
        type(vector_t), allocatable :: pos(:)  ! positions of body points 
        type(vector_t), allocatable :: vel(:)  ! velocities of body point
        type(nbrs_t), allocatable   :: nbrs(:) ! neighboring mesh pos
    contains
        private
        procedure, public  :: num_pos
        procedure, public  :: update_nbrs
        procedure, public  :: check_pos
    end type body_t


    interface body_t
        procedure :: body_constructor
    end interface body_t


contains


    function body_constructor(type_id, pos) result(body)
        integer(ip), intent(in)    :: type_id
        type(vector_t), intent(in) :: pos(:)
        type(body_t)               :: body

        body % type_id = type_id
        body % pos  = pos

        allocate(body % vel(size(pos)))
        body % vel = vector_t(0.0_wp, 0.0_wp)

        allocate(body % nbrs(size(pos)))
        body % nbrs = nbrs_t(0_ip, 0_ip, 0_ip)

        if (.not. body % check_pos()) then
            print *, 'body is not a simple curve'
            stop
        end if
    end function body_constructor


    elemental function num_pos(this) result(num)
        class(body_t), intent(in) :: this
        integer(ip)               :: num
        if (.not. allocated(this % pos)) then
            num = 0_ip
        else
            num = size(this % pos)
        end if
    end function num_pos


    subroutine update_nbrs(this, num_x, num_y, ds)
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

        this % nbrs  = nbrs_t(0_ip, 0_ip, 0_ip)
        do k = 1, size(this % pos)

            i_min = ceiling((this % pos(k) % x - 2.0_wp*ds)/ds)
            j_min = ceiling((this % pos(k) % y - 2.0_wp*ds)/ds)
            i_max = floor((this % pos(k) % x + 2.0_wp*ds)/ds)
            j_max = floor((this % pos(k) % y + 2.0_wp*ds)/ds)

            i_min = min(max(i_min,1), num_x)
            j_min = min(max(j_min,1), num_y)
            i_max = min(max(i_max,1), num_x)
            j_max = min(max(j_max,1), num_y)

            do i = i_min, i_max
                do j = j_min, j_max
                    num = this % nbrs(k) % num + 1
                    this % nbrs(k) % num  = num
                    this % nbrs(k) % ix(num) = i
                    this % nbrs(k) % iy(num) = j
                end do
            end do 
        end do


    end subroutine update_nbrs


    elemental function check_pos(this) result(ok)
        class(body_t), intent(in) :: this
        type(line_seg_t)          :: seg1
        type(line_seg_t)          :: seg2
        integer(ip)               :: i
        integer(ip)               :: j
        logical                   :: ok

        ! Check positions to make sure that body is simple curve.
        ok = .true.
        do i = 1, size(this % pos)-2
            do j = i+1, size(this % pos)-1
                seg1 = line_seg_t(this % pos(i), this % pos(i+1))
                seg2 = line_seg_t(this % pos(j), this % pos(j+1))
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
    end function check_pos


end module slbm_2d_body
