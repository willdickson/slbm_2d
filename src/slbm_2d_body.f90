module slbm_2d_body

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_const,    only : PI
    use slbm_2d_const,    only : VECTOR_ZERO
    use slbm_2d_const,    only : BODY_TYPE_OPEN
    use slbm_2d_const,    only : BODY_TYPE_CLOSED
    use slbm_2d_const,    only : BODY_TYPE_UNKNOWN
    use slbm_2d_vector,   only : vector_t
    use slbm_2d_nbrs,     only : nbrs_t
    use slbm_2d_state,    only : state_t
    use slbm_2d_mesh,     only : mesh_t
    use slbm_2d_line_seg, only : line_seg_t
    use slbm_2d_line_seg, only : intersect
    use slbm_2d_line_seg, only : is_chain

    implicit none
    private


    type, public :: body_t
        integer(ip)                 :: type_id = BODY_TYPE_UNKNOWN
        !type(motion_t)              :: motion    !
        type(vector_t), allocatable :: pos(:)    ! positions of body points 
        type(vector_t), allocatable :: vel(:)    ! velocities of body point
        type(nbrs_t),   allocatable :: nbrs(:)   ! neighboring mesh pos
        real(wp),       allocatable :: rho(:)    ! density at body points
    contains
        private
        procedure, public  :: update
        procedure          :: update_pos_and_vel
        procedure          :: update_nbrs_and_rho
        procedure, public  :: num_pos
        procedure          :: check_pos
    end type body_t


    interface body_t
        procedure :: body_constructor
    end interface body_t


contains


    function body_constructor(type_id, pos) result(body)
        integer(ip),    intent(in) :: type_id
        type(vector_t), intent(in) :: pos(:)
        type(body_t)               :: body
        integer(ip)                :: i

        body % type_id = type_id
        body % pos  = pos
        allocate(body % vel(size(pos)), source=VECTOR_ZERO)
        allocate(body % rho(size(pos)), source=0.0_wp)
        allocate(body % nbrs(size(pos)))
        do i = 1, size(pos)
            call body % nbrs(i) % set_to_zero()
        end do
        if (.not. body % check_pos()) then
            print *, 'body is not a simple curve'
            stop
        end if
    end function body_constructor

    
    subroutine update(this, state, mesh, ds, time)
        class(body_t), intent(inout) :: this    ! the current body 
        type(state_t), intent(in)    :: state   ! fluid state: ux, uy, rho
        type(mesh_t),  intent(in)    :: mesh    ! x and y meshes
        real(wp),      intent(in)    :: ds      ! mesh spacing
        real(wp),      intent(in)    :: time    ! simulation time 
        call update_pos_and_vel(this, state, mesh, ds, time)
        call update_nbrs_and_rho(this, state, mesh, ds, time)
    end subroutine update


    subroutine update_pos_and_vel(this, state, mesh, ds, time)
        class(body_t), intent(inout) :: this    ! the current body 
        type(state_t), intent(in)    :: state   ! fluid state: ux, uy, rho
        type(mesh_t),  intent(in)    :: mesh    ! x and y meshes
        real(wp),      intent(in)    :: ds      ! mesh spacing
        real(wp),      intent(in)    :: time    ! simulation time 
    end subroutine update_pos_and_vel


    subroutine update_nbrs_and_rho(this, state, mesh, ds, time)
        class(body_t), intent(inout) :: this    ! the current body 
        type(state_t), intent(in)    :: state   ! fluid state: ux, uy, rho
        type(mesh_t),  intent(in)    :: mesh    ! x and y meshes
        real(wp),      intent(in)    :: ds      ! mesh spacing
        real(wp),      intent(in)    :: time    ! simulation time 

        real(wp)                     :: r       ! radial dist. from pos to mesh pt
        real(wp)                     :: big_r   ! a big radial dist. start value 
        real(wp)                     :: min_r   ! current minimal radial dist. 
        real(wp)                     :: xm      ! mesh point x
        real(wp)                     :: ym      ! mesh point y
        real(wp)                     :: xb      ! body point x
        real(wp)                     :: yb      ! body point y
        integer(ip)                  :: num_x   ! mesh size in x-dim
        integer(ip)                  :: num_y   ! mesh size in y-dim
        integer(ip)                  :: nnbrs   ! body point # of neighbors
        integer(ip)                  :: i_min   ! min i index for search
        integer(ip)                  :: i_max   ! max i index for search
        integer(ip)                  :: j_min   ! min j index for search
        integer(ip)                  :: j_max   ! max j index for search
        integer(ip)                  :: i, j, k ! looping indices

        num_x = size(state % rho, 1)
        num_y = size(state % rho, 2)
        big_r = ds*sqrt((num_x-1.0_wp)**2 + (num_y - 1.0_wp)**2)

        ! For each body pint find mesh point neighbors closest mesh cell 
        do k = 1, size(this % pos)

            call this % nbrs(k) % set_to_zero()

            ! Get search indices for neighbors which are inside mesh
            i_min = ceiling((this % pos(k) % x - 2.0_wp*ds)/ds) + 1_wp
            j_min = ceiling((this % pos(k) % y - 2.0_wp*ds)/ds) + 1_wp
            i_max = floor((this % pos(k) % x + 2.0_wp*ds)/ds) + 1_wp
            j_max = floor((this % pos(k) % y + 2.0_wp*ds)/ds) + 1_wp
            i_min = min(max(i_min,1), num_x)
            j_min = min(max(j_min,1), num_y)
            i_max = min(max(i_max,1), num_x)
            j_max = min(max(j_max,1), num_y)

            min_r = big_r
            do i = i_min, i_max
                do j = j_min, j_max
                    ! Position of mesh point
                    xm = mesh % x(i,j)
                    ym = mesh % y(i,j)

                    ! Add neighbors to array
                    nnbrs = this % nbrs(k) % num + 1
                    this % nbrs(k) % num  = nnbrs
                    this % nbrs(k) % ix(nnbrs) = i
                    this % nbrs(k) % jy(nnbrs) = j
                    this % nbrs(k) % pos(nnbrs) = vector_t(xm, ym)
                    this % nbrs(k) % u(nnbrs) = state % u(i,j) 

                    ! Find closest cell to body point
                    xb = this % pos(k) % x
                    yb = this % pos(k) % y
                    r = sqrt((xm-xb)**2 + (ym-yb)**2)
                    if (r < min_r) then
                        min_r = r
                        this % rho(k) = state % rho(i,j) 
                    end if
                end do
            end do 
        end do
    end subroutine update_nbrs_and_rho


    elemental function num_pos(this) result(num)
        class(body_t), intent(in) :: this
        integer(ip)               :: num
        if (.not. allocated(this % pos)) then
            num = 0_ip
        else
            num = size(this % pos)
        end if
    end function num_pos


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
