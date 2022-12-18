module slbm_2d_ibsol

    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_const,  only : VECTOR_ZERO
    use slbm_2d_const,  only : BODY_TYPE_OPEN
    use slbm_2d_vector, only : vector_t
    use slbm_2d_body,   only : body_t
    use slbm_2d_mesh,   only : mesh_t
    use slbm_2d_spmat,  only : spmat_t
    use slbm_2d_state,  only : state_t
    use slbm_2d_funcs,  only : kernel
    use slbm_2d_funcs,  only : constrain
    use minres,         only : minres_ez_t
    use minres,         only : minres_info_t

    implicit none
    private

    type, public :: ibsol_t
        type(body_t), allocatable :: body(:)     ! Array of immersed bodies
        integer(ip)               :: num_pos     ! total number of pos pts (all bodies)
        type(spmat_t)             :: a           ! A matrix for velocity correction 
        real(wp), allocatable     :: bx(:)       ! b vector for velocity corr. x comp.
        real(wp), allocatable     :: by(:)       ! b vector for velocity corr. y comp. 
        real(wp), allocatable     :: vx(:)       ! velocity correction x component 
        real(wp), allocatable     :: vy(:)       ! velocity correction y component
        integer(ip), allocatable  :: id(:,:)     ! marker ids markers for mesh points 
                                                 !   id = 0         for fluid
                                                 !   id = 1, 2, ... for body interiors
    contains
        private
        procedure, public :: update    => ibsol_update
        procedure, public :: update_id => ibsol_update_id
        procedure, public :: corrector => ibsol_corrector
        procedure, public :: num_body  => ibsol_num_body
    end type ibsol_t

    interface ibsol_t
        procedure :: ibsol_constructor
    end interface ibsol_t

contains

    function ibsol_constructor(body, num_x, num_y) result(ibsol)
        type(body_t), intent(in) :: body(:)
        integer(ip), intent(in)  :: num_x
        integer(ip), intent(in)  :: num_y
        type(ibsol_t)            :: ibsol
        integer(ip)              :: num_pos
        integer(ip)              :: i
        allocate(ibsol % body(size(body)))
        ibsol % body = body
        num_pos = 0_ip
        do i = 1, size(ibsol % body) 
            num_pos = num_pos + ibsol % body(i) % num_pos()
        end do
        ibsol % num_pos = num_pos
        ibsol % a = spmat_t(num_pos)
        allocate(ibsol % bx(num_pos), source=0.0_wp)
        allocate(ibsol % by(num_pos), source=0.0_wp)
        allocate(ibsol % vx(num_pos), source=0.0_wp)
        allocate(ibsol % vy(num_pos), source=0.0_wp)
        allocate(ibsol % id(num_x, num_y), source=0_ip)
    end function ibsol_constructor


    subroutine ibsol_update(this, state, mesh, ds, time)
        class(ibsol_t), intent(inout) :: this    ! immersed boundry solver
        type(state_t), intent(inout)  :: state   ! fluid state: ux, uy, rho
        type(mesh_t), intent(in)      :: mesh    ! x and y meshes
        real(wp), intent(in)          :: ds      ! mesh spacing
        real(wp), intent(in)          :: time    ! simulation time 
        integer(ip)                   :: i
        do i = 1, this % num_body()
            call this % body(i) % update_pos_and_vel(state, mesh, ds, time)
        end do
        call this % update_id(state, mesh, ds)
        do i = 1, this % num_body()
            call this % body(i) % update_nbrs_and_rho(state, mesh, this % id,  ds, time)
        end do
    end subroutine ibsol_update


    subroutine ibsol_update_id(this, state, mesh, ds)
        class(ibsol_t), intent(inout), target :: this    ! immersed boundry solver
        type(state_t),  intent(inout)         :: state   ! fluid state: ux, uy, rho
        type(mesh_t),   intent(in)            :: mesh    ! x and y meshes
        real(wp),       intent(in)            :: ds      ! mesh spacing

        type(body_t), pointer  :: body ! alias for current body
        type(vector_t)         :: p
        integer(ip)            :: imin ! minimum i search index 
        integer(ip)            :: imax ! maximum i search index
        integer(ip)            :: jmin ! minimum j search index
        integer(ip)            :: jmax ! maximum j search index
        integer(ip)            :: n    ! index for looping over bodies
        integer(ip)            :: i    ! index for looping over x grid cells
        integer(ip)            :: j    ! index for looping over y grid cells

        ! Reset all ids
        this % id = 0
        ! Loop over all immersed bodies
        do n = 1, this % num_body()
            body => this % body(n)
            ! Skip open bodies - no interior
            if (body % type_id == BODY_TYPE_OPEN) then
                cycle
            end if
            ! Get range of i,j indices for to search
            call bounding_ind(body, mesh, ds, imin, imax, jmin, jmax)
            do i = imin, imax
                do j = jmin, jmax
                    ! check if point is inside body
                    p = vector_t(mesh % x(i,j), mesh % y(i,j))
                    if (body % is_interior(p) ) then
                        this % id(i,j) = n
                        state % u(i,j) = VECTOR_ZERO
                    end if
                end do
            end do
        end do
    end subroutine ibsol_update_id


    subroutine ibsol_corrector(this, ds, u)
        class(ibsol_t), intent(inout), target :: this    ! immersed boundry solver
        real(wp),       intent(in)            :: ds      ! mesh spacing
        type(vector_t), intent(inout)         :: u(:,:)  ! velocity mesh to be corrected


        type(body_t), pointer   :: body         ! alias for current body
        type(minres_ez_t)       :: minres_ez    ! sparse linear solver
        type(minres_info_t)     :: minres_info  ! linear solver info
        real(wp)                :: kerni        ! kernel value for ith pos
        real(wp)                :: kernj        ! kernel value for jth pos
        real(wp)                :: aij          ! A matrix value at i,j
        real(wp)                :: ds2          ! square of mesh spacing ds**2
        integer(ip)             :: cnt          ! element counter 
        integer(ip)             :: n            ! loop index for bodies
        integer(ip)             :: i            ! loop index for body pos points
        integer(ip)             :: j            ! loop index for body pos points
        integer(ip)             :: k            ! loop index for body nbr points
        integer(ip)             :: ix           ! x coord. mesh index of nbr pt
        integer(ip)             :: jy           ! y coord. mesh index of nbr pt

        ! square of mesh spacing
        ds2 = ds**2

        ! Create A matrix for finding velocity corrections, Ax=b
        this % a % nnz = 0_ip
        cnt = 0_ip

        do n = 1, this % num_body()
            body => this % body(n)
            do i = 1, body % num_pos()
                do j = 1, body % num_pos()
                    aij = 0.0_wp
                    do k = 1, body % nbrs(i) % num 
                        kerni = kernel(body % nbrs(i) % pos(k), body % pos(i), ds)
                        kernj = kernel(body % nbrs(i) % pos(k), body % pos(j), ds)
                        aij = aij + kerni*kernj*ds2
                    end do
                    if (aij > 0.0_wp) then
                        cnt = cnt + 1_ip
                        this % a % ix(cnt) = i
                        this % a % jy(cnt) = j
                        this % a % val(cnt) = aij
                    end if
                end do
            end do
        end do
        this % a % nnz = cnt

        ! Create b vector for finding velocity corrections
        cnt = 0_ip
        do n = 1, this % num_body()
            body => this % body(n)
            do i = 1, body % num_pos()
                cnt = cnt + 1
                this % bx(cnt) = body % vel(i) % x
                this % by(cnt) = body % vel(i) % y
                do k = 1, body % nbrs(i) % num
                    kerni = kernel(body % nbrs(i) % pos(k), body % pos(i), ds)
                    this % bx(cnt) = this % bx(cnt) - kerni * ds2*body % nbrs(i) % u(k) % x 
                    this % by(cnt) = this % by(cnt) - kerni * ds2*body % nbrs(i) % u(k) % y 
                end do
            end do
        end do

        ! Configure minres linear solver
        minres_ez % checka = .true.
        minres_ez % precon = .true.
        !call minres_ez % print()

        ! Solve linear system A * vx = bx
        call minres_ez % solve(   & 
            this % a % ix,        & 
            this % a % jy,        &
            this % a % val,       &
            this % bx,            &
            this % vx,            &
            minres_info,          &
            this % a % nnz        &
            )
        !call minres_info % print()
        if (minres_info % istop >= 2) then
            print *, 'minres istop > 2 for dux velocity correction'
            stop 
        end if

        ! Solve linear system A * vy = by
        call minres_ez % solve(   & 
            this % a % ix,        & 
            this % a % jy,        &
            this % a % val,       &
            this % by,            &
            this % vy,            &
            minres_info,          &
            this % a % nnz        &
            )
        !call minres_info % print()
        if (minres_info % istop >= 2) then
            print *, 'minres istop > 2 for  duy velocity correction'
            stop 
        end if

        ! Apply velocity corrections to mesh
        cnt = 0_ip
        do n = 1, this % num_body()
            body => this % body(n)
            do i = 1, body % num_pos()
                cnt = cnt + 1
                do k = 1, body % nbrs(i) % num
                    kerni = kernel(body % nbrs(i) % pos(k), body % pos(i), ds)
                    ix = body % nbrs(i) % ix(k)
                    jy = body % nbrs(i) % jy(k)
                    u(ix, jy) % x = u(ix,jy) % x + (kerni * this % vx(cnt))
                    u(ix, jy) % y = u(ix,jy) % y + (kerni * this % vy(cnt))
                end do
            end do
        end do

    end subroutine ibsol_corrector


    function ibsol_num_body(this) result(num)
        class(ibsol_t), intent(in) :: this
        integer(ip)                :: num
        if (.not. allocated(this % body)) then
            num = 0_ip
        else
            num = size(this % body)
        end if
    end function ibsol_num_body


    subroutine bounding_ind(body, mesh, ds, imin, imax, jmin, jmax)
        type(body_t), intent(in)  :: body
        type(mesh_t), intent(in)  :: mesh
        real(wp),     intent(in)  :: ds
        integer(ip),  intent(out) :: imin
        integer(ip),  intent(out) :: imax
        integer(ip),  intent(out) :: jmin
        integer(ip),  intent(out) :: jmax
        real(wp)                  :: xmin
        real(wp)                  :: xmax
        real(wp)                  :: ymin
        real(wp)                  :: ymax
        integer(ip)               :: num_x
        integer(ip)               :: num_y
        num_x = mesh % num_x()
        num_y = mesh % num_y()
        call body % bounding_box(xmin, xmax, ymin, ymax)
        imin = floor(xmin/ds) + 1
        jmin = floor(ymin/ds) + 1
        imax = ceiling(xmax/ds) + 1
        jmax = ceiling(ymax/ds) + 1
        imin = constrain(imin, 1, num_x)
        imax = constrain(imax, 1, num_x)
        jmin = constrain(jmin, 1, num_y)
        jmax = constrain(jmax, 1, num_y)
    end subroutine bounding_ind


end module slbm_2d_ibsol
