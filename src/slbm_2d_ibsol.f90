module slbm_2d_ibsol

    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_vector, only : vector_t
    use slbm_2d_body,   only : body_t
    use slbm_2d_mesh,   only : mesh_t
    use slbm_2d_spmat,  only : spmat_t
    use slbm_2d_state,  only : state_t
    use minres,         only : minres_ez_t
    use minres,         only : minres_info_t

    implicit none
    private

    type, public :: ibsol_t
        type(body_t), allocatable :: body(:)   ! Array of immersed bodies
        type(spmat_t)             :: a         ! A matrix for velocity correction 
        real(wp), allocatable     :: bx(:)     ! b vector for velocity corr. x comp.
        real(wp), allocatable     :: by(:)     ! b vector for velocity corr. y comp. 
        real(wp), allocatable     :: vx(:)     ! velocity correction x component 
        real(wp), allocatable     :: vy(:)     ! velocity correction y component
        integer(ip)               :: num_pos   ! total number of pos pts (all bodies)
    contains
        private
        procedure, public :: update    => ibsol_update
        procedure, public :: corrector => ibsol_corrector
        procedure, public :: num_body  => ibsol_num_body
    end type ibsol_t

    interface ibsol_t
        procedure :: ibsol_constructor
    end interface ibsol_t

contains

    function ibsol_constructor(body) result(ibsol)
        type(body_t), intent(in) :: body(:)
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
    end function ibsol_constructor


    subroutine ibsol_update(this, state, mesh, ds, time)
        class(ibsol_t), intent(inout) :: this    ! immersed boundry solver
        type(state_t), intent(in)     :: state   ! fluid state: ux, uy, rho
        type(mesh_t), intent(in)      :: mesh    ! x and y meshes
        real(wp), intent(in)          :: ds      ! mesh spacing
        real(wp), intent(in)          :: time    ! simulation time 
        integer(ip)                   :: i
        do i = 1, this % num_body()
            call this % body(i) % update(state, mesh, ds, time)
        end do
    end subroutine ibsol_update


    subroutine ibsol_corrector(this, mesh, ds, du)
        class(ibsol_t), intent(inout) :: this    ! immersed boundry solver
        type(mesh_t), intent(in)      :: mesh    ! x and  y meshes
        real(wp), intent(in)          :: ds      ! mesh spacing
        type(vector_t), intent(out)   :: du(:)   ! correction velocities

        type(minres_ez_t)             :: minres_ez
        type(minres_info_t)           :: minres_info
        real(wp)                      :: kval
        real(wp)                      :: kval1
        real(wp)                      :: kval2
        real(wp)                      :: aij
        integer(ip)                   :: cnt
        integer(ip)                   :: i,j,k

        this % a % nnz = 0_ip
        cnt = 0_ip



    end subroutine ibsol_corrector


    function ibsol_num_body(this) result(num)
        class(ibsol_t), intent(in) :: this
        integer(ip)                :: num
        num = size(this % body)
    end function ibsol_num_body

end module slbm_2d_ibsol
