module slbm_2d

    use slbm_2d_kinds,      only : wp, ip
    use slbm_2d_config,     only : config_t
    use slbm_2d_simulation, only : simulation_t

    use slbm_2d_ibsol,      only : ibsol_t
    use slbm_2d_body,       only : body_t
    use slbm_2d_const,      only : VECTOR_ZERO
    use slbm_2d_const,      only : BODY_TYPE_OPEN
    use slbm_2d_const,      only : PI
    use slbm_2d_line_seg,   only : line_seg_t
    use slbm_2d_line_seg,   only : intersect 
    use slbm_2d_line_seg,   only : orientation 
    use slbm_2d_line_seg,   only : is_chain 
    use slbm_2d_vector,     only : vector_t

    implicit none
    private

    public :: sim_test
    public :: body_test
    public :: line_seg_test

contains

    subroutine sim_test
        character(:), allocatable  :: filename
        type(config_t)             :: config
        type(simulation_t)         :: sim

        filename = './example/config.toml'
        config = config_t(filename)
        call config % pprint()
        sim = simulation_t(config)
        call sim % run()
    end subroutine sim_test


    subroutine body_test
        real(wp), parameter         :: s0 = 0.0_wp
        real(wp), parameter         :: s1 = 0.5_wp
        real(wp), parameter         :: x0 = 0.5_wp
        real(wp), parameter         :: y0 = 0.5_wp
        real(wp), parameter         :: amp_x = 0.2_wp
        real(wp), parameter         :: amp_y = 0.2_wp
        integer(ip), parameter      :: num_pts = 50_ip 
        integer(ip), parameter      :: num_body = 1_ip

        character(:), allocatable   :: filename
        type(config_t)              :: config
        type(simulation_t)          :: sim
        type(body_t), allocatable   :: body(:)
        type(vector_t), allocatable :: pts(:)
        type(vector_t), allocatable :: du(:)
        real(wp), allocatable       :: s(:)
        real(wp), allocatable       :: a(:,:)
        real(wp)                    :: aij
        integer(ip)                 :: ix, jy
        integer(ip)                 :: i,j,k

        filename = './example/config.toml'
        config = config_t(filename)
        call config % pprint()
        sim = simulation_t(config)

        s = [((s0 + (s1 - s0)*k/(num_pts - 1.0_wp)), k = 0, num_pts-1)]
        allocate(pts(size(s)))
        pts % x = x0 + amp_x*cos(2.0_wp*PI*s)
        pts % y = y0 + amp_y*sin(2.0_wp*PI*s)
        allocate(body(num_body))
        body(1) = body_t(BODY_TYPE_OPEN, pts)
        sim % ibsol  = ibsol_t(body)

        allocate(du(size(pts)), source=VECTOR_ZERO) 

        call sim % ibsol % update(sim % curr, sim % mesh, config % ds, 0.0_wp)
        call sim % ibsol % corrector(config % ds, sim % curr % u)

        !print *, 'nnz = ', sim % ibsol % a % nnz
        !print *, 'density = ', sim % ibsol % a % density()

        !do k = 1, sim % ibsol % a % nnz
        !    ix  = sim % ibsol % a % ix(k)
        !    jy  = sim % ibsol % a % jy(k)
        !    aij = sim % ibsol % a % val(k)
        !    print *, ix, jy, aij 
        !end do


        !a = sim % ibsol % a % as_dense_mat()

        !do i = 1, size(a,1) 
        !    do j = 1, size(a,2)
        !        print *, i, j, a(i,j), a(i,j) - a(j,i)
        !    end do
        !    print *, ''
        !end do


    end subroutine body_test


    subroutine line_seg_test
        type(vector_t)   :: a = vector_t(0.0, 0.0)
        type(vector_t)   :: b = vector_t(1.0, 0.0)
        type(vector_t)   :: c = vector_t(0.0, 0.0)
        type(vector_t)   :: d = vector_t(1.0, 1.0)
        type(line_seg_t) :: s1
        type(line_seg_t) :: s2

        s1 = line_seg_t(a,b)
        s2 = line_seg_t(c,d)
        print *, intersect(s1, s2), is_chain(s1, s2)
    end subroutine line_seg_test


end module slbm_2d
