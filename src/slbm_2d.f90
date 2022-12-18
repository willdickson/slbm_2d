module slbm_2d

    use slbm_2d_kinds,      only : wp, ip
    use slbm_2d_config,     only : config_t
    use slbm_2d_simulation, only : simulation_t

    use slbm_2d_ibsol,      only : ibsol_t
    use slbm_2d_body,       only : body_t
    use slbm_2d_const,      only : VECTOR_ZERO
    use slbm_2d_const,      only : BODY_TYPE_OPEN
    use slbm_2d_const,      only : BODY_TYPE_CLOSED
    use slbm_2d_const,      only : PI
    use slbm_2d_line_seg,   only : line_seg_t
    use slbm_2d_line_seg,   only : intersect 
    use slbm_2d_line_seg,   only : orientation 
    use slbm_2d_line_seg,   only : is_chain 
    use slbm_2d_vector,     only : vector_t
    use stdlib_io_npy,      only : save_npy

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
        real(wp), parameter         :: s0 =  0.0_wp
        real(wp), parameter         :: s1 =  1.0_wp
        real(wp), parameter         :: x0 =  2.0_wp
        real(wp), parameter         :: y0 =  1.0_wp
        real(wp), parameter         :: amp_x = 0.1_wp
        real(wp), parameter         :: amp_y = 0.1_wp
        integer(ip), parameter      :: num_pts = 50_ip 
        integer(ip), parameter      :: num_body = 1_ip

        character(:), allocatable   :: filename
        type(config_t)              :: config
        type(simulation_t)          :: sim
        type(body_t), allocatable   :: body(:)
        type(vector_t), allocatable :: pts(:)
        real(wp), allocatable       :: s(:)
        real(wp), allocatable       :: a(:,:)
        real(wp)                    :: aij
        integer(ip)                 :: ix, jy
        integer(ip)                 :: i,j,k

        filename = './example/config.toml'
        config = config_t(filename)
        call config % pprint()
        sim = simulation_t(config)

        s = [((s0 + (s1 - s0)*k/(num_pts - 1.0_wp)), k = 0, num_pts-2)]
        allocate(pts(size(s)))
        pts % x = x0 + amp_x*cos(2.0_wp*PI*s)
        pts % y = y0 + amp_y*sin(2.0_wp*PI*s)

        allocate(body(num_body))
        body(1) = body_t(BODY_TYPE_CLOSED, pts)
        !body(1) % vel = vector_t(0.0_wp, 0.0_wp)
        sim % ibsol  = ibsol_t(body, config % num_x, config % num_y)
        call sim % ibsol % update(sim % curr, sim % mesh, config % ds, 0.0_wp)

        !call sim % ibsol % corrector(config % ds, sim % curr % u)
        !call save_npy('id.npy', sim % ibsol % id)

        call sim % run()


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
