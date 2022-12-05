module slbm_2d

    use slbm_2d_kinds,      only : wp, ip
    use slbm_2d_config,     only : config_t
    use slbm_2d_simulation, only : simulation_t

    use slbm_2d_body,       only : body_t
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
        integer(ip), parameter      :: num_pts = 100

        character(:), allocatable   :: filename
        type(config_t)              :: config
        type(simulation_t)          :: sim
        type(body_t)                :: body
        type(vector_t), allocatable :: pts(:)
        real(wp), allocatable       :: s(:)
        integer(ip)                 :: i,j,k

        filename = './example/config.toml'
        config = config_t(filename)
        call config % pprint()
        sim = simulation_t(config)

        s = [((s0 + (s1 - s0)*k/(num_pts - 1.0_wp)), k = 0, num_pts-1)]
        allocate(pts(size(s)))
        pts % x = x0 + amp_x*cos(2.0_wp*PI*s)
        pts % y = y0 + amp_y*sin(2.0_wp*PI*s)
        body = body_t(BODY_TYPE_OPEN, pts)

        call body % update(sim % curr, sim % mesh, config % ds, 0.0_wp)
        do k = 1, body % num_pos()
            print *, k, body % nbrs(k) % num
            print *, body % pos(k) % x / config % ds, body % pos(k) % y / config % ds
            do i = 1, body % nbrs(k) % num
                print *, '      ', body % nbrs(k) % ix(i), body % nbrs(k) % iy(i)
            end do
        end do 

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
