module slbm_2d

    use slbm_2d_config,     only : config_t
    use slbm_2d_simulation, only : simulation_t

    use slbm_2d_line_seg,   only : line_seg_t
    use slbm_2d_line_seg,   only : intersect 
    use slbm_2d_line_seg,   only : ccw 
    use slbm_2d_vector,     only : vector_t

    implicit none
    private

    public :: sim_test
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


    subroutine line_seg_test
        type(vector_t)   :: a = vector_t(0.0, 0.0)
        type(vector_t)   :: b = vector_t(1.0, 1.0)
        type(vector_t)   :: c = vector_t(1.0, 0.0)
        type(vector_t)   :: d = vector_t(0.0, 1.0)
        type(line_seg_t) :: s1
        type(line_seg_t) :: s2

        s1 = line_seg_t(a,b)
        s2 = line_seg_t(c,d)

        print *, intersect(s1,s2)




    end subroutine line_seg_test

end module slbm_2d
