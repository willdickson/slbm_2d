module slbm_2d
    use slbm_2d_kinds,      only : wp, ip
    use slbm_2d_config,     only : config_t
    use slbm_2d_vector,     only : vector_t
    use slbm_2d_simulation, only : simulation_t

    implicit none
    private
    public :: testing

contains

    subroutine testing

        character(:), allocatable  :: filename
        type(config_t)             :: config
        type(simulation_t)         :: sim


        filename = './example/config.toml'
        config = config_t(filename)
        call config % pprint()

        sim = simulation_t(config)
        call sim % run()


        !print *, sim % config % init % id
        !print *, sim % config % bndry_left % id
        !print *, sim % config % bndry_array(1) % ptr % id


    end subroutine testing

end module slbm_2d
