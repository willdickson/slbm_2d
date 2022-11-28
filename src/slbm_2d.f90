module slbm_2d
    use slbm_2d_config,     only : config_t
    use slbm_2d_simulation, only : simulation_t
    use m_system,           only : system_dir

    implicit none
    private
    public :: testing

contains

    subroutine testing

        integer                    :: i
        character(:), allocatable  :: filelist(:)

        filelist = system_dir('src', pattern='*.f90')
        print *, size(filelist)
        do i = 1, size(filelist)
            print *, filelist(i)
        end do

        !character(:), allocatable  :: filename
        !type(config_t)             :: config
        !type(simulation_t)         :: sim

        !filename = './example/config.toml'
        !config = config_t(filename)
        !call config % pprint()

        !sim = simulation_t(config)
        !call sim % run()

    end subroutine testing

end module slbm_2d
