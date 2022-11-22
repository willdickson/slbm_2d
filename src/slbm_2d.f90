module slbm_2d
    use slbm_2d_kinds,      only : wp, ip
    use slbm_2d_config,     only : config_t
    use slbm_2d_simulation, only : simulation_t

    implicit none
    private
    public :: testing

contains

    subroutine testing
        !integer(ip)        :: num_x = 101_ip
        !integer(ip)        :: num_y = 101_ip
        !real(wp)           :: ds    = 0.01_wp/5.0_wp
        !real(wp)           :: kvisc = 1.0e-4_wp
        !real(wp)           :: density = 1.0_wp
        !real(wp)           :: tstop = 25.0_wp

        character(:), allocatable :: filename
        type(config_t)            :: config
        type(simulation_t)        :: sim

        filename = './example/config.toml'
        config = config_t(filename)
        call config % pprint()

        !config  = config_t(num_x, num_y, ds, kvisc, density, tstop)
        !call config % pprint()

        !sim = simulation_t(config)
        !print *, sim % domain % velocity % num_x, sim % domain % velocity % num_y

    end subroutine testing

end module slbm_2d
