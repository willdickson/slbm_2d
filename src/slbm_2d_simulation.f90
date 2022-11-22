module slbm_2d_simulation
    use slbm_2d_kinds,   only : wp, ip
    use slbm_2d_config,  only : config_t
    use slbm_2d_domain,  only : domain_t
    implicit none
    private

    type, public :: simulation_t
        type(config_t)  :: config  ! configuration data
        type(domain_t)  :: domain  ! simulation domain
    contains
        private
        !procedure, public  :: run
        !procedure, public :: predictor_update
        !procedure, public :: corrector_update
    end type simulation_t

    interface simulation_t
        procedure simulation_constructor
    end interface simulation_t

contains

    function simulation_constructor(config) result(simulation)
        type(config_t), intent(in) :: config
        type(simulation_t)         :: simulation
        simulation % config = config
        simulation % domain = domain_t(config)
    end function simulation_constructor
    
end module slbm_2d_simulation
