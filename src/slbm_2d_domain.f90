module slbm_2d_domain

    use slbm_2d_kinds,     only : wp, ip
    use slbm_2d_config,    only : config_t
    use slbm_2d_vector,    only : vector_t
    use slbm_2d_density,   only : density_t
    use slbm_2d_distrib,   only : distrib_t
    use slbm_2d_velocity,  only : velocity_t
    implicit none
    private

    type, public :: domain_t
        type(distrib_t)     :: equilib 
        type(density_t)     :: density 
        type(velocity_t)    :: velocity 
    contains
        private
    end type domain_t

    interface domain_t
        procedure domain_constructor
    end interface domain_t
        
contains

    function domain_constructor(config) result(domain)
        type(config_t), intent(in) :: config
        type(domain_t)             :: domain
        domain % equilib  = distrib_t(config)
        domain % density  = density_t(config) 
        domain % velocity = velocity_t(config) 
    end function domain_constructor

end module slbm_2d_domain
