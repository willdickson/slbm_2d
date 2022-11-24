module slbm_2d_distrib
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_const,  only : LATTICE_Q
    use slbm_2d_config, only : config_t
    implicit none
    private

    type, public :: distrib_t
        real(wp), allocatable :: val(:,:,:)
    contains
        private
        procedure, public :: deallocate => distrib_deallocate
        final             :: distrib_destructor
    end type distrib_t

    interface distrib_t
        procedure distrib_constructor
    end interface distrib_t

contains

    function distrib_constructor(config) result(distrib)
        type(config_t), intent(in) :: config
        type(distrib_t)            :: distrib
        allocate( distrib % val(LATTICE_Q, config % num_x, config % num_y) )
        distrib % val = 0.0_wp
    end function distrib_constructor

    subroutine distrib_deallocate(this)
        class(distrib_t), intent(inout) :: this
        if ( allocated(this % val) ) then
            deallocate(this % val)
        end if
    end subroutine distrib_deallocate

    subroutine distrib_destructor(this)
        type(distrib_t), intent(inout) :: this
        call this % deallocate()
    end subroutine distrib_destructor

end module slbm_2d_distrib
