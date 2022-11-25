module slbm_2d_density
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_config, only : config_t
    use slbm_2d_vector, only : vector_t
    implicit none
    private

    type, public :: density_t
        real(wp), allocatable :: curr(:,:)
        real(wp), allocatable :: pred(:,:)

    contains
        private
        procedure         :: set_initial_cond
        procedure, public :: deallocate => density_deallocate
        final    :: density_destructor 
    end type density_t

    interface density_t
        procedure density_constructor
    end interface density_t


contains

    function density_constructor(config) result(density)
        type(config_t), intent(in) :: config
        type(density_t)            :: density
        allocate(density % curr(config % num_x, config % num_y))
        allocate(density % pred(config % num_x, config % num_y))
        call density % set_initial_cond(config)
    end function density_constructor


    subroutine set_initial_cond(this, config)
        class(density_t), intent(inout) :: this
        type(config_t), intent(in)      :: config
        this % curr = config % density
        this % pred = 0.0_wp
    end subroutine set_initial_cond


    subroutine density_deallocate(this)
        class(density_t), intent(inout) :: this
        if ( allocated(this % curr) ) then
            deallocate(this % curr)
        end if
        if ( allocated(this % pred) ) then
            deallocate(this % pred)
        end if
    end subroutine density_deallocate


    subroutine density_destructor(this)
        type(density_t), intent(inout) :: this
        call this % deallocate()
    end subroutine density_destructor


end module slbm_2d_density
