module slbm_2d_velocity

    use slbm_2d_kinds,  only : wp, ip

    use slbm_2d_const,  only : NUM_BNDRY
    use slbm_2d_const,  only : BNDRY_SIDE_LEFT
    use slbm_2d_const,  only : BNDRY_SIDE_RIGHT
    use slbm_2d_const,  only : BNDRY_SIDE_TOP
    use slbm_2d_const,  only : BNDRY_SIDE_BOTTOM

    use slbm_2d_const,  only : BNDRY_COND_INFLOW
    use slbm_2d_const,  only : BNDRY_COND_MOVING
    use slbm_2d_const,  only : BNDRY_COND_OUTFLOW
    use slbm_2d_const,  only : BNDRY_COND_NOSLIP
    use slbm_2d_const,  only : BNDRY_COND_SLIP

    use slbm_2d_init,   only : init_file_t
    use slbm_2d_init,   only : init_const_t
    use slbm_2d_vector, only : vector_t
    use slbm_2d_config, only : config_t
    use slbm_2d_bndry,  only : bndry_t

    implicit none
    private

    type, public :: velocity_t
        type(vector_t), allocatable :: curr(:,:)
        type(vector_t), allocatable :: pred(:,:)
    contains
        private
        procedure, public  :: get_bndry_ind
        procedure, public  :: set_initial_cond
        procedure, public  :: deallocate => velocity_deallocate
        final              :: velocity_destructor
    end type velocity_t

    interface velocity_t
        procedure velocity_constructor
    end interface velocity_t

contains

    function velocity_constructor(config) result(velocity)
        type(config_t), intent(in) :: config
        type(velocity_t)        :: velocity
        allocate(velocity % curr(config % num_x, config % num_y))
        allocate(velocity % pred(config % num_x, config % num_y))
        call velocity % set_initial_cond(config)
    end function velocity_constructor


    subroutine set_initial_cond(this, config)
        class(velocity_t), intent(inout) :: this
        type(config_t), intent(in)       :: config
        integer(ip)                      :: i

        ! Initialize velocity field
        select type( init => config % init )
        class is ( init_const_t )
            this % curr = init % velocity
            this % pred = vector_t(0.0_wp, 0.0_wp)
        class is ( init_file_t )
            print *, "setting initial cond. from file net implemented yet"
            stop
        class default
            print *, "error unknown method for settign initial cond."
            stop
        end select

        ! Ensure boundry values conform with boundary cond.
        block
            type(bndry_t), pointer :: bndry
            type(vector_t)         :: vel
            integer(ip)            :: ix1, ixn  
            integer(ip)            :: iy1, iyn 

            ! Set Inflow and slip boundaries first
            do i=1,NUM_BNDRY
                bndry => config % bndry_cond(i) % ptr
                call this % get_bndry_ind(bndry % side_id, ix1, ixn, iy1, iyn)
                select case(bndry % cond_id)
                case ( BNDRY_COND_INFLOW )
                    this % curr(ix1:ixn,iy1:iyn) = bndry % velocity
                case ( BNDRY_COND_SLIP )
                    select case(bndry % side_id)
                    case (BNDRY_SIDE_LEFT, BNDRY_SIDE_RIGHT )
                        this % curr(ix1:ixn, iy1:iyn) % x = 0.0_wp
                    case (BNDRY_SIDE_TOP, BNDRY_SIDE_BOTTOM )
                        this % curr(ix1:ixn, iy1:iyn) % y = 0.0_wp
                    end select
                end select
            end do

            ! Set noslip and moving boundaries second 
            do i=1,NUM_BNDRY
                bndry => config % bndry_cond(i) % ptr
                call this % get_bndry_ind(bndry % side_id, ix1, ixn, iy1, iyn)
                select case(bndry % cond_id)
                case ( BNDRY_COND_NOSLIP, BNDRY_COND_MOVING )
                    this % curr(ix1:ixn,iy1:iyn) = bndry % velocity
                end select

            end do
        end block

    end subroutine set_initial_cond


    subroutine get_bndry_ind(this, loc_id, ix1, ixn, iy1, iym)
        class(velocity_t), intent(in) :: this
        integer(ip), intent(in)       :: loc_id 
        integer(ip), intent(out)      :: ix1
        integer(ip), intent(out)      :: ixn
        integer(ip), intent(out)      :: iy1 
        integer(ip), intent(out)      :: iym 
        integer(ip)                   :: num_x
        integer(ip)                   :: num_y
        num_x = size(this % curr, 1)
        num_y = size(this % curr, 2)
        select case(loc_id)
        case ( BNDRY_SIDE_LEFT )
            ix1 = 1_ip
            ixn = 1_ip
            iy1 = 1_ip
            iym = num_y
        case ( BNDRY_SIDE_RIGHT )
            ix1 = num_x
            ixn = num_x 
            iy1 = 1_ip
            iym = num_y
        case ( BNDRY_SIDE_TOP )
            ix1 = 1_ip
            ixn = num_x
            iy1 = 1_ip
            iym = 1_ip 
        case ( BNDRY_SIDE_BOTTOM )
            ix1 = 1_ip
            ixn = num_x
            iy1 = num_y
            iym = num_y
        case default
            print *, 'get_bndry_indices, unknown boundary loc_id ', loc_id
            stop
        end select
    end subroutine get_bndry_ind


    subroutine velocity_deallocate(this)
        class(velocity_t), intent(inout) :: this
        if ( allocated(this % curr) ) then
            deallocate(this % curr)
        end if
        if ( allocated(this % pred) ) then
            deallocate(this % pred)
        end if
    end subroutine velocity_deallocate


    subroutine velocity_destructor(this)
        type(velocity_t), intent(inout) :: this
        call this % deallocate()
    end subroutine velocity_destructor


end module slbm_2d_velocity
