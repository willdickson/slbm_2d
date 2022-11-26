module slbm_2d_bndry

    use slbm_2d_kinds,  only : wp, ip

    use slbm_2d_const,  only : BNDRY_SIDE_LEFT
    use slbm_2d_const,  only : BNDRY_SIDE_RIGHT
    use slbm_2d_const,  only : BNDRY_SIDE_TOP
    use slbm_2d_const,  only : BNDRY_SIDE_BOTTOM
    use slbm_2d_const,  only : BNDRY_SIDE_UNKNOWN

    use slbm_2d_const,  only : BNDRY_COND_INFLOW
    use slbm_2d_const,  only : BNDRY_COND_MOVING
    use slbm_2d_const,  only : BNDRY_COND_OUTFLOW
    use slbm_2d_const,  only : BNDRY_COND_NOSLIP
    use slbm_2d_const,  only : BNDRY_COND_SLIP
    use slbm_2d_const,  only : BNDRY_COND_UNKNOWN

    use slbm_2d_vector, only : vector_t

    implicit none
    private

    type, public :: bndry_t
        integer(ip)    :: side_id = BNDRY_SIDE_UNKNOWN
        integer(ip)    :: cond_id = BNDRY_COND_UNKNOWN
        type(vector_t) :: velocity
    contains
        private
        procedure, public :: get_indices
        procedure, public :: get_offsets
    end type bndry_t

    interface bndry_t
        procedure bndry_constructor
    end interface bndry_t

    type, public :: bndry_ptr_t
        type(bndry_t), pointer :: ptr 
    end type bndry_ptr_t

    public :: cond_id_from_name
    public :: side_id_from_name
    public :: cond_id_to_name
    public :: side_id_to_name

contains

    function bndry_constructor(side_id, cond_id, velocity) result(bndry)
        integer(ip), intent(in)    :: side_id
        integer(ip), intent(in)    :: cond_id
        type(vector_t), intent(in) :: velocity
        type(bndry_t)              :: bndry
        bndry % side_id = side_id
        bndry % cond_id = cond_id
        bndry % velocity = velocity
    end function bndry_constructor


    subroutine get_indices(this, num_x, num_y, i1, i2, j1, j2)
        class(bndry_t), intent(in) :: this   ! the bndry object
        integer(ip), intent(in)    :: num_x  ! mesh size in x dim
        integer(ip), intent(in)    :: num_y  ! mesh size in y dim
        integer(ip), intent(out)   :: i1     ! 1st index of x range
        integer(ip), intent(out)   :: i2     ! 2nd index of x range
        integer(ip), intent(out)   :: j1     ! 1st index of y range
        integer(ip), intent(out)   :: j2     ! 2nd index of y range

        select case (this % side_id)
        case (BNDRY_SIDE_LEFT)
            i1 = 1_ip
            i2 = 1_ip
            j1 = 1_ip
            j2 = num_y
        case (BNDRY_SIDE_RIGHT)
            i1 = num_x
            i2 = num_x 
            j1 = 1_ip
            j2 = num_y
        case (BNDRY_SIDE_TOP)
            i1 = 1_ip
            i2 = num_x
            j1 = num_y
            j2 = num_y
        case (BNDRY_SIDE_BOTTOM)
            i1 = 1_ip
            i2 = num_x
            j1 = 1_ip
            j2 = 1_ip 
        case default
            print *, 'get_indices, unknown boundary side_id = ', &
                this % side_id
            stop
        end select

        select case (this % cond_id)
        case (BNDRY_COND_INFLOW, BNDRY_COND_OUTFLOW)
            select case (this % side_id)
            case (BNDRY_SIDE_LEFT, BNDRY_SIDE_RIGHT)
                j1 = j1 + 1
                j2 = j2 - 1
            case (BNDRY_SIDE_TOP, BNDRY_SIDE_BOTTOM)
                i1 = i1 + 1
                i2 = i2 - 1
            end select
        end select

    end subroutine get_indices


    subroutine get_offsets(this, ik, jk)
        class(bndry_t), intent(in) :: this  ! the bndry object
        integer(ip), intent(out)   :: ik    ! offset in x direction
        integer(ip), intent(out)   :: jk    ! offset in y direction
        select case (this % side_id)
        case (BNDRY_SIDE_LEFT)
            ik =  1_ip
            jk =  0_ip
        case (BNDRY_SIDE_RIGHT)
            ik = -1_ip
            jk =  0_ip
        case (BNDRY_SIDE_TOP)
            ik =  0_ip
            jk = -1_ip
        case (BNDRY_SIDE_BOTTOM)
            ik =  0_ip
            jk =  1_ip
        case default
            print *, 'get_offsets, unknown boundary side_id = ', & 
                this % side_id
            stop
        end select
    end subroutine get_offsets


    function cond_id_from_name(cond_name) result(id)
        character(*), intent(in) :: cond_name
        integer(ip)              :: id
        select case(cond_name)
        case ('inflow')
            id = BNDRY_COND_INFLOW
        case ('moving')
            id = BNDRY_COND_MOVING
        case ('outflow')
            id = BNDRY_COND_OUTFLOW
        case ('noslip')
            id = BNDRY_COND_NOSLIP
        case ('slip')
            id = BNDRY_COND_SLIP
        case default
            id = BNDRY_COND_UNKNOWN
        end select
    end function cond_id_from_name


    function cond_id_to_name(type_id) result(cond_name)
        integer(ip), intent(in)   :: type_id
        character(:), allocatable :: cond_name
        select case(type_id)
        case (BNDRY_COND_INFLOW)
            cond_name = 'inflow'
        case (BNDRY_COND_MOVING)
            cond_name = 'moving'
        case (BNDRY_COND_OUTFLOW)
            cond_name = 'outflow'
        case (BNDRY_COND_NOSLIP)
            cond_name = 'noslip'
        case (BNDRY_COND_SLIP)
            cond_name = 'slip'
        case default
            cond_name = 'unknown'
        end select
    end function cond_id_to_name


    function side_id_to_name(side_id) result(side_name)
        integer(ip), intent(in)   :: side_id
        character(:), allocatable :: side_name
        select case(side_id)
        case ( BNDRY_SIDE_LEFT )
            side_name = 'left'
        case ( BNDRY_SIDE_RIGHT )
            side_name = 'right'
        case ( BNDRY_SIDE_TOP )
            side_name = 'top'
        case ( BNDRY_SIDE_BOTTOM )
            side_name = 'bottom'
        case default
            side_name = 'unknown'
        end select
    end function side_id_to_name


    function side_id_from_name(side_name) result(side_id)
        character(*), intent(in) :: side_name
        integer(ip)              :: side_id
        select case(side_name)
        case ('left')
            side_id = BNDRY_SIDE_LEFT
        case ('right')
            side_id = BNDRY_SIDE_RIGHT
        case ('top')
            side_id = BNDRY_SIDE_TOP
        case ('bottom')
            side_id = BNDRY_SIDE_BOTTOM
        case default
            side_id = BNDRY_SIDE_UNKNOWN
        end select
    end function side_id_from_name


end module slbm_2d_bndry
