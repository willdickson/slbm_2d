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

    function cond_id_from_name(cond_name) result(id)
        character(*), intent(in) :: cond_name
        integer(ip)              :: id
        select case(cond_name)
        case ( 'inflow' )
            id = BNDRY_COND_INFLOW
        case ( 'moving' )
            id = BNDRY_COND_MOVING
        case ( 'outflow' )
            id = BNDRY_COND_OUTFLOW
        case ( 'noslip' )
            id = BNDRY_COND_NOSLIP
        case ( 'slip' )
            id = BNDRY_COND_SLIP
        case default
            id = BNDRY_COND_UNKNOWN
        end select
    end function cond_id_from_name


    function cond_id_to_name(type_id) result(cond_name)
        integer(ip), intent(in)   :: type_id
        character(:), allocatable :: cond_name
        select case(type_id)
        case ( BNDRY_COND_INFLOW )
            cond_name = 'inflow'
        case ( BNDRY_COND_MOVING )
            cond_name = 'moving'
        case ( BNDRY_COND_OUTFLOW )
            cond_name = 'outflow'
        case ( BNDRY_COND_NOSLIP )
            cond_name = 'noslip'
        case ( BNDRY_COND_SLIP )
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
        case ( 'left' )
            side_id = BNDRY_SIDE_LEFT
        case ( 'right' )
            side_id = BNDRY_SIDE_RIGHT
        case ( 'top' )
            side_id = BNDRY_SIDE_TOP
        case ( 'bottom' )
            side_id = BNDRY_SIDE_BOTTOM
        case default
            side_id = BNDRY_SIDE_UNKNOWN
        end select
    end function side_id_from_name


end module slbm_2d_bndry
