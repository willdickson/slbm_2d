module slbm_2d_init

    use slbm_2d_kinds,  only : wp, ip

    use slbm_2d_const,  only : INIT_COND_CONST
    use slbm_2d_const,  only : INIT_COND_FILE
    use slbm_2d_const,  only : INIT_COND_UNKNOWN

    use slbm_2d_vector, only : vector_t

    implicit none
    private

    type, public :: init_t
        character(:), allocatable :: name 
        integer(ip)               :: id
    end type init_t

    type, public, extends(init_t) :: init_const_t
        type(vector_t)  :: velocity
    end type init_const_t

    type, public, extends(init_t) :: init_file_t
        character(:), allocatable :: filename
    end type init_file_t

    public :: init_id_from_string
    public :: init_id_to_string

contains

    function init_id_from_string(init_string) result(init_id)
        character(*), intent(in) :: init_string
        integer(ip)              :: init_id
        select case (init_string)
        case ( 'constant' )
            init_id = INIT_COND_CONST
        case ( 'file' )
            init_id = INIT_COND_FILE
        case default
            init_id = INIT_COND_UNKNOWN
        end select
    end function init_id_from_string


    function init_id_to_string(init_id) result(init_string)
        integer(ip), intent(in)   :: init_id
        character(:), allocatable :: init_string
        select case (init_id)
        case ( INIT_COND_CONST )
            init_string = 'constant'
        case ( INIT_COND_FILE )
            init_string = 'file'
        case default
            init_string = 'unknown'
        end select
    end function init_id_to_string

end module slbm_2d_init
