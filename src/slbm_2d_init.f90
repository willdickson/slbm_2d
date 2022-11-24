module slbm_2d_init
    use slbm_2d_kinds,  only : wp, ip
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


end module slbm_2d_init
