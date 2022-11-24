module slbm_2d_bndry
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_const,  only : BNDRY_COND_UNKNOWN
    use slbm_2d_vector, only : vector_t
    implicit none
    private

    type, public :: bndry_t
        integer(ip)    :: id  = BNDRY_COND_UNKNOWN
        type(vector_t) :: velocity
    end type bndry_t

    type, public :: bndry_ptr_t
        type(bndry_t), pointer :: ptr 
    end type bndry_ptr_t

    interface bndry_t
        procedure bndry_constructor
    end interface bndry_t

contains

    function bndry_constructor(id, velocity) result(bndry)
        integer(ip), intent(in)    :: id
        type(vector_t), intent(in) :: velocity
        type(bndry_t)              :: bndry
        bndry % id = id
        bndry % velocity = velocity
    end function bndry_constructor

end module slbm_2d_bndry
