module slbm_2d_vector
    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    type, public :: vector_t
        real(wp) :: x = 0.0_wp
        real(wp) :: y = 0.0_wp
    contains
        private
        procedure, pass(this) :: add
        generic, public       :: operator(+) => add
    end type

contains

    elemental function add(this,v) result(w)
        class(vector_t), intent(in) :: this
        type(vector_t),  intent(in) :: v 
        type(vector_t)              :: w
        w % x = this % x + v % x
        w % y = this % y + v % y
    end function add

end module slbm_2d_vector
