module slbm_2d_vector
    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    type, public :: vector_t
        real(wp) :: x = 0.0_wp
        real(wp) :: y = 0.0_wp
    contains
        private
        procedure, pass(this) :: add_vector
        procedure, pass(this) :: add_scalar
        procedure, pass(this) :: scalar_add
        procedure, pass(this) :: sub_vector
        procedure, pass(this) :: sub_scalar
        procedure, pass(this) :: scalar_sub
        procedure, pass(this) :: mul_vector
        procedure, pass(this) :: mul_scalar 
        procedure, pass(this) :: scalar_mul 
        procedure, pass(this) :: div_vector 
        procedure, pass(this) :: div_scalar 
        procedure, pass(this) :: assign_vector
        procedure, pass(this) :: assign_scalar
        
        generic, public ::   operator(+) => add_vector, &
                                            scalar_add, &
                                            add_scalar
        generic, public ::   operator(-) => sub_vector, &
                                            scalar_sub, &
                                            sub_scalar
        generic, public ::   operator(*) => mul_vector, &
                                            scalar_mul, &
                                            mul_scalar

        generic, public ::   operator(/) => div_vector, &
                                            div_scalar

        generic, public :: assignment(=) => assign_vector, &
                                            assign_scalar
    end type

    public :: dot
    public :: mag

contains

    elemental function add_vector(this,v) result(w)
        class(vector_t), intent(in) :: this
        type(vector_t),  intent(in) :: v 
        type(vector_t)              :: w
        w % x = this % x + v % x
        w % y = this % y + v % y
    end function add_vector


    elemental function scalar_add(a, this) result(w)
        real(wp), intent(in)        :: a
        class(vector_t), intent(in) :: this
        type(vector_t)              :: w
        w % x = this % x + a
        w % y = this % y + a
    end function scalar_add


    elemental function add_scalar(this,a) result(w)
        class(vector_t), intent(in) :: this
        real(wp), intent(in)        :: a
        type(vector_t)              :: w
        w % x = this % x + a
        w % y = this % y + a
    end function add_scalar


    elemental function sub_vector(this,v) result(w)
        class(vector_t), intent(in) :: this
        type(vector_t),  intent(in) :: v 
        type(vector_t)              :: w
        w % x = this % x - v % x
        w % y = this % y - v % y
    end function sub_vector


    elemental function scalar_sub(a, this) result(w)
        real(wp), intent(in)        :: a
        class(vector_t), intent(in) :: this
        type(vector_t)              :: w
        w % x = a - this % x
        w % y = a - this % y 
    end function scalar_sub


    elemental function sub_scalar(this,a) result(w)
        class(vector_t), intent(in) :: this
        real(wp), intent(in)        :: a
        type(vector_t)              :: w
        w % x = this % x - a
        w % y = this % y - a
    end function sub_scalar


    elemental function mul_vector(this,v) result(w)
        class(vector_t), intent(in) :: this
        type(vector_t),  intent(in) :: v 
        type(vector_t)              :: w
        w % x = this % x * v % x
        w % y = this % y * v % y
    end function mul_vector


    elemental function scalar_mul(a, this) result(w)
        real(wp), intent(in)        :: a
        class(vector_t), intent(in) :: this
        type(vector_t)              :: w
        w % x = this % x * a
        w % y = this % y * a
    end function scalar_mul


    elemental function mul_scalar(this,a) result(w)
        class(vector_t), intent(in) :: this
        real(wp), intent(in)        :: a
        type(vector_t)              :: w
        w % x = this % x * a
        w % y = this % y * a
    end function mul_scalar


    elemental function div_vector(this,v) result(w)
        class(vector_t), intent(in) :: this
        type(vector_t),  intent(in) :: v 
        type(vector_t)              :: w
        w % x = this % x / v % x
        w % y = this % y / v % y
    end function div_vector


    elemental function div_scalar(this,d) result(w)
        class(vector_t), intent(in) :: this
        real(wp), intent(in)        :: d
        type(vector_t)              :: w
        w % x = this % x / d
        w % y = this % y / d
    end function div_scalar


    pure subroutine assign_vector(this,v) 
        class(vector_t), intent(inout) :: this
        type(vector_t),  intent(in)    :: v 
        this % x = v % x
        this % y = v % y
    end subroutine assign_vector


    pure subroutine assign_scalar(this,a)
        class(vector_t), intent(inout) :: this
        real(wp), intent(in)           :: a
        this % x = a
        this % y = a
    end subroutine assign_scalar


    elemental function dot(u,v) result(val)
        type(vector_t), intent(in) :: u
        type(vector_t), intent(in) :: v
        real(wp)                   :: val 
        val = (u % x)*(v % x) + (u % y)*(v % y)
    end function dot


    elemental function mag(u) result(val)
        type(vector_t), intent(in) :: u
        real(wp)                   :: val 
        val = sqrt(dot(u,u))
    end function mag


end module slbm_2d_vector
