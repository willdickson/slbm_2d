module slbm_2d_funcs

    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_const,  only : PI
    use slbm_2d_const,  only : CS2
    use slbm_2d_const,  only : CS4
    use slbm_2d_const,  only : LATTICE_E
    use slbm_2d_const,  only : LATTICE_W
    use slbm_2d_vector, only : vector_t
    use slbm_2d_vector, only : dot

    implicit none
    private

    public :: kernel
    public :: equilib_func

    ! Constants used in equilibrium function calculation
    real(wp), parameter  :: A1 = 1.0_wp/CS2
    real(wp), parameter  :: A2 = 1.0_wp/(2.0_wp*CS4)
    real(wp), parameter  :: A3 = 1.0_wp/(2.0_wp*CS2)

contains


    function equilib_func(rho, u, k) result(feq)
        real(wp),       intent(in)   :: rho ! fluid density
        type(vector_t), intent(in)   :: u   ! fluid velocity vector
        integer(ip),    intent(in)   :: k   ! lattice index
        real(wp)                     :: feq ! kth comp. of equilibrum dist.

        real(wp)   :: wt    ! k-th lattice weight
        real(wp)   :: uu    ! squared magnitude of velocity
        real(wp)   :: eu    ! lattice velocity (ex,ey) to velocity (ux,uy)
        real(wp)   :: eu2   ! square of eu

        wt  = LATTICE_W(k)
        uu  = dot(u,u)
        eu  = dot(LATTICE_E(k),u)
        eu2 = eu**2 
        feq = rho*wt*(1.0_wp + A1*eu + A2*eu2 - A3*uu)
    end function equilib_func


    function kernel(p, q, ds) result(val)
        type(vector_t), intent(in) :: p    
        type(vector_t), intent(in) :: q
        real(wp), intent(in)       :: ds
        real(wp)                   :: val
        real(wp)                   :: dx
        real(wp)                   :: dy 
        real(wp)                   :: kx 
        real(wp)                   :: ky 
        dx = abs(p % x - q % x)/ds
        dy = abs(p % y - q % y)/ds
        if ((dx >= 2.0_wp .or. dy >= 2.0_wp)) then
            val = 0.0_wp
        else
            kx  = 0.25_wp*(1.0_wp + cos(0.5_wp*PI*dx))
            ky  = 0.25_wp*(1.0_wp + cos(0.5_wp*PI*dy))
            val = (kx * ky) / (ds**2)
        end if
    end function kernel


end module slbm_2d_funcs
