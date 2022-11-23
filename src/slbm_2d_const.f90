module slbm_2d_const

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_vector,   only : vector_t
    implicit none
    private

    public :: PI                ! mathematical constant (3.145...)

    ! Stop condition ids
    public :: STOP_COND_TIME    ! Run time (duration) stop condition 
    public :: STOP_COND_STEADY  ! Steady state condition 
    public :: STOP_COND_UNKNOWN ! Steady state condition 

    ! Boundary condition ids
    public :: BNDRY_COND_INFLOW
    public :: BNDRY_COND_OUTFLOW
    public :: BNDRY_COND_NOSLIP
    public :: BNDRY_COND_SLIP
    public :: BNDRY_COND_UNKNOWN

    ! Latice constants
    public :: CS                ! lbm speed of sound 
    public :: CS2               ! lbm square of sound speed
    public :: CS4               ! lbm 4th power of sound speed
    public :: LATTICE_Q         ! Number of lattice velocity vectors
    public :: LATTICE_W         ! Lattice velocity vector weights
    public :: LATTICE_E         ! Lettice velocity vectors

    ! ------------------------------------------------------------------
    real(wp),    parameter :: PI  = 4_wp*atan(1.0_wp)

    ! Stop condition ids
    integer(ip), parameter :: STOP_COND_TIME    = 1
    integer(ip), parameter :: STOP_COND_STEADY  = 2
    integer(ip), parameter :: STOP_COND_UNKNOWN = 3

    ! Boundary condition ids
    integer(ip), parameter :: BNDRY_COND_INFLOW   = 1
    integer(ip), parameter :: BNDRY_COND_OUTFLOW  = 2
    integer(ip), parameter :: BNDRY_COND_NOSLIP   = 3
    integer(ip), parameter :: BNDRY_COND_SLIP     = 4
    integer(ip), parameter :: BNDRY_COND_UNKNOWN  = 1
    

    ! Latice constants
    real(wp),    parameter :: CS  = 1.0_wp/sqrt(3.0_wp)
    real(wp),    parameter :: CS2 = 1.0_wp/3.0_wp
    real(wp),    parameter :: CS4 = CS2*CS2

    integer(ip), parameter :: LATTICE_Q = 9
    real(wp),    parameter :: LATTICE_W(LATTICE_Q) = [ &
        4.0_wp / 9.0_wp,   & 
        1.0_wp / 9.0_wp,   &
        1.0_wp / 9.0_wp,   &
        1.0_wp / 9.0_wp,   & 
        1.0_wp / 9.0_wp,   &
        1.0_wp / 36.0_wp,  &
        1.0_wp / 36.0_wp,  &
        1.0_wp / 36.0_wp,  &
        1.0_wp / 36.0_wp   &
        ]

    type(vector_t), parameter    :: LATTICE_E(LATTICE_Q) = [ &
        vector_t( 0.0_wp,  0.0_wp),  &
        vector_t( 1.0_wp,  0.0_wp),  &
        vector_t( 0.0_wp,  1.0_wp),  &
        vector_t(-1.0_wp,  0.0_wp),  &
        vector_t( 0.0_wp, -1.0_wp),  &
        vector_t( 1.0_wp,  1.0_wp),  &
        vector_t(-1.0_wp,  1.0_wp),  &
        vector_t(-1.0_wp, -1.0_wp),  &
        vector_t( 1.0_wp, -1.0_wp)   &
        ]

end module slbm_2d_const
