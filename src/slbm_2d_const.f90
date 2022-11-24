module slbm_2d_const

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_vector,   only : vector_t
    implicit none
    private

    public :: PI                   ! mathematical constant (3.145...)

    ! Stop condition ids
    public :: STOP_COND_TIME       ! Run time (duration) stop condition 
    public :: STOP_COND_STEADY     ! Steady state condition 
    public :: STOP_COND_UNKNOWN    ! Unknown stop condition (error) 

    ! Boundary information
    public :: NUM_BNDRY            ! Number of boundries for rect. region
    public :: BNDRY_NAMES          ! Array of boundary names

    ! Boundary condition ids       
    public :: BNDRY_COND_INFLOW    ! Inflow boundary condition
    public :: BNDRY_COND_MOVING    ! Moving wall boundary condition
    public :: BNDRY_COND_OUTFLOW   ! Outflow boundary condition
    public :: BNDRY_COND_NOSLIP    ! No slip boundary condition
    public :: BNDRY_COND_SLIP      ! Slip boundary condition
    public :: BNDRY_COND_UNKNOWN   ! Unknown boundary condition (error)

    ! Initial condition ids
    public :: INIT_COND_CONST      ! Constant velocity initial condition
    public :: INIT_COND_FILE       ! Initial condition from file
    public :: INIT_COND_UNKNOWN    ! Unknown initial condition

    ! Latice constants
    public :: CS                   ! lbm speed of sound 
    public :: CS2                  ! lbm square of sound speed
    public :: CS4                  ! lbm 4th power of sound speed
    public :: LATTICE_Q            ! Number of lattice velocity vectors
    public :: LATTICE_W            ! Lattice velocity vector weights
    public :: LATTICE_E            ! Lettice velocity vectors

    ! ------------------------------------------------------------------
    real(wp),    parameter :: PI  = 4_wp*atan(1.0_wp)

    ! Stop condition ids
    integer(ip), parameter :: STOP_COND_TIME    = 1
    integer(ip), parameter :: STOP_COND_STEADY  = 2
    integer(ip), parameter :: STOP_COND_UNKNOWN = 3

    ! Boundry number and names
    integer(ip), parameter :: NUM_BNDRY = 4
    character(len=10)      :: BNDRY_NAMES(NUM_BNDRY) = & 
        [character(len=10) :: 'left', 'right', 'top', 'bottom']

    ! Boundary condition ids
    integer(ip), parameter :: BNDRY_COND_INFLOW   = 1
    integer(ip), parameter :: BNDRY_COND_MOVING   = 2
    integer(ip), parameter :: BNDRY_COND_OUTFLOW  = 3
    integer(ip), parameter :: BNDRY_COND_NOSLIP   = 4
    integer(ip), parameter :: BNDRY_COND_SLIP     = 5
    integer(ip), parameter :: BNDRY_COND_UNKNOWN  = 6

    ! Initial condition ids
    integer(ip), parameter :: INIT_COND_CONST    = 1
    integer(ip), parameter :: INIT_COND_FILE     = 2
    integer(ip), parameter :: INIT_COND_UNKNOWN  = 3

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
