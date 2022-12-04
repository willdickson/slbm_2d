module slbm_2d_const

    use slbm_2d_kinds,    only : wp, ip
    use slbm_2d_vector,   only : vector_t
    implicit none
    private

    public :: PI                   ! mathematical constant (3.145...)

    ! Body type ids
    public :: BODY_TYPE_OPEN       ! Open curve body type
    public :: BODY_TYPE_CLOSED     ! Closed curve body type
    public :: BODY_TYPE_UNKNOWN    ! Unknown body type

    public :: BODY_PT_MAX_NBRS     ! Max number (mesh) neighbors for body pt.

    ! Stop condition ids
    public :: STOP_COND_TIME       ! Run time (duration) stop condition 
    public :: STOP_COND_STEADY     ! Steady state condition 
    public :: STOP_COND_UNKNOWN    ! Unknown stop condition (error) 

    ! Boundary information
    public :: NUM_BNDRY            ! Number of boundries for rect. region
    public :: BNDRY_SIDE_LEFT      ! left boundary id
    public :: BNDRY_SIDE_RIGHT     ! right boundary id
    public :: BNDRY_SIDE_TOP       ! top boundary id
    public :: BNDRY_SIDE_BOTTOM    ! bottom boundary id
    public :: BNDRY_SIDE_UNKNOWN   ! unknwn boundary id
    public :: BNDRY_SIDE_NAMES     ! Array of boundary names
    public :: BNDRY_SIDE_IDS       ! Array of boundary names

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

    ! Mesh constants
    public :: MESH_LEN_ETOL        ! Relative tolerance for mesh when set
                                   ! using len_x and len_y


    ! ------------------------------------------------------------------
    real(wp),    parameter :: PI  = 4_wp*atan(1.0_wp)

    ! Body type ids
    integer(ip), parameter :: BODY_TYPE_OPEN    = 1 
    integer(ip), parameter :: BODY_TYPE_CLOSED  = 2 
    integer(ip), parameter :: BODY_TYPE_UNKNOWN = 3

    integer(ip), parameter :: BODY_PT_MAX_NBRS  = 25   

    ! Stop condition ids
    integer(ip), parameter :: STOP_COND_TIME    = 1
    integer(ip), parameter :: STOP_COND_STEADY  = 2
    integer(ip), parameter :: STOP_COND_UNKNOWN = 3

    ! Boundry number and names
    integer(ip), parameter :: NUM_BNDRY = 4
    integer(ip), parameter :: BNDRY_SIDE_LEFT    = 1
    integer(ip), parameter :: BNDRY_SIDE_RIGHT   = 2
    integer(ip), parameter :: BNDRY_SIDE_TOP     = 3
    integer(ip), parameter :: BNDRY_SIDE_BOTTOM  = 4
    integer(ip), parameter :: BNDRY_SIDE_UNKNOWN = 5
    character(len=10)      :: BNDRY_SIDE_NAMES(NUM_BNDRY) = & 
        [character(len=10) :: 'left', 'right', 'top', 'bottom']
    integer(ip), parameter :: BNDRY_SIDE_IDS(NUM_BNDRY) = [ & 
        BNDRY_SIDE_LEFT,    & 
        BNDRY_SIDE_RIGHT,   &
        BNDRY_SIDE_TOP,     &
        BNDRY_SIDE_BOTTOM   &
        ]

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
    real(wp),    parameter :: CS4 = CS2**2

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

    real(wp), parameter    :: MESH_LEN_ETOL = 1.0e-6

end module slbm_2d_const
