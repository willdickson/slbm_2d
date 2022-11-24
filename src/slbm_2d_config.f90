module slbm_2d_config

    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_value 
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_const,  only : NUM_BNDRY
    use slbm_2d_const,  only : BNDRY_NAMES
    use slbm_2d_const,  only : STOP_COND_TIME
    use slbm_2d_const,  only : STOP_COND_STEADY
    use slbm_2d_const,  only : STOP_COND_UNKNOWN
    use slbm_2d_const,  only : BNDRY_COND_INFLOW
    use slbm_2d_const,  only : BNDRY_COND_MOVING
    use slbm_2d_const,  only : BNDRY_COND_OUTFLOW
    use slbm_2d_const,  only : BNDRY_COND_NOSLIP
    use slbm_2d_const,  only : BNDRY_COND_SLIP
    use slbm_2d_const,  only : BNDRY_COND_UNKNOWN
    use slbm_2d_const,  only : INIT_COND_CONST
    use slbm_2d_const,  only : INIT_COND_FILE
    use slbm_2d_const,  only : INIT_COND_UNKNOWN
    use slbm_2d_const,  only : CS
    use slbm_2d_const,  only : CS2
    use slbm_2d_bndry,  only : bndry_t
    use slbm_2d_bndry,  only : bndry_ptr_t
    use slbm_2d_vector, only : vector_t
    use slbm_2d_init,   only : init_t
    use slbm_2d_init,   only : init_const_t
    use tomlf, only : toml_table
    use tomlf, only : toml_error
    use tomlf, only : toml_parse
    use tomlf, only : get_value

    implicit none
    private

    type, public :: config_t

        ! Grid parameters
        integer(ip) :: num_x      = 0_ip   ! number of x grid points
        integer(ip) :: num_y      = 0_ip   ! number of y grid points
        real(wp)    :: ds         = 0.0_wp ! mesh step size  ds=dx=dy

        ! Fluid parameters
        real(wp)    :: kvisc      = 0.0_wp ! kinematic viscosity 
        real(wp)    :: density    = 0.0_wp ! reference density


        ! Stopping conditon parameters
        real(wp)    :: stop_time  = 0.0_wp ! simulation stop time
        real(wp)    :: stop_etol  = 0.0_wp ! simulation (steady) stop err tol.
        integer(ip) :: stop_cond  = STOP_COND_TIME 

        ! Boundary condition parameters
        type(bndry_t) :: bndry_left         ! left   boundary condition
        type(bndry_t) :: bndry_right        ! right  boundary condition
        type(bndry_t) :: bndry_top          ! top    boundary condition 
        type(bndry_t) :: bndry_bottom       ! bottom boundary condition

        ! Initial condition parameters
        class(init_t), pointer  :: init     ! initial condition

        ! Derived configuration parameters
        integer(ip) :: nstep   = 0_ip   ! number of time steps
        real(wp)    :: len_x   = 0.0_wp ! length of mesh in x dim
        real(wp)    :: len_y   = 0.0_wp ! length of mesh in y dim
        real(wp)    :: dt      = 0.0_wp ! time step
        real(wp)    :: tau     = 0.0_wp ! relaxation parameter

    contains
        private
        procedure, public :: pprint => config_pprint
    end type config_t 

    interface config_t
        procedure config_from_toml
        procedure config_constructor
    end interface config_t


contains

    ! config_t constructors
    ! -------------------------------------------------------------------------

    function config_from_toml(filename) result(config)
        character(*),     intent(in)          :: filename ! name of config file
        type(config_t)                        :: config   
        type(toml_table), allocatable, target :: table
        type(toml_table), pointer             :: table_ptr

        ! Check to make sure the configuration file exists
        block
            logical :: file_exists
            inquire(file=filename, exist=file_exists)
            if ( .not. file_exists ) then
                print *, "unable to continue - " // filename // " doesn't exist"
                stop 
            end if
        end block

        ! Parse the configuration file - exit on failure. 
        block 
            integer                       :: io
            type(toml_error), allocatable :: parse_error
            open(newunit=io, file=filename)
            call toml_parse(table, io, parse_error)
            close(unit=io)
            ! If there was an error parsing the file exit
            if ( allocated(parse_error) ) then 
                print *, "unable to parse file " // filename
                print *, parse_error % message
                stop
            end if
            table_ptr => table
        end block

        !Read values from  sub-tables
        call read_mesh_table(config, table_ptr)
        call read_fluid_table(config, table_ptr)
        call read_init_table(config, table_ptr)
        call read_bndry_table(config, table_ptr)
        call read_stop_table(config, table_ptr)

        ! Set derivied values
        config % dt = config % ds 
        config % nstep = nint(config % stop_time / (config % dt))
        config % len_x = (config % num_x - 1.0_wp) * config % ds 
        config % len_y = (config % num_y - 1.0_wp) * config % ds
        config % tau = 0.5_wp + (config % kvisc) / (CS2 * config % dt)

    end function config_from_toml

    function config_constructor( &
            num_x,        &  ! number of z grid points
            num_y,        &  ! number of y grid points
            ds,           &  ! mesh step size
            kvisc,        &  ! kinematic viscosity
            density,      &  ! reference density
            stop_time,    &  ! simulation stop time  
            stop_etol,    &  ! error tolerance for steady stop
            stop_cond,    &  ! stop condition id
            bndry_left,   &  ! left boundary condition  (optional)
            bndry_right,  &  ! right boundary condition (optional)
            bndry_top,    &  ! top boundary condition   (optional)
            bndry_bottom, &  ! bottom boundary contiion (optional)
            init          &  ! initial conditions
            ) result(config)
        integer(ip),   intent(in)                    :: num_x
        integer(ip),   intent(in)                    :: num_y
        real(wp),      intent(in)                    :: ds
        real(wp),      intent(in)                    :: kvisc
        real(wp),      intent(in), optional          :: density
        real(wp),      intent(in), optional          :: stop_time
        real(wp),      intent(in), optional          :: stop_etol
        integer(ip),   intent(in), optional          :: stop_cond
        type(bndry_t), intent(in), optional          :: bndry_left
        type(bndry_t), intent(in), optional          :: bndry_right
        type(bndry_t), intent(in), optional          :: bndry_top
        type(bndry_t), intent(in), optional          :: bndry_bottom
        class(init_t), intent(in), optional, pointer :: init
        type(config_t)                               :: config
        type(bndry_t)                                :: noslip_bndry
        type(init_const_t), pointer                  :: init_const_zero

        ! Set input parameters
        config % num_x = abs(num_x)
        config % num_y = abs(num_y)
        config % ds = abs(ds)
        config % kvisc = abs(kvisc)
        config % stop_time = abs(stop_time)

        ! Set derivied values
        config % dt = ds
        config % nstep = nint(stop_time / (config % dt))
        config % len_x = (num_x - 1.0_wp) * ds 
        config % len_y = (num_y - 1.0_wp) * ds
        config % tau = 0.5_wp + kvisc / (CS2 * config % dt)

        ! Set default values
        config % density = 1.0_wp
        config % stop_time = 1.0_wp
        config % stop_etol = 1.0e-6_wp
        config % stop_cond = STOP_COND_TIME

        noslip_bndry % id = BNDRY_COND_NOSLIP
        noslip_bndry % velocity = vector_t(0.0_wp, 0.0_wp)

        config % bndry_left = noslip_bndry
        config % bndry_right = noslip_bndry
        config % bndry_top = noslip_bndry
        config % bndry_bottom = noslip_bndry

        allocate(init_const_zero)
        init_const_zero % id = INIT_COND_CONST
        init_const_zero % velocity = vector_t(0.0_wp, 0.0_wp)
        config % init => init_const_zero

        ! Set optional values
        if ( present(density) ) then
            config % density = abs(density)
        end if
        if ( present(stop_time) ) then
            config % stop_time = abs(stop_time)
        end if
        if ( present(stop_etol) ) then
            config % stop_etol = abs(stop_etol)
        end if
        if ( present(stop_cond) ) then
            block
                logical :: ok = .false.
                select case(stop_cond)
                case (STOP_COND_TIME, STOP_COND_STEADY)
                    ok = .true.
                end select
            end block
            config % stop_cond = stop_cond
        end if
        ! bndrys are NOT CHECKED ... do this in bndry constructor?
        if ( present(bndry_left) ) then
            config % bndry_left = bndry_left
        end if
        if ( present(bndry_right) ) then
            config % bndry_right = bndry_right
        end if
        if ( present(bndry_top) ) then
            config % bndry_top = bndry_top
        end if
        if ( present(bndry_bottom) ) then
            config % bndry_bottom = bndry_bottom
        end if

        ! NOT CHECKED .... do this in init constructor?
        if ( present(init) ) then
            config % init => init
        end if


    end function config_constructor


    !! Initial condition parameters
    !class(init_t), pointer  :: init     ! initial condition


    subroutine config_pprint(this)
        class(config_t), intent(in), target :: this
        print *, ''
        print *, 'grid parameters'
        print *, '---------------------------------------------'
        print *, 'num_x     = ', this % num_x
        print *, 'num_y     = ', this % num_y
        print *, 'ds        = ', this % ds
        print *, ''
        print *, 'fluid parameters'
        print *, '---------------------------------------------'
        print *, 'kvisc     = ', this % kvisc
        print *, 'density   = ', this % density
        print *, ''
        print *, 'stop parameters'
        print *, '---------------------------------------------'
        print *, 'stop_cond = ', stop_cond_to_string(this % stop_cond)
        print *, 'stop_time = ', this % stop_time
        if ( this % stop_cond == STOP_COND_STEADY ) then
            print *, 'stop_etol = ', this % stop_etol
        end if
        print *, ''
        print *, 'boundry parameters'
        print *, '---------------------------------------------'
        block
            integer(ip)               :: i
            character(:), allocatable :: bndry_name
            character(:), allocatable :: type_string
            type(bndry_ptr_t)         :: bndry_array(NUM_BNDRY)
            bndry_array(1) % ptr => this % bndry_left
            bndry_array(2) % ptr => this % bndry_right
            bndry_array(3) % ptr => this % bndry_top
            bndry_array(4) % ptr => this % bndry_bottom
            do i=1,NUM_BNDRY
                bndry_name = 'bndry.'//trim(BNDRY_NAMES(i))
                type_string = bndry_to_string(bndry_array(i) % ptr % id)
                print *, [character(15)::bndry_name]// 'type', ' = ', type_string 
                print *, [character(16)::bndry_name//'.id'], & 
                    '=', bndry_array(i) % ptr % id
                print *, [character(24)::bndry_name//'.velocity.x'], & 
                    '=',  bndry_array(i) % ptr % velocity % x
                print *, [character(24)::bndry_name//'.velocity.y'], & 
                    '=',  bndry_array(i) % ptr % velocity % y
                print *, ''
            end do
        end block
        print *, 'initial condition parameters'
        print *, '---------------------------------------------'
        print *, 'name      = ', this % init % name
        print *, 'id        = ', this % init % id
        select type (init => this % init)
        class is (init_const_t)
            print *, 'velcity x = ', init % velocity % x
            print *, 'velcity y = ', init % velocity % y
        class default
            print *, 'value display not implemented'
        end select
        print *, ''
        print *, 'derived parameters'
        print *, '---------------------------------------------'
        print *, 'nstep     = ', this % nstep
        print *, 'len_x     = ', this % len_x
        print *, 'len_y     = ', this % len_y
        print *, 'dt        = ', this % dt
        print *, 'tau       = ', this % tau
        print *, ''
    end subroutine config_pprint


    ! utiliy functions & subroutines
    ! -------------------------------------------------------------------------

    subroutine read_mesh_table(config, table)
        type(config_t), intent(inout)          :: config
        type(toml_table), pointer, intent(in)  :: table
        type(toml_table), pointer              :: mesh_table
        integer(ip)                            :: stat
        real(wp)                               :: nan
        nan = ieee_value(nan, ieee_quiet_nan)

        call get_value(table, "mesh", mesh_table, .false.)
        if (.not. associated(mesh_table) ) then
            print *, "config .toml missing mesh section"
            stop
        endif

        config % ds = nan
        call get_value(mesh_table, "ds", config % ds, nan, stat)
        if ( (stat < 0) .or. ( ieee_is_nan(config % ds) ) ) then 
            print *, "config .toml mesh is missing ds or value is not a real" 
            print *, "ds = ", config % ds
            stop
        end if
        config % ds = abs(config % ds)

        config % num_x = 0_ip
        call get_value(mesh_table, "num_x", config % num_x, -1_ip, stat)
        if ( (stat < 0) .or. (config % num_x < 0) ) then
            print *, "config .toml mesh is missing num_x or value is negative"
            print *, "num_x = ", config % num_x
            stop
        end if

        config % num_y = 0_ip
        call get_value(mesh_table, "num_y", config % num_y, -1_ip, stat)
        if ( (stat < 0) .or. (config % num_y < 0) ) then
            print *, "config .toml mesh is missing num_y or value is negative"
            print *, "num_y = ", config % num_y
            stop
        end if
    end subroutine read_mesh_table


    subroutine read_fluid_table(config, table)
        type(config_t), intent(inout)          :: config
        type(toml_table), pointer, intent(in)  :: table
        type(toml_table), pointer              :: fluid_table
        integer(ip)                            :: stat
        real(wp)                               :: nan
        nan = ieee_value(nan, ieee_quiet_nan)

        call get_value(table, "fluid", fluid_table, .false.)
        if (.not. associated(fluid_table) ) then
            print *, "config .toml missing fluid section"
            stop
        endif

        config % kvisc = nan
        call get_value(fluid_table, "kvisc", config % kvisc, nan, stat)
        if ( (stat < 0) .or. (ieee_is_nan(config % kvisc)) ) then
            print *, "config .toml fluid is missing kvisc or is not a real"
            print *, "kvisc = ", config % kvisc
            stop
        end if
        config % kvisc = abs(config % kvisc)

        config % density = nan
        call get_value(fluid_table, "density", config % density, nan, stat)
        if ( (stat < 0) .or. (ieee_is_nan(config % density)) ) then
            print *, "config .toml fluid is missing density or is not a real"
            print *, "kvisc = ", config % density 
            stop
        end if
        config % density = abs(config % density)
    end subroutine read_fluid_table


    subroutine read_init_table(config, table)
        type(config_t), intent(inout)          :: config
        type(toml_table), pointer, intent(in)  :: table
        type(toml_table), pointer              :: init_table
        integer(ip)                            :: stat
        integer(ip)                            :: init_id
        character(:), allocatable              :: init_name

        call get_value(table, "init", init_table, .false.)
        if (.not. associated(init_table) ) then
            print *, "config .toml missing init section"
            stop
        endif

        call get_value(init_table, "type", init_name, init_name) 
        if ( .not. allocated(init_name) ) then
            print *, "config .toml init is missing type or it is not a string"
            stop
        end if

        select case(init_name)
        case ('constant') 
            block
                type(init_const_t), pointer :: init
                call create_init_const(init_table, init_name, init)
                config % init => init
            end block
        case ('file')
            print *, "config .toml init type file not implement yet"
            stop
        case default
            print *, "config .toml init unknown type"
            print *, "type = ", init_name
            stop
        end select

    end subroutine read_init_table

    subroutine create_init_const(init_table, init_name, init) 
        type(toml_table), pointer, intent(in)    :: init_table
        character(*), intent(in)                 :: init_name
        type(init_const_t), pointer, intent(out) :: init
        type(toml_table), pointer                :: velo_table
        integer(ip)                              :: stat
        real(wp)                                 :: nan
        nan = ieee_value(nan, ieee_quiet_nan)

        allocate(init)
        init % name = init_name
        init % id = INIT_COND_CONST

        call get_value(init_table, "velocity", velo_table, .false.)
        if (.not. associated(velo_table) ) then
            print *, "config .toml init, type=constant missing velocity"
            stop
        endif

        init % velocity % x = nan
        call get_value(velo_table, "x", init % velocity % x, nan, stat)
        if ( (stat < 0) .or. (ieee_is_nan(init % velocity % x)) ) then
            print *, "config .toml init velocity  missing x or is not a real"
            print *, "x = ", init % velocity % x
            stop
        end if

        init % velocity % y = nan
        call get_value(velo_table, "y", init % velocity % y, nan, stat)
        if ( (stat < 0) .or. (ieee_is_nan(init % velocity % y)) ) then
            print *, "config .toml init velocity  missing y or is not a real"
            print *, "y = ", init % velocity % y
            stop
        end if

    end subroutine create_init_const


    subroutine read_stop_table(config, table)
        type(config_t), intent(inout)          :: config
        type(toml_table), pointer, intent(in)  :: table
        type(toml_table), pointer              :: stop_table
        character(:), allocatable              :: cond
        integer(ip)                            :: stat
        real(wp)                               :: nan
        nan = ieee_value(nan, ieee_quiet_nan)

        call get_value(table, "stop", stop_table, .false.)
        if (.not. associated(stop_table) ) then
            print *, "config .toml missing stop section"
            stop
        endif

        call get_value(stop_table, "cond", cond, cond) 
        if ( .not. allocated(cond) ) then
            print *, "config .toml stop is missing cond or not a string"
            stop
        end if
        config % stop_cond = stop_cond_from_string(cond)
        if (config % stop_cond == STOP_COND_UNKNOWN) then
            print *, "config .toml stop, unknown stop condition"
            print *, "cond = ", cond
            stop
        end if

        config % stop_time = nan
        call get_value(stop_table, "time", config % stop_time, nan, stat)
        if ( (stat < 0) .or. (ieee_is_nan(config % stop_time)) ) then
            print *, "config .toml stop is missing time or is not a real"
            print *, "stop_time = ", config % stop_time
            stop
        end if
        config % stop_time = abs(config % stop_time)

        config % stop_etol = nan
        call get_value(stop_table, "etol", config % stop_etol, nan, stat)
        if ( (stat < 0) .or. (ieee_is_nan(config % stop_etol)) ) then
            print *, "config .toml stop is missing time or is not a real"
            print *, "stop_etol = ", config % stop_etol
            stop
        end if
        config % stop_etol = abs(config % stop_etol)
    end subroutine read_stop_table


    subroutine read_bndry_table(config, table)
        type(config_t), target, intent(inout)  :: config
        type(toml_table), pointer, intent(in)  :: table
        type(toml_table), pointer              :: bndry_table
        type(toml_table), pointer              :: side_table
        type(bndry_ptr_t)                      :: bndry_array(NUM_BNDRY)
        character(:), allocatable              :: type_name
        integer(ip)                            :: type_id
        integer(ip)                            :: stat
        integer(ip)                            :: i
        real(wp)                               :: nan
        real(wp)                               :: val 
        nan = ieee_value(nan, ieee_quiet_nan)

        call get_value(table, "bndry", bndry_table, .false.)
        if (.not. associated(bndry_table) ) then
            print *, "config .toml missing bndry section"
            stop
        endif

        bndry_array(1) % ptr => config % bndry_left
        bndry_array(2) % ptr => config % bndry_right
        bndry_array(3) % ptr => config % bndry_top
        bndry_array(4) % ptr => config % bndry_bottom

        do i=1,NUM_BNDRY
            call get_value(bndry_table, BNDRY_NAMES(i), side_table, .false.)
            if (.not. associated(side_table) ) then
                print *, "config .toml is missing boundary." & 
                    // trim(BNDRY_NAMES(i))
                stop
            endif

            call get_value(side_table, "type", type_name, type_name, stat)
            if ( .not. allocated(type_name) ) then
                print *, "config .toml boundry."//trim(BNDRY_NAMES(i))& 
                    // "type missing"
                stop
            end if
            type_id = bndry_from_string(type_name)
            if (type_id == BNDRY_COND_UNKNOWN) then
                print *, "config .toml boundry."//trim(BNDRY_NAMES(i))& 
                    // " type unknown condition"
                print *, "type = ", type_name
                stop
            end if
            bndry_array(i) % ptr % id = type_id

            select case (type_id)
            case (BNDRY_COND_INFLOW, BNDRY_COND_MOVING)
                val = nan
                call get_value(side_table, "value", val, nan, stat) 
                if ( (stat < 0) .or. (ieee_is_nan(val)) ) then
                    print *, "config .toml boundary."//trim(BNDRY_NAMES(i))& 
                        // " is missing "//type_name//" value or is not a real"
                    print *, "value = ", val 
                    stop
                end if
            case default
                val = 0.0_wp
            end select 

            if ((type_id == BNDRY_COND_INFLOW) .and. (val < 0.0_wp)) then
                print *, "config .toml boundary."//trim(BNDRY_NAMES(i))&
                    // " value is < 0"
                stop
            end if

            ! Set boundary velocity values
            block
                type(vector_t) :: v
                v = get_bndry_velocity(BNDRY_NAMES(i), type_id, val)
                bndry_array(i) % ptr % velocity = v 
            end block
        end do

    end subroutine read_bndry_table


    function stop_cond_from_string(stop_cond_string) result(cstop)
        character(*), intent(in) :: stop_cond_string
        integer(ip)              :: cstop
        select case (stop_cond_string)
        case ('time')
            cstop = STOP_COND_TIME
        case ('steady')
            cstop = STOP_COND_STEADY
        case default
            cstop = STOP_COND_UNKNOWN
        end select
    end function stop_cond_from_string


    function stop_cond_to_string(cstop) result(stop_cond_string)
        integer(ip), intent(in)   :: cstop
        character(:), allocatable :: stop_cond_string
        select case (cstop)
        case (STOP_COND_TIME)
            stop_cond_string = "time"
        case (STOP_COND_STEADY)
            stop_cond_string = "steady"
        case default
            stop_cond_string = "unknown"
        end select
    end function stop_cond_to_string


    function bndry_from_string(type_string) result(id)
        character(*), intent(in) :: type_string
        integer(ip)              :: id
        select case(type_string)
        case ("inflow")
            id = BNDRY_COND_INFLOW
        case ("moving")
            id = BNDRY_COND_MOVING
        case ("outflow")
            id = BNDRY_COND_OUTFLOW
        case ("noslip")
            id = BNDRY_COND_NOSLIP
        case ("slip")
            id = BNDRY_COND_SLIP
        case default
            id = BNDRY_COND_UNKNOWN
        end select
    end function bndry_from_string


    function bndry_to_string(type_id) result(type_string)
        integer(ip), intent(in)   :: type_id
        character(:), allocatable :: type_string
        select case(type_id)
        case (BNDRY_COND_INFLOW)
            type_string = "inflow"
        case (BNDRY_COND_MOVING)
            type_string = "moving"
        case (BNDRY_COND_OUTFLOW)
            type_string = "outflow"
        case (BNDRY_COND_NOSLIP)
            type_string = "noslip"
        case (BNDRY_COND_SLIP)
            type_string = "slip"
        case default
            type_string = "unknown"
        end select
    end function bndry_to_string


    function get_bndry_velocity(bndry_name, type_id, val) result(velocity)
        character(*), intent(in) :: bndry_name
        integer(ip), intent(in)  :: type_id
        real(wp), intent(in)     :: val
        type(vector_t)           :: velocity 
        velocity = vector_t(0.0_wp, 0.0_wp)
        select case (type_id)
        case (BNDRY_COND_INFLOW)
            select case(bndry_name)
            case ("left")
                velocity % x = abs(val)
                velocity % y = 0.0_wp
            case ("right")
                velocity % x = -abs(val)
                velocity % y = 0.0_wp
            case ("top")
                velocity % x = 0.0_wp 
                velocity % y = -abs(val)
            case ("bottom")
                velocity % x = 0.0_wp  
                velocity % y = abs(val)
            end select
        case (BNDRY_COND_MOVING)
            select case(bndry_name)
            case ("left", "right")
                velocity % x = 0.0_wp 
                velocity % y = val 
            case ("top", "bottom")
                velocity % x = val 
                velocity % y = 0.0_wp 
            end select
        end select
    end function get_bndry_velocity


    function init_id_from_string(init_string) result(init_id)
        character(*), intent(in) :: init_string
        integer(ip)              :: init_id
        select case (init_string)
        case ('constant')
            init_id = INIT_COND_CONST
        case ('file')
            init_id = INIT_COND_FILE
        case default
            init_id = INIT_COND_UNKNOWN
        end select
    end function init_id_from_string


    function init_id_to_string(init_id) result(init_string)
        integer(ip), intent(in)   :: init_id
        character(:), allocatable :: init_string
        select case (init_id)
        case (INIT_COND_CONST)
            init_string = 'constant'
        case (INIT_COND_FILE)
            init_string = 'file'
        case default
            init_string = 'unknown'
        end select
    end function init_id_to_string


end module slbm_2d_config
