module slbm_2d_config
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: ieee_arithmetic, only : ieee_value 
    use slbm_2d_kinds, only : wp, ip
    use slbm_2d_const, only : STOP_COND_TIME
    use slbm_2d_const, only : STOP_COND_STEADY
    use slbm_2d_const, only : STOP_COND_UNKNOWN
    use slbm_2d_const, only : CS
    use slbm_2d_const, only : CS2
    use slbm_2d_bndry, only : bndry_t
    use tomlf, only : toml_table
    use tomlf, only : toml_error
    use tomlf, only : toml_parse
    use tomlf, only : get_value
    implicit none
    private


    type, public :: config_t
        ! Input configuration parameters
        integer(ip) :: num_x      = 0_ip   ! number of x grid points
        integer(ip) :: num_y      = 0_ip   ! number of y grid points
        real(wp)    :: ds         = 0.0_wp ! mesh step size  ds=dx=dy
        real(wp)    :: kvisc      = 0.0_wp ! kinematic viscosity 
        real(wp)    :: density    = 0.0_wp ! reference density
        real(wp)    :: stop_time  = 0.0_wp ! simulation stop time
        real(wp)    :: stop_etol  = 0.0_wp ! simulation (steady) stop err tol.
        integer(ip) :: stop_cond  = STOP_COND_TIME 

        ! Boundary condition parameters
        type(bndry_t) :: bndry_left
        type(bndry_t) :: bndry_right
        type(bndry_t) :: bndry_top
        type(bndry_t) :: bndry_bottom

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
        character(*),     intent(in)          :: filename
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
        call read_stop_table(config, table_ptr)
        call read_bndry_table(config, table_ptr)

    end function config_from_toml

    function config_constructor( &
            num_x,       &  ! number of z grid points
            num_y,       &  ! number of y grid points
            ds,          &  ! mesh step size
            kvisc,       &  ! kinematic viscosity
            density,     &  ! reference density
            stop_time    &  ! simulation stop time  
            ) result(config)
        integer(ip), intent(in) :: num_x
        integer(ip), intent(in) :: num_y
        real(wp), intent(in)    :: ds
        real(wp), intent(in)    :: kvisc
        real(wp), intent(in)    :: density
        real(wp), intent(in)    :: stop_time
        type(config_t)          :: config

        ! Set input parameters
        config % num_x = num_x
        config % num_y = num_y
        config % ds = ds
        config % kvisc = kvisc
        config % stop_time = stop_time

        ! Set derivied values
        config % dt = ds
        config % nstep = nint(stop_time / (config % dt))
        config % len_x = (num_x - 1.0_wp) * ds 
        config % len_y = (num_y - 1.0_wp) * ds
        config % tau = 0.5_wp + kvisc / (CS2 * config % dt)
    end function config_constructor


    ! config_t type bound procedures
    ! -------------------------------------------------------------------------

    subroutine config_pprint(this)
        class(config_t), intent(in) :: this
        print *, ''
        print *, 'config input parameters'
        print *, '---------------------------------------------'
        print *, 'num_x     = ', this % num_x
        print *, 'num_y     = ', this % num_y
        print *, 'ds        = ', this % ds
        print *, 'kvisc     = ', this % kvisc
        print *, 'density   = ', this % density
        print *, 'stop_time = ', this % stop_time
        print *, 'stop_cond = ', stop_cond_to_string(this % stop_cond)
        print *, 'stop_etol = ', this % stop_etol

        print *, ''
        print *, 'config derived parameters'
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

        call get_value(stop_table, "cond", cond, cond, stat) 
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
        type(config_t), intent(inout)          :: config
        type(toml_table), pointer, intent(in)  :: table
        type(toml_table), pointer              :: bndry_table
        integer(ip)                            :: stat
        real(wp)                               :: nan
        nan = ieee_value(nan, ieee_quiet_nan)

        call get_value(table, "boundary", bndry_table, .false.)
        if (.not. associated(bndry_table) ) then
            print *, "config .toml missing boundary section"
            stop
        endif

        block
            type(toml_table), pointer  :: left_table
            call get_value(bndry_table, "left", left_table, .false.)
            if (.not. associated(bndry_table) ) then
                print *, "config .toml missing boundary.left section"
                stop
            endif
        end block

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


end module slbm_2d_config
