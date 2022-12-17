module slbm_2d_motion

    use slbm_2d_kinds, only : wp, ip
    implicit none
    private

    ! Abstract class for body motions
    type, public, abstract :: motion_t
    contains
        private
        procedure(iface_update), deferred :: update
    end type motion_t


    abstract interface
        subroutine iface_update(this, time)
            import :: motion_t, wp
            class(motion_t), intent(inout) :: this
            real(wp),        intent(in)    :: time
        end subroutine iface_update
    end interface


    ! Prescribed rigid body motions
    type, public, extends(motion_t) :: motion_pre_rigid_t
    contains
        private
        procedure, public :: update => update_pre_rigid
    end type motion_pre_rigid_t


    ! Prescribed motions from points file
    type, public, extends(motion_t) :: motion_pre_point_t
    contains
        private
        procedure, public :: update => update_pre_point
    end type motion_pre_point_t


    !! Dyanmic rigid body motion 
    !type, public, extends(motion_t) :: motion_dyn_rigid_t
    !end type motion_dyn_rigid_t


contains

    subroutine update_pre_rigid(this, time)
        class(motion_pre_rigid_t), intent(inout) :: this
        real(wp), intent(in)                     :: time
    end subroutine update_pre_rigid


    subroutine update_pre_point(this, time)
        class(motion_pre_point_t), intent(inout) :: this
        real(wp), intent(in)                     :: time
    end subroutine update_pre_point


end module slbm_2d_motion
