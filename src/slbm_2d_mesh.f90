module slbm_2d_mesh
    use slbm_2d_kinds,  only : wp, ip
    use slbm_2d_config, only : config_t
    implicit none
    private


    type, public :: mesh_t
        real(wp), allocatable    :: x(:,:)
        real(wp), allocatable    :: y(:,:)
        logical,  allocatable    :: is_object(:,:)
    contains
        private
        procedure, public :: num_x => mesh_num_x
        procedure, public :: num_y => mesh_num_y
        procedure, public :: dx    => mesh_dx
        procedure, public :: dy    => mesh_dy
        procedure, public :: xmin  => mesh_xmin
        procedure, public :: xmax  => mesh_xmax
        procedure, public :: ymin  => mesh_ymin
        procedure, public :: ymax  => mesh_ymax
        procedure, public :: bounding_box => mesh_bounding_box
    end type mesh_t

    interface mesh_t
        procedure :: mesh_constructor
    end interface mesh_t

contains

    function mesh_constructor(config) result(mesh)
        type(config_t), intent(in) :: config
        type(mesh_t)               :: mesh

        real(wp), allocatable      :: x(:)
        real(wp), allocatable      :: y(:)
        integer(ip)                :: i
        integer(ip)                :: j

        ! Create x mesh
        allocate(mesh % x(config % num_x, config % num_y))
        x = [(config % ds*(i - 1.0_wp), i = 1, config % num_x)]
        do j = 1, config % num_y
            mesh % x(:,j) = x
        end do 

        ! Create y mesh
        allocate(mesh % y(config % num_x, config % num_y))
        y = [(config % ds*(j - 1.0_wp), j = 1, config % num_y)]
        do i = 1, config % num_x
            mesh % y(i,:) = y
        end do

        ! Create id markers
        allocate(mesh % is_object(config % num_x, config % num_y))
        mesh % is_object = .false.

    end function mesh_constructor


    function mesh_num_x(this) result(num_x)
        class(mesh_t), intent(in) :: this
        integer(ip)               :: num_x
        if (allocated(this % x)) then
            num_x = size(this % x, 1)
        else
            num_x = 0
        end if
    end function mesh_num_x


    function mesh_num_y(this) result(num_y)
        class(mesh_t), intent(in) :: this
        integer(ip)               :: num_y
        if (allocated(this % y)) then
            num_y = size(this % y, 2)
        else
            num_y = 0
        end if
    end function mesh_num_y


    function mesh_dx(this) result(dx)
        class(mesh_t), intent(in) :: this
        real(wp)                  :: dx
        if (allocated(this % x)) then
            dx = this % x(2,1) - this % x(1,1)
        else
            dx = 0.0_wp
        end if
    end function mesh_dx


    function mesh_dy(this) result(dy)
        class(mesh_t), intent(in) :: this
        real(wp)                  :: dy
        if (allocated(this % y)) then
            dy = this % y(1,2) - this % y(1,1)
        else
            dy = 0.0_wp
        end if
    end function mesh_dy


    function mesh_xmin(this) result(xmin)
        class(mesh_t), intent(in) :: this
        real(wp)                  :: xmin
        xmin = 0.0_wp
    end function mesh_xmin


    function mesh_xmax(this) result(xmax)
        class(mesh_t), intent(in) :: this
        real(wp)                  :: xmax
        xmax = this % dx() * (this % num_x() - 1.0_wp)
    end function mesh_xmax


    function mesh_ymin(this) result(ymin)
        class(mesh_t), intent(in) :: this
        real(wp)                  :: ymin
        ymin = 0.0_wp
    end function mesh_ymin


    function mesh_ymax(this) result(ymax)
        class(mesh_t), intent(in) :: this
        real(wp)                  :: ymax
        ymax = this % dy() * (this % num_y() - 1.0_wp)
    end function mesh_ymax


    subroutine mesh_bounding_box(this, xmin, xmax, ymin, ymax)
        class(mesh_t), intent(in) :: this
        real(wp), intent(out)     :: xmin
        real(wp), intent(out)     :: xmax
        real(wp), intent(out)     :: ymin
        real(wp), intent(out)     :: ymax
        xmin = this % xmin()
        xmax = this % xmax()
        ymin = this % ymin()
        ymax = this % ymax()
    end subroutine mesh_bounding_box

end module slbm_2d_mesh
