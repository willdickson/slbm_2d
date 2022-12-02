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

end module slbm_2d_mesh
