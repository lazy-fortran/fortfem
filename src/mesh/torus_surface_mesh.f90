module fortfem_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: generate_torus_surface_mesh

contains

    subroutine generate_torus_surface_mesh( &
            major_radius, minor_radius, poloidal_count, toroidal_count, &
            vertices, triangles, parameter_coordinates)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: poloidal_count, toroidal_count
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)
        real(dp), allocatable, intent(out), optional :: parameter_coordinates(:, :)

        real(dp) :: phi, theta
        integer :: cell, phi_index, theta_index
        integer :: lower_left, lower_right, upper_left, upper_right
        integer :: vertex

        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) then
            error stop "Torus mesh requires major_radius > minor_radius > 0"
        end if
        if (poloidal_count < 3 .or. toroidal_count < 3) then
            error stop "Torus mesh requires at least three nodes per direction"
        end if

        allocate(vertices(3, poloidal_count*toroidal_count))
        allocate(triangles(3, 2*poloidal_count*toroidal_count))
        if (present(parameter_coordinates)) then
            allocate(parameter_coordinates(2, poloidal_count*toroidal_count))
        end if
        do phi_index = 0, toroidal_count - 1
            phi = 2.0_dp*acos(-1.0_dp)*real(phi_index, dp)/ &
                real(toroidal_count, dp)
            do theta_index = 0, poloidal_count - 1
                theta = 2.0_dp*acos(-1.0_dp)*real(theta_index, dp)/ &
                    real(poloidal_count, dp)
                vertex = node_index(theta_index, phi_index)
                vertices(:, vertex) = [ &
                    (major_radius + minor_radius*cos(theta))*cos(phi), &
                    (major_radius + minor_radius*cos(theta))*sin(phi), &
                    minor_radius*sin(theta)]
                if (present(parameter_coordinates)) then
                    parameter_coordinates(:, vertex) = [theta, phi]
                end if
            end do
        end do

        cell = 0
        do phi_index = 0, toroidal_count - 1
            do theta_index = 0, poloidal_count - 1
                lower_left = node_index(theta_index, phi_index)
                lower_right = node_index( &
                    theta_index, modulo(phi_index + 1, toroidal_count))
                upper_left = node_index( &
                    modulo(theta_index + 1, poloidal_count), phi_index)
                upper_right = node_index( &
                    modulo(theta_index + 1, poloidal_count), &
                    modulo(phi_index + 1, toroidal_count))
                cell = cell + 1
                triangles(:, 2*cell - 1) = &
                    [lower_left, lower_right, upper_left]
                triangles(:, 2*cell) = &
                    [upper_left, lower_right, upper_right]
            end do
        end do

    contains

        pure integer function node_index(theta_node, phi_node) result(index)
            integer, intent(in) :: theta_node, phi_node

            index = phi_node*poloidal_count + theta_node + 1
        end function node_index
    end subroutine generate_torus_surface_mesh

end module fortfem_torus_surface_mesh
