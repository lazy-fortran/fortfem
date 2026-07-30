module fortfem_solid_torus_tetra_mesh
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: generate_solid_torus_tetra_mesh

contains

    subroutine generate_solid_torus_tetra_mesh( &
            major_radius, minor_radius, radial_layers, poloidal_count, &
            toroidal_count, vertices, tetrahedra, boundary_triangles, &
            parameter_coordinates)
        real(dp), intent(in) :: major_radius, minor_radius
        integer, intent(in) :: radial_layers, poloidal_count, toroidal_count
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: tetrahedra(:, :)
        integer, allocatable, intent(out) :: boundary_triangles(:, :)
        real(dp), allocatable, intent(out) :: parameter_coordinates(:, :)

        integer, allocatable :: disk_triangles(:, :)
        integer :: boundary_cell, disk_triangle, layer, local_vertex
        integer :: phi_index, tetrahedron, theta_index
        integer :: lower(3), sorted(3), upper(3), vertex
        integer :: vertices_per_slice
        real(dp) :: phi, radius, theta

        if (major_radius <= minor_radius .or. minor_radius <= 0.0_dp) then
            error stop "Solid torus requires major_radius > minor_radius > 0"
        end if
        if (radial_layers < 1) then
            error stop "Solid torus requires at least one radial layer"
        end if
        if (poloidal_count < 3 .or. toroidal_count < 3) then
            error stop "Solid torus requires at least three angular nodes"
        end if
        vertices_per_slice = 1 + radial_layers*poloidal_count
        allocate(vertices(3, vertices_per_slice*toroidal_count))
        allocate(parameter_coordinates(2, size(vertices, 2)))
        do phi_index = 0, toroidal_count - 1
            phi = 2.0_dp*acos(-1.0_dp)*real(phi_index, dp)/ &
                real(toroidal_count, dp)
            vertex = volume_vertex(phi_index, 1)
            vertices(:, vertex) = [ &
                major_radius*cos(phi), major_radius*sin(phi), 0.0_dp]
            parameter_coordinates(:, vertex) = [0.0_dp, phi]
            do layer = 1, radial_layers
                radius = minor_radius*real(layer, dp)/real(radial_layers, dp)
                do theta_index = 0, poloidal_count - 1
                    theta = 2.0_dp*acos(-1.0_dp)*real(theta_index, dp)/ &
                        real(poloidal_count, dp)
                    local_vertex = ring_vertex(layer, theta_index)
                    vertex = volume_vertex(phi_index, local_vertex)
                    vertices(:, vertex) = [ &
                        (major_radius + radius*cos(theta))*cos(phi), &
                        (major_radius + radius*cos(theta))*sin(phi), &
                        radius*sin(theta)]
                    parameter_coordinates(:, vertex) = [theta, phi]
                end do
            end do
        end do

        call build_disk_triangles( &
            radial_layers, poloidal_count, disk_triangles)
        allocate(tetrahedra( &
            4, 3*size(disk_triangles, 2)*toroidal_count))
        tetrahedron = 0
        do phi_index = 0, toroidal_count - 1
            do disk_triangle = 1, size(disk_triangles, 2)
                sorted = disk_triangles(:, disk_triangle)
                call sort_three(sorted)
                lower = [ &
                    volume_vertex(phi_index, sorted(1)), &
                    volume_vertex(phi_index, sorted(2)), &
                    volume_vertex(phi_index, sorted(3))]
                upper = [ &
                    volume_vertex(modulo(phi_index + 1, toroidal_count), &
                    sorted(1)), &
                    volume_vertex(modulo(phi_index + 1, toroidal_count), &
                    sorted(2)), &
                    volume_vertex(modulo(phi_index + 1, toroidal_count), &
                    sorted(3))]
                tetrahedron = tetrahedron + 1
                tetrahedra(:, tetrahedron) = &
                    [lower(1), lower(2), lower(3), upper(3)]
                tetrahedron = tetrahedron + 1
                tetrahedra(:, tetrahedron) = &
                    [lower(1), lower(2), upper(2), upper(3)]
                tetrahedron = tetrahedron + 1
                tetrahedra(:, tetrahedron) = &
                    [lower(1), upper(1), upper(2), upper(3)]
                call orient_positive( &
                    vertices, tetrahedra(:, tetrahedron - 2))
                call orient_positive( &
                    vertices, tetrahedra(:, tetrahedron - 1))
                call orient_positive(vertices, tetrahedra(:, tetrahedron))
            end do
        end do

        allocate(boundary_triangles( &
            3, 2*poloidal_count*toroidal_count))
        boundary_cell = 0
        do phi_index = 0, toroidal_count - 1
            do theta_index = 0, poloidal_count - 1
                boundary_cell = boundary_cell + 1
                boundary_triangles(:, 2*boundary_cell - 1) = [ &
                    volume_vertex(phi_index, &
                    ring_vertex(radial_layers, theta_index)), &
                    volume_vertex(modulo(phi_index + 1, toroidal_count), &
                    ring_vertex(radial_layers, theta_index)), &
                    volume_vertex(phi_index, &
                    ring_vertex(radial_layers, &
                    modulo(theta_index + 1, poloidal_count)))]
                boundary_triangles(:, 2*boundary_cell) = [ &
                    volume_vertex(phi_index, &
                    ring_vertex(radial_layers, &
                    modulo(theta_index + 1, poloidal_count))), &
                    volume_vertex(modulo(phi_index + 1, toroidal_count), &
                    ring_vertex(radial_layers, theta_index)), &
                    volume_vertex(modulo(phi_index + 1, toroidal_count), &
                    ring_vertex(radial_layers, &
                    modulo(theta_index + 1, poloidal_count)))]
            end do
        end do

    contains

        pure integer function volume_vertex(phi_node, local_node) result(index)
            integer, intent(in) :: phi_node, local_node

            index = phi_node*vertices_per_slice + local_node
        end function volume_vertex

        pure integer function ring_vertex( &
                radial_layer, theta_node) result(index)
            integer, intent(in) :: radial_layer, theta_node

            index = 2 + (radial_layer - 1)*poloidal_count + theta_node
        end function ring_vertex

    end subroutine generate_solid_torus_tetra_mesh

    subroutine build_disk_triangles( &
            radial_layers, poloidal_count, triangles)
        integer, intent(in) :: radial_layers, poloidal_count
        integer, allocatable, intent(out) :: triangles(:, :)

        integer :: cell, layer, theta_index
        integer :: inner, inner_next, outer, outer_next

        allocate(triangles(3, poloidal_count*(2*radial_layers - 1)))
        cell = 0
        do theta_index = 0, poloidal_count - 1
            cell = cell + 1
            triangles(:, cell) = [ &
                1, disk_ring_vertex(1, theta_index, poloidal_count), &
                disk_ring_vertex( &
                1, modulo(theta_index + 1, poloidal_count), poloidal_count)]
        end do
        do layer = 2, radial_layers
            do theta_index = 0, poloidal_count - 1
                inner = &
                    disk_ring_vertex(layer - 1, theta_index, poloidal_count)
                inner_next = disk_ring_vertex( &
                    layer - 1, modulo(theta_index + 1, poloidal_count), &
                    poloidal_count)
                outer = disk_ring_vertex(layer, theta_index, poloidal_count)
                outer_next = disk_ring_vertex( &
                    layer, modulo(theta_index + 1, poloidal_count), &
                    poloidal_count)
                cell = cell + 1
                triangles(:, cell) = [inner, outer, inner_next]
                cell = cell + 1
                triangles(:, cell) = [inner_next, outer, outer_next]
            end do
        end do
    end subroutine build_disk_triangles

    pure integer function disk_ring_vertex( &
            layer, theta_index, poloidal_count) result(index)
        integer, intent(in) :: layer, theta_index, poloidal_count

        index = 2 + (layer - 1)*poloidal_count + theta_index
    end function disk_ring_vertex

    pure subroutine sort_three(values)
        integer, intent(inout) :: values(3)

        integer :: first, second, temporary

        do first = 1, 2
            do second = first + 1, 3
                if (values(second) < values(first)) then
                    temporary = values(first)
                    values(first) = values(second)
                    values(second) = temporary
                end if
            end do
        end do
    end subroutine sort_three

    pure subroutine orient_positive(vertices, tetrahedron)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(inout) :: tetrahedron(4)

        integer :: temporary
        real(dp) :: determinant

        determinant = dot_product( &
            vertices(:, tetrahedron(2)) - vertices(:, tetrahedron(1)), &
            cross_product( &
            vertices(:, tetrahedron(3)) - vertices(:, tetrahedron(1)), &
            vertices(:, tetrahedron(4)) - vertices(:, tetrahedron(1))))
        if (determinant < 0.0_dp) then
            temporary = tetrahedron(3)
            tetrahedron(3) = tetrahedron(4)
            tetrahedron(4) = temporary
        end if
    end subroutine orient_positive

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_solid_torus_tetra_mesh
