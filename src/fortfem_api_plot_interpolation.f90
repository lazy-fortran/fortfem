module fortfem_api_plot_interpolation
    use fortfem_kinds, only: dp
    use fortfem_api_types, only: function_t, vector_function_t
    use fortfem_nedelec_field_2d, only: evaluate_nedelec_field_2d
    use fortfem_rt_field_2d, only: evaluate_rt_field_2d
    implicit none

    private

    public :: compute_scalar_plot_grid
    public :: compute_vector_plot_grid

contains

    subroutine compute_scalar_plot_grid(uh, nx, ny, x_grid, y_grid, z_grid)
        type(function_t), intent(in) :: uh
        integer, intent(in) :: nx, ny
        real(dp), intent(out) :: x_grid(:)
        real(dp), intent(out) :: y_grid(:)
        real(dp), intent(out) :: z_grid(:, :)

        real(dp) :: x_min, x_max, y_min, y_max, dx_grid, dy_grid
        integer :: i, j

        x_min = minval(uh%space%mesh%data%vertices(1, :))
        x_max = maxval(uh%space%mesh%data%vertices(1, :))
        y_min = minval(uh%space%mesh%data%vertices(2, :))
        y_max = maxval(uh%space%mesh%data%vertices(2, :))

        dx_grid = (x_max - x_min)/real(nx, dp)
        dy_grid = (y_max - y_min)/real(ny, dp)

        do i = 1, nx + 1
            x_grid(i) = x_min + real(i - 1, dp)*dx_grid
        end do

        do j = 1, ny + 1
            y_grid(j) = y_min + real(j - 1, dp)*dy_grid
        end do

        if (uh%space%mesh%data%n_triangles > 0) then
            call interpolate_to_grid(uh, x_grid(1:nx), y_grid(1:ny), z_grid)
        else if (uh%space%mesh%data%n_quads > 0) then
            call interpolate_quad_to_grid(uh, x_grid(1:nx), y_grid(1:ny), &
                z_grid)
        else
            z_grid = 0.0_dp
        end if
    end subroutine compute_scalar_plot_grid

    subroutine compute_vector_plot_grid(Eh, nx, ny, x_grid, y_grid, u_grid, &
            v_grid)
        type(vector_function_t), intent(in) :: Eh
        integer, intent(in) :: nx, ny
        real(dp), intent(out) :: x_grid(:)
        real(dp), intent(out) :: y_grid(:)
        real(dp), intent(out) :: u_grid(:, :)
        real(dp), intent(out) :: v_grid(:, :)

        real(dp) :: x_min, x_max, y_min, y_max, dx_grid, dy_grid
        integer :: i, j

        x_min = minval(Eh%space%mesh%data%vertices(1, :))
        x_max = maxval(Eh%space%mesh%data%vertices(1, :))
        y_min = minval(Eh%space%mesh%data%vertices(2, :))
        y_max = maxval(Eh%space%mesh%data%vertices(2, :))

        dx_grid = (x_max - x_min)/real(nx - 1, dp)
        dy_grid = (y_max - y_min)/real(ny - 1, dp)

        do i = 1, nx
            x_grid(i) = x_min + real(i - 1, dp)*dx_grid
        end do

        do j = 1, ny
            y_grid(j) = y_min + real(j - 1, dp)*dy_grid
        end do

        call interpolate_vector_to_grid(Eh, x_grid, y_grid, u_grid, v_grid)
    end subroutine compute_vector_plot_grid

    subroutine interpolate_to_grid(uh, x_grid, y_grid, z_grid)
        type(function_t), intent(in) :: uh
        real(dp), intent(in) :: x_grid(:), y_grid(:)
        real(dp), intent(out) :: z_grid(:, :)

        integer :: i, j, e, v1, v2, v3
        real(dp) :: x, y, x1, y1, x2, y2, x3, y3
        real(dp) :: lambda1, lambda2, lambda3, val
        logical :: found

        do i = 1, size(x_grid)
            do j = 1, size(y_grid)
                x = x_grid(i)
                y = y_grid(j)
                found = .false.

                do e = 1, uh%space%mesh%data%n_triangles
                    if (found) exit

                    v1 = uh%space%mesh%data%triangles(1, e)
                    v2 = uh%space%mesh%data%triangles(2, e)
                    v3 = uh%space%mesh%data%triangles(3, e)

                    x1 = uh%space%mesh%data%vertices(1, v1)
                    y1 = uh%space%mesh%data%vertices(2, v1)
                    x2 = uh%space%mesh%data%vertices(1, v2)
                    y2 = uh%space%mesh%data%vertices(2, v2)
                    x3 = uh%space%mesh%data%vertices(1, v3)
                    y3 = uh%space%mesh%data%vertices(2, v3)

                    call barycentric_coordinates(x, y, x1, y1, x2, y2, &
                        x3, y3, lambda1, &
                        lambda2, lambda3)

                    if (lambda1 >= -1.0e-10_dp .and. &
                        lambda2 >= -1.0e-10_dp .and. &
                        lambda3 >= -1.0e-10_dp) then
                        val = lambda1*uh%values(v1) + &
                            lambda2*uh%values(v2) + &
                            lambda3*uh%values(v3)
                        z_grid(i, j) = val
                        found = .true.
                    end if
                end do

                if (.not. found) then
                    z_grid(i, j) = find_nearest_value(uh, x, y)
                end if
            end do
        end do
    end subroutine interpolate_to_grid

    subroutine interpolate_quad_to_grid(uh, x_grid, y_grid, z_grid)
        type(function_t), intent(in) :: uh
        real(dp), intent(in) :: x_grid(:), y_grid(:)
        real(dp), intent(out) :: z_grid(:, :)

        integer :: i, j

        z_grid = 0.0_dp

        do i = 1, size(x_grid)
            do j = 1, size(y_grid)
                z_grid(i, j) = interpolate_quad_value(uh, x_grid(i), &
                    y_grid(j))
            end do
        end do
    end subroutine interpolate_quad_to_grid

    subroutine interpolate_vector_to_grid(Eh, x_grid, y_grid, u_grid, v_grid)
        type(vector_function_t), intent(in) :: Eh
        real(dp), intent(in) :: x_grid(:), y_grid(:)
        real(dp), intent(out) :: u_grid(:, :), v_grid(:, :)

        complex(dp), allocatable :: dofs(:)
        complex(dp) :: complex_value(2)
        real(dp) :: value(2)
        integer :: i, j, triangle
        logical :: is_nedelec, is_rt, is_nodal

        is_nedelec = trim(Eh%space%element_family) == "Nedelec" .or. &
            trim(Eh%space%element_family) == "Nedelec1" .or. &
            trim(Eh%space%element_family) == "Edge"
        is_rt = trim(Eh%space%element_family) == "RT" .or. &
            trim(Eh%space%element_family) == "Raviart-Thomas"
        is_nodal = size(Eh%values, 1) == Eh%space%mesh%data%n_vertices

        if (.not. is_nedelec .and. .not. is_rt .and. .not. is_nodal) then
            error stop "vector plotting cannot reconstruct this FE space"
        end if
        if (is_nedelec .or. is_rt) then
            if (size(Eh%values, 1) /= Eh%space%mesh%data%n_edges) then
                error stop "vector plotting needs lowest-order edge DOFs"
            end if
            allocate(dofs(size(Eh%values, 1)))
            ! The solver stores the real FE coefficient vector in column one;
            ! column two is reserved for the second source component.
            dofs = cmplx(Eh%values(:, 1), 0.0_dp, dp)
        end if

        u_grid = 0.0_dp
        v_grid = 0.0_dp
        do i = 1, size(x_grid)
            do j = 1, size(y_grid)
                call locate_triangle(Eh%space%mesh%data, x_grid(i), &
                    y_grid(j), triangle)
                if (triangle > 0) then
                    if (is_nedelec) then
                        call evaluate_nedelec_field_2d( &
                            Eh%space%mesh%data, triangle, x_grid(i), &
                            y_grid(j), dofs, complex_value)
                        value = real(complex_value, dp)
                    else if (is_rt) then
                        call evaluate_rt_field_2d( &
                            Eh%space%mesh%data, triangle, x_grid(i), &
                            y_grid(j), dofs, complex_value)
                        value = real(complex_value, dp)
                    else
                        call interpolate_nodal_vector( &
                            Eh, triangle, x_grid(i), y_grid(j), value)
                    end if
                    u_grid(i, j) = value(1)
                    v_grid(i, j) = value(2)
                else if (is_nodal) then
                    call nearest_nodal_vector(Eh, x_grid(i), y_grid(j), value)
                    u_grid(i, j) = value(1)
                    v_grid(i, j) = value(2)
                end if
            end do
        end do

        if (allocated(dofs)) deallocate(dofs)
    end subroutine interpolate_vector_to_grid

    subroutine locate_triangle(mesh, x, y, triangle)
        use fortfem_mesh_2d, only: mesh_2d_t
        type(mesh_2d_t), intent(in) :: mesh
        real(dp), intent(in) :: x, y
        integer, intent(out) :: triangle

        real(dp) :: lambda1, lambda2, lambda3
        integer :: candidate, v1, v2, v3

        triangle = 0
        do candidate = 1, mesh%n_triangles
            v1 = mesh%triangles(1, candidate)
            v2 = mesh%triangles(2, candidate)
            v3 = mesh%triangles(3, candidate)
            call barycentric_coordinates( &
                x, y, mesh%vertices(1, v1), mesh%vertices(2, v1), &
                mesh%vertices(1, v2), mesh%vertices(2, v2), &
                mesh%vertices(1, v3), mesh%vertices(2, v3), &
                lambda1, lambda2, lambda3)
            if (lambda1 >= -1.0e-10_dp .and. &
                lambda2 >= -1.0e-10_dp .and. &
                lambda3 >= -1.0e-10_dp) then
                triangle = candidate
                return
            end if
        end do
    end subroutine locate_triangle

    subroutine interpolate_nodal_vector(Eh, triangle, x, y, value)
        type(vector_function_t), intent(in) :: Eh
        integer, intent(in) :: triangle
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        integer :: vertex(3), component
        real(dp) :: lambda(3)

        vertex = Eh%space%mesh%data%triangles(:, triangle)
        call barycentric_coordinates( &
            x, y, Eh%space%mesh%data%vertices(1, vertex(1)), &
            Eh%space%mesh%data%vertices(2, vertex(1)), &
            Eh%space%mesh%data%vertices(1, vertex(2)), &
            Eh%space%mesh%data%vertices(2, vertex(2)), &
            Eh%space%mesh%data%vertices(1, vertex(3)), &
            Eh%space%mesh%data%vertices(2, vertex(3)), &
            lambda(1), lambda(2), lambda(3))
        do component = 1, 2
            value(component) = dot_product(lambda, &
                Eh%values(vertex, component))
        end do
    end subroutine interpolate_nodal_vector

    subroutine nearest_nodal_vector(Eh, x, y, value)
        type(vector_function_t), intent(in) :: Eh
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        integer :: component, vertex, nearest
        real(dp) :: distance, minimum, dx, dy

        nearest = 1
        minimum = huge(1.0_dp)
        do vertex = 1, Eh%space%mesh%data%n_vertices
            dx = x - Eh%space%mesh%data%vertices(1, vertex)
            dy = y - Eh%space%mesh%data%vertices(2, vertex)
            distance = dx*dx + dy*dy
            if (distance < minimum) then
                minimum = distance
                nearest = vertex
            end if
        end do
        do component = 1, 2
            value(component) = Eh%values(nearest, component)
        end do
    end subroutine nearest_nodal_vector

    subroutine barycentric_coordinates(x, y, x1, y1, x2, y2, x3, y3, &
            lambda1, lambda2, lambda3)
        real(dp), intent(in) :: x, y, x1, y1, x2, y2, x3, y3
        real(dp), intent(out) :: lambda1, lambda2, lambda3

        real(dp) :: denom

        denom = (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3)

        if (abs(denom) < 1.0e-14_dp) then
            lambda1 = -1.0_dp
            lambda2 = -1.0_dp
            lambda3 = -1.0_dp
        else
            lambda1 = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3))/denom
            lambda2 = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3))/denom
            lambda3 = 1.0_dp - lambda1 - lambda2
        end if
    end subroutine barycentric_coordinates

    pure function find_nearest_value(uh, x, y) result(val)
        type(function_t), intent(in) :: uh
        real(dp), intent(in) :: x, y
        real(dp) :: val

        integer :: i, nearest_vertex
        real(dp) :: min_dist, dist, vx, vy

        min_dist = huge(1.0_dp)
        nearest_vertex = 1

        do i = 1, uh%space%mesh%data%n_vertices
            vx = uh%space%mesh%data%vertices(1, i)
            vy = uh%space%mesh%data%vertices(2, i)
            dist = (x - vx)**2 + (y - vy)**2

            if (dist < min_dist) then
                min_dist = dist
                nearest_vertex = i
            end if
        end do

        val = uh%values(nearest_vertex)
    end function find_nearest_value

    pure subroutine quad_shape_functions(xi_ref, eta_ref, N)
        real(dp), intent(in) :: xi_ref, eta_ref
        real(dp), intent(out) :: N(4)

        N(1) = 0.25_dp*(1.0_dp - xi_ref)*(1.0_dp - eta_ref)
        N(2) = 0.25_dp*(1.0_dp + xi_ref)*(1.0_dp - eta_ref)
        N(3) = 0.25_dp*(1.0_dp + xi_ref)*(1.0_dp + eta_ref)
        N(4) = 0.25_dp*(1.0_dp - xi_ref)*(1.0_dp + eta_ref)
    end subroutine quad_shape_functions

    pure function interpolate_quad_value(uh, x, y) result(val)
        type(function_t), intent(in) :: uh
        real(dp), intent(in) :: x, y
        real(dp) :: val

        integer :: q, k, vi
        integer :: v_ids(4)
        real(dp) :: x1, y1, x2, y2
        real(dp) :: xc, yc, xi_ref, eta_ref
        real(dp) :: N(4)
        logical :: found

        val = 0.0_dp
        found = .false.

        do q = 1, uh%space%mesh%data%n_quads
            v_ids = uh%space%mesh%data%quads(:, q)

            x1 = uh%space%mesh%data%vertices(1, v_ids(1))
            y1 = uh%space%mesh%data%vertices(2, v_ids(1))
            x2 = uh%space%mesh%data%vertices(1, v_ids(3))
            y2 = uh%space%mesh%data%vertices(2, v_ids(3))

            if (x >= x1 .and. x <= x2 .and. y >= y1 .and. y <= y2) then
                if (x2 > x1 .and. y2 > y1) then
                    xc = 0.5_dp*(x1 + x2)
                    yc = 0.5_dp*(y1 + y2)

                    xi_ref = 2.0_dp*(x - xc)/(x2 - x1)
                    eta_ref = 2.0_dp*(y - yc)/(y2 - y1)

                    call quad_shape_functions(xi_ref, eta_ref, N)

                    val = 0.0_dp
                    do k = 1, 4
                        vi = v_ids(k)
                        val = val + N(k)*uh%values(vi)
                    end do

                    found = .true.
                end if
            end if

            if (found) exit
        end do

        if (.not. found) then
            val = find_nearest_value(uh, x, y)
        end if
    end function interpolate_quad_value

end module fortfem_api_plot_interpolation
