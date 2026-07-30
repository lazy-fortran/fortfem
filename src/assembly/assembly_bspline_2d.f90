module fortfem_assembly_bspline_2d
    !! Sparse scalar isogeometric assembly on one rational tensor patch.
    use fortfem_bspline_feec, only: &
        evaluate_bspline_basis, evaluate_nurbs_surface_geometry
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_bspline_h1_operator_csc

contains

    subroutine assemble_bspline_h1_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient, &
            mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: quadrature_weights_x(:)
        real(dp), allocatable :: quadrature_weights_y(:), triplet_values(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: determinant, geometry_jacobian(2, 2), geometry_point(2)
        real(dp) :: basis_column, basis_row, gradient_column(2), gradient_row(2)
        real(dp) :: inverse(2, 2)
        real(dp) :: mass_weight, physical_weight, stiffness_weight
        integer :: entry, local_column, local_count, local_row
        integer :: local_status, max_entries, nx, ny, point_x, point_y
        integer :: span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric H1 assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*(degree_y + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), &
            local_dofs(local_count), local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            quadrature_weights_x(quadrature_order), &
            quadrature_weights_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, quadrature_weights_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, quadrature_weights_x)
                call build_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), value_y, &
                        derivative_y, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), value_x, &
                            derivative_x, local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant* &
                            quadrature_weights_x(point_x)* &
                            quadrature_weights_y(point_y)
                        do local_column = 1, local_count
                            call local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, value_x, derivative_x, value_y, &
                                derivative_y, basis_column, gradient_column)
                            gradient_column = &
                                matmul(transpose(inverse), gradient_column)
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                gradient_row = &
                                    matmul(transpose(inverse), gradient_row)
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(stiffness_weight* &
                                    dot_product(gradient_row, gradient_column) + &
                                    mass_weight*basis_row*basis_column)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            nx*ny, nx*ny, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_h1_operator_csc

    subroutine local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, value_x, &
            derivative_x, value_y, derivative_y, basis_value, gradient)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: value_x(:), derivative_x(:)
        real(dp), intent(in) :: value_y(:), derivative_y(:)
        real(dp), intent(out) :: basis_value
        real(dp), intent(out) :: gradient(2)

        integer :: basis_x, basis_y, offset_x, offset_y

        offset_x = modulo(local_dof - 1, degree_x + 1)
        offset_y = (local_dof - 1)/(degree_x + 1)
        basis_x = span_x - degree_x + offset_x
        basis_y = span_y - degree_y + offset_y
        basis_value = value_x(basis_x)*value_y(basis_y)
        gradient = [ &
            derivative_x(basis_x)*value_y(basis_y), &
            value_x(basis_x)*derivative_y(basis_y)]
    end subroutine local_basis_data

    pure subroutine build_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof

        local_dof = 0
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*nx
            end do
        end do
    end subroutine build_local_dofs

    pure integer function positive_span_count( &
            knots, degree, basis_count) result(count)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree, basis_count
        integer :: span

        count = 0
        do span = degree + 1, basis_count
            if (knots(span + 1) > knots(span)) count = count + 1
        end do
    end function positive_span_count

    pure subroutine inverse_2d(matrix, inverse, determinant, status)
        real(dp), intent(in) :: matrix(2, 2)
        real(dp), intent(out) :: inverse(2, 2), determinant
        integer, intent(out) :: status

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        status = 1
        inverse = 0.0_dp
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(matrix))**2)) return
        inverse = reshape([ &
            matrix(2, 2), -matrix(2, 1), &
            -matrix(1, 2), matrix(1, 1)], [2, 2])/determinant
        status = 0
    end subroutine inverse_2d

end module fortfem_assembly_bspline_2d
