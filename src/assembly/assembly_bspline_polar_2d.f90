module fortfem_assembly_bspline_polar_2d
    !! Sparse Galerkin restriction from periodic tensor to polar bases.
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_bspline_polar, only: &
        build_bspline_polar_feec_2d_extractions, &
        build_bspline_polar_h1_extraction, &
        evaluate_periodic_bspline_basis
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: &
        csc_from_triplet, csc_matmul, csc_t, csc_transpose, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    public :: restrict_bspline_polar_operator_csc
    public :: assemble_bspline_polar_h1_operator_csc
    public :: assemble_bspline_polar_hcurl_operator_csc
    public :: assemble_bspline_polar_l2_mass_csc

contains

    subroutine assemble_bspline_polar_hcurl_operator_csc( &
            radial_knots, radial_degree, azimuth_count, azimuth_degree, &
            control_points, weights, quadrature_order, matrix, status, &
            curl_coefficient, mass_coefficient)
        !! Covariant-Piola H(curl) Hodge and curl-curl operator.
        real(dp), intent(in) :: radial_knots(:)
        integer, intent(in) :: radial_degree, azimuth_count, azimuth_degree
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: quadrature_order
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        real(dp), allocatable :: h1_extraction(:, :)
        real(dp), allocatable :: hcurl_extraction(:, :), l2_extraction(:, :)
        type(csc_t) :: tensor_matrix
        real(dp) :: curl_weight, mass_weight
        integer :: local_status, radial_count

        curl_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        call assemble_tensor_polar_form( &
            1, radial_knots, radial_degree, azimuth_count, azimuth_degree, &
            control_points, weights, quadrature_order, tensor_matrix, status, &
            curl_weight, mass_weight)
        if (status%code /= 0) return
        radial_count = size(radial_knots) - radial_degree - 1
        call build_bspline_polar_feec_2d_extractions( &
            azimuth_count, radial_count, h1_extraction, hcurl_extraction, &
            l2_extraction, local_status)
        if (local_status /= 0) return
        call restrict_bspline_polar_operator_csc( &
            hcurl_extraction, tensor_matrix, matrix, status)
    end subroutine assemble_bspline_polar_hcurl_operator_csc

    subroutine assemble_bspline_polar_l2_mass_csc( &
            radial_knots, radial_degree, azimuth_count, azimuth_degree, &
            control_points, weights, quadrature_order, matrix, status)
        !! Contravariant-density L2 Hodge operator.
        real(dp), intent(in) :: radial_knots(:)
        integer, intent(in) :: radial_degree, azimuth_count, azimuth_degree
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: quadrature_order
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: h1_extraction(:, :)
        real(dp), allocatable :: hcurl_extraction(:, :), l2_extraction(:, :)
        type(csc_t) :: tensor_matrix
        integer :: local_status, radial_count

        call assemble_tensor_polar_form( &
            2, radial_knots, radial_degree, azimuth_count, azimuth_degree, &
            control_points, weights, quadrature_order, tensor_matrix, status, &
            0.0_dp, 1.0_dp)
        if (status%code /= 0) return
        radial_count = size(radial_knots) - radial_degree - 1
        call build_bspline_polar_feec_2d_extractions( &
            azimuth_count, radial_count, h1_extraction, hcurl_extraction, &
            l2_extraction, local_status)
        if (local_status /= 0) return
        call restrict_bspline_polar_operator_csc( &
            l2_extraction, tensor_matrix, matrix, status)
    end subroutine assemble_bspline_polar_l2_mass_csc

    subroutine assemble_tensor_polar_form( &
            form_degree, radial_knots, radial_degree, azimuth_count, &
            azimuth_degree, control_points, weights, quadrature_order, &
            matrix, status, differential_weight, mass_weight)
        integer, intent(in) :: form_degree
        real(dp), intent(in) :: radial_knots(:)
        integer, intent(in) :: radial_degree, azimuth_count, azimuth_degree
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: quadrature_order
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in) :: differential_weight, mass_weight

        real(dp), allocatable :: angular_d(:), angular_n(:), angular_n_d(:)
        real(dp), allocatable :: dense(:, :), radial_d(:), radial_n(:)
        real(dp), allocatable :: ignored_derivative(:)
        real(dp), allocatable :: radial_n_d(:), reduced_knots(:)
        real(dp), allocatable :: reference_curl(:), reference_value(:, :)
        real(dp), allocatable :: scalar_value(:)
        real(dp), allocatable :: radial_nodes(:), radial_weights(:)
        real(dp) :: angular_nodes(quadrature_order)
        real(dp) :: angular_weights(quadrature_order)
        real(dp) :: determinant, inverse(2, 2), jacobian(2, 2)
        real(dp) :: physical_column(2), physical_row(2), quadrature_weight
        integer :: angular, angular_cell, angular_point, column, dof_count
        integer :: local_status, radial, radial_count, radial_point, span
        integer :: row, tensor_radial_count

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Physical polar form assembly requires a regular mapped patch")
        radial_count = size(radial_knots) - radial_degree - 1
        if (form_degree < 1 .or. form_degree > 2) return
        if (radial_degree < 1 .or. azimuth_degree < 1) return
        if (quadrature_order < 1) return
        if (radial_count < 5 .or. azimuth_count < azimuth_degree + 1) return
        if (size(control_points, 1) /= 2) return
        if (any(shape(control_points(1, :, :)) /= &
            [azimuth_count, radial_count])) return
        if (any(shape(weights) /= [azimuth_count, radial_count])) return
        if (any(weights <= 0.0_dp)) return

        tensor_radial_count = (radial_count - 1)*azimuth_count
        if (form_degree == 1) then
            dof_count = tensor_radial_count + radial_count*azimuth_count
            allocate(reference_value(2, dof_count), reference_curl(dof_count))
        else
            dof_count = tensor_radial_count
            allocate(scalar_value(dof_count))
        end if
        allocate(dense(dof_count, dof_count))
        allocate(reduced_knots(size(radial_knots) - 2))
        allocate(radial_nodes(quadrature_order), radial_weights(quadrature_order))
        reduced_knots = radial_knots(2:size(radial_knots) - 1)
        dense = 0.0_dp
        do span = radial_degree + 1, radial_count
            if (radial_knots(span + 1) <= radial_knots(span)) cycle
            call gauss_legendre_ab( &
                quadrature_order, radial_knots(span), &
                radial_knots(span + 1), radial_nodes, radial_weights)
            do angular_cell = 1, azimuth_count
                call gauss_legendre_ab( &
                    quadrature_order, &
                    real(angular_cell - 1, dp)/real(azimuth_count, dp), &
                    real(angular_cell, dp)/real(azimuth_count, dp), &
                    angular_nodes, angular_weights)
                do angular_point = 1, quadrature_order
                    call evaluate_periodic_bspline_basis( &
                        azimuth_count, azimuth_degree, &
                        angular_nodes(angular_point), angular_n, angular_n_d, &
                        local_status)
                    if (local_status /= 0) return
                    call evaluate_periodic_bspline_basis( &
                        azimuth_count, azimuth_degree - 1, &
                        angular_nodes(angular_point), angular_d, &
                        ignored_derivative, local_status)
                    if (local_status /= 0) return
                    angular_d = real(azimuth_count, dp)*angular_d
                    do radial_point = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            radial_knots, radial_degree, &
                            radial_nodes(radial_point), radial_n, radial_n_d, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_bspline_basis( &
                            reduced_knots, radial_degree - 1, &
                            radial_nodes(radial_point), radial_d, &
                            ignored_derivative, &
                            local_status)
                        if (local_status /= 0) return
                        do radial = 1, radial_count - 1
                            radial_d(radial) = radial_d(radial)* &
                                real(radial_degree, dp)/( &
                                radial_knots(radial + radial_degree + 1) - &
                                radial_knots(radial + 1))
                        end do
                        call mapped_geometry_jacobian( &
                            radial_n, radial_n_d, angular_n, angular_n_d, &
                            radial_count, azimuth_count, control_points, &
                            weights, jacobian)
                        call invert_2d(jacobian, inverse, determinant, local_status)
                        if (local_status /= 0) return
                        quadrature_weight = radial_weights(radial_point)* &
                            angular_weights(angular_point)
                        if (form_degree == 1) then
                            reference_value = 0.0_dp
                            reference_curl = 0.0_dp
                            do radial = 1, radial_count - 1
                                do angular = 1, azimuth_count
                                    row = angular + (radial - 1)*azimuth_count
                                    reference_value(1, row) = &
                                        radial_d(radial)*angular_n(angular)
                                    reference_curl(row) = -radial_d(radial)* &
                                        angular_n_d(angular)
                                end do
                            end do
                            do radial = 1, radial_count
                                do angular = 1, azimuth_count
                                    row = tensor_radial_count + angular + &
                                        (radial - 1)*azimuth_count
                                    reference_value(2, row) = &
                                        radial_n(radial)*angular_d( &
                                        modulo(angular, azimuth_count) + 1)
                                    reference_curl(row) = &
                                        radial_n_d(radial)*angular_d( &
                                        modulo(angular, azimuth_count) + 1)
                                end do
                            end do
                            do column = 1, dof_count
                                physical_column = matmul( &
                                    transpose(inverse), &
                                    reference_value(:, column))
                                do row = 1, dof_count
                                    physical_row = matmul( &
                                        transpose(inverse), &
                                        reference_value(:, row))
                                    dense(row, column) = dense(row, column) + &
                                        quadrature_weight*( &
                                        mass_weight*determinant*dot_product( &
                                        physical_row, physical_column) + &
                                        differential_weight* &
                                        reference_curl(row)* &
                                        reference_curl(column)/determinant)
                                end do
                            end do
                        else
                            do radial = 1, radial_count - 1
                                do angular = 1, azimuth_count
                                    row = angular + (radial - 1)*azimuth_count
                                    scalar_value(row) = &
                                        radial_d(radial)*angular_d( &
                                        modulo(angular, azimuth_count) + 1)
                                end do
                            end do
                            do column = 1, dof_count
                                do row = 1, dof_count
                                    dense(row, column) = dense(row, column) + &
                                        quadrature_weight*mass_weight* &
                                        scalar_value(row)* &
                                        scalar_value(column)/determinant
                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do
        call dense_to_sparse(dense, matrix, status)
    end subroutine assemble_tensor_polar_form

    subroutine mapped_geometry_jacobian( &
            radial_value, radial_derivative, angular_value, &
            angular_derivative, radial_count, azimuth_count, control_points, &
            weights, jacobian)
        real(dp), intent(in) :: radial_value(:), radial_derivative(:)
        real(dp), intent(in) :: angular_value(:), angular_derivative(:)
        integer, intent(in) :: radial_count, azimuth_count
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        real(dp), intent(out) :: jacobian(2, 2)
        real(dp) :: denominator, denominator_derivative(2), factor
        real(dp) :: point(2), numerator(2), numerator_derivative(2, 2)
        integer :: angular, coordinate, radial

        denominator = 0.0_dp
        denominator_derivative = 0.0_dp
        numerator = 0.0_dp
        numerator_derivative = 0.0_dp
        do radial = 1, radial_count
            do angular = 1, azimuth_count
                factor = weights(angular, radial)*radial_value(radial)* &
                    angular_value(angular)
                denominator = denominator + factor
                numerator = numerator + factor*control_points(:, angular, radial)
                factor = weights(angular, radial)*radial_derivative(radial)* &
                    angular_value(angular)
                denominator_derivative(1) = denominator_derivative(1) + factor
                numerator_derivative(:, 1) = numerator_derivative(:, 1) + &
                    factor*control_points(:, angular, radial)
                factor = weights(angular, radial)*radial_value(radial)* &
                    angular_derivative(angular)
                denominator_derivative(2) = denominator_derivative(2) + factor
                numerator_derivative(:, 2) = numerator_derivative(:, 2) + &
                    factor*control_points(:, angular, radial)
            end do
        end do
        point = numerator/denominator
        do coordinate = 1, 2
            jacobian(:, coordinate) = ( &
                numerator_derivative(:, coordinate) - &
                point*denominator_derivative(coordinate))/denominator
        end do
    end subroutine mapped_geometry_jacobian

    subroutine invert_2d(jacobian, inverse, determinant, status)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp), intent(out) :: inverse(2, 2), determinant
        integer, intent(out) :: status

        determinant = jacobian(1, 1)*jacobian(2, 2) - &
            jacobian(1, 2)*jacobian(2, 1)
        status = 1
        inverse = 0.0_dp
        if (determinant <= 128.0_dp*epsilon(1.0_dp)) return
        inverse = reshape([ &
            jacobian(2, 2), -jacobian(2, 1), &
            -jacobian(1, 2), jacobian(1, 1)], [2, 2])/determinant
        status = 0
    end subroutine invert_2d

    subroutine dense_to_sparse(dense, matrix, status)
        real(dp), intent(in) :: dense(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        integer :: column, entry, nonzero_count, row

        nonzero_count = count(dense /= 0.0_dp)
        allocate(rows(nonzero_count), columns(nonzero_count), &
            values(nonzero_count))
        entry = 0
        do column = 1, size(dense, 2)
            do row = 1, size(dense, 1)
                if (dense(row, column) == 0.0_dp) cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = dense(row, column)
            end do
        end do
        call csc_from_triplet( &
            size(dense, 1), size(dense, 2), rows, columns, values, &
            matrix, status)
    end subroutine dense_to_sparse

    subroutine assemble_bspline_polar_h1_operator_csc( &
            radial_knots, radial_degree, azimuth_count, azimuth_degree, &
            control_points, weights, quadrature_order, matrix, status, &
            stiffness_coefficient, mass_coefficient)
        !! Physical H1 operator on a radial x periodic-angular polar patch.
        real(dp), intent(in) :: radial_knots(:)
        integer, intent(in) :: radial_degree, azimuth_count, azimuth_degree
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        integer, intent(in) :: quadrature_order
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: angular_derivative(:), angular_value(:)
        real(dp), allocatable :: extraction(:, :), local_matrix(:, :)
        real(dp), allocatable :: radial_derivative(:), radial_nodes(:)
        real(dp), allocatable :: radial_value(:), radial_weights(:)
        real(dp), allocatable :: triplet_values(:)
        real(dp) :: angular_nodes(quadrature_order)
        real(dp) :: angular_weights(quadrature_order)
        real(dp) :: determinant, gradient_column(2), gradient_row(2)
        real(dp) :: inverse(2, 2), jacobian(2, 2)
        real(dp) :: basis_column, basis_row, physical_weight
        real(dp) :: mass_weight, stiffness_weight
        integer :: angular_cell, angular_point, entry, local_column, local_count
        integer :: local_row, local_status, radial_count, radial_point, span
        type(csc_t) :: tensor_matrix

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Physical polar H1 assembly requires a regular mapped patch")
        radial_count = size(radial_knots) - radial_degree - 1
        if (radial_degree < 1 .or. azimuth_degree < 1 .or. &
            quadrature_order < 1) return
        if (radial_count < 5 .or. azimuth_count < azimuth_degree + 1) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= &
            [azimuth_count, radial_count])) return
        if (any(shape(weights) /= [azimuth_count, radial_count])) return
        if (any(weights <= 0.0_dp)) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) &
            stiffness_weight = stiffness_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (radial_degree + 1)*(azimuth_degree + 1)
        allocate( &
            rows(radial_span_count()*azimuth_count*local_count**2), &
            columns(radial_span_count()*azimuth_count*local_count**2), &
            triplet_values( &
            radial_span_count()*azimuth_count*local_count**2), &
            local_dofs(local_count), local_matrix(local_count, local_count), &
            radial_nodes(quadrature_order), radial_weights(quadrature_order))
        entry = 0
        do span = radial_degree + 1, radial_count
            if (radial_knots(span + 1) <= radial_knots(span)) cycle
            call gauss_legendre_ab( &
                quadrature_order, radial_knots(span), &
                radial_knots(span + 1), radial_nodes, radial_weights)
            do angular_cell = 1, azimuth_count
                call gauss_legendre_ab( &
                    quadrature_order, &
                    real(angular_cell - 1, dp)/real(azimuth_count, dp), &
                    real(angular_cell, dp)/real(azimuth_count, dp), &
                    angular_nodes, angular_weights)
                call local_dof_indices(span, angular_cell, local_dofs)
                local_matrix = 0.0_dp
                do angular_point = 1, quadrature_order
                    call evaluate_periodic_bspline_basis( &
                        azimuth_count, azimuth_degree, &
                        angular_nodes(angular_point), angular_value, &
                        angular_derivative, local_status)
                    if (local_status /= 0) return
                    do radial_point = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            radial_knots, radial_degree, &
                            radial_nodes(radial_point), radial_value, &
                            radial_derivative, local_status)
                        if (local_status /= 0) return
                        call geometry_jacobian( &
                            radial_value, radial_derivative, angular_value, &
                            angular_derivative, jacobian)
                        call inverse_jacobian( &
                            jacobian, inverse, determinant, local_status)
                        if (local_status /= 0) return
                        physical_weight = determinant* &
                            radial_weights(radial_point)* &
                            angular_weights(angular_point)
                        do local_column = 1, local_count
                            call basis_data( &
                                local_dofs(local_column), radial_value, &
                                radial_derivative, angular_value, &
                                angular_derivative, basis_column, &
                                gradient_column)
                            gradient_column = &
                                matmul(transpose(inverse), gradient_column)
                            do local_row = 1, local_count
                                call basis_data( &
                                    local_dofs(local_row), radial_value, &
                                    radial_derivative, angular_value, &
                                    angular_derivative, basis_row, gradient_row)
                                gradient_row = &
                                    matmul(transpose(inverse), gradient_row)
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*( &
                                    stiffness_weight*dot_product( &
                                    gradient_row, gradient_column) + &
                                    mass_weight*basis_row*basis_column)
                            end do
                        end do
                    end do
                end do
                call append_local_matrix( &
                    local_dofs, local_matrix, rows, columns, triplet_values, &
                    entry)
            end do
        end do
        call csc_from_triplet( &
            radial_count*azimuth_count, radial_count*azimuth_count, &
            rows(:entry), columns(:entry), triplet_values(:entry), &
            tensor_matrix, status)
        if (status%code /= 0) return
        call build_bspline_polar_h1_extraction( &
            azimuth_count, radial_count, extraction, local_status)
        if (local_status /= 0) return
        call restrict_bspline_polar_operator_csc( &
            extraction, tensor_matrix, matrix, status)

    contains

        integer function radial_span_count() result(count)
            integer :: radial_span

            count = 0
            do radial_span = radial_degree + 1, radial_count
                if (radial_knots(radial_span + 1) > &
                    radial_knots(radial_span)) count = count + 1
            end do
        end function radial_span_count

        subroutine local_dof_indices(radial_span, angular_span, dofs)
            integer, intent(in) :: radial_span, angular_span
            integer, intent(out) :: dofs(:)
            integer :: angular_basis, local, radial_basis

            local = 0
            do radial_basis = radial_span - radial_degree, radial_span
                do angular_basis = angular_span - azimuth_degree, angular_span
                    local = local + 1
                    dofs(local) = modulo(angular_basis - 1, azimuth_count) + &
                        1 + (radial_basis - 1)*azimuth_count
                end do
            end do
        end subroutine local_dof_indices

        subroutine basis_data( &
                dof, rv, rd, av, ad, value, gradient)
            integer, intent(in) :: dof
            real(dp), intent(in) :: rv(:), rd(:), av(:), ad(:)
            real(dp), intent(out) :: value, gradient(2)
            integer :: angular_basis, radial_basis

            angular_basis = modulo(dof - 1, azimuth_count) + 1
            radial_basis = (dof - 1)/azimuth_count + 1
            value = rv(radial_basis)*av(angular_basis)
            gradient = [ &
                rd(radial_basis)*av(angular_basis), &
                rv(radial_basis)*ad(angular_basis)]
        end subroutine basis_data

        subroutine geometry_jacobian(rv, rd, av, ad, jacobian)
            real(dp), intent(in) :: rv(:), rd(:), av(:), ad(:)
            real(dp), intent(out) :: jacobian(2, 2)
            real(dp) :: denominator, denominator_derivative(2), factor
            real(dp) :: point(2), numerator(2), numerator_derivative(2, 2)
            integer :: angular_basis, coordinate, radial_basis

            denominator = 0.0_dp
            denominator_derivative = 0.0_dp
            numerator = 0.0_dp
            numerator_derivative = 0.0_dp
            do radial_basis = 1, radial_count
                do angular_basis = 1, azimuth_count
                    factor = weights(angular_basis, radial_basis)* &
                        rv(radial_basis)*av(angular_basis)
                    denominator = denominator + factor
                    numerator = numerator + factor* &
                        control_points(:, angular_basis, radial_basis)
                    factor = weights(angular_basis, radial_basis)* &
                        rd(radial_basis)*av(angular_basis)
                    denominator_derivative(1) = &
                        denominator_derivative(1) + factor
                    numerator_derivative(:, 1) = &
                        numerator_derivative(:, 1) + factor* &
                        control_points(:, angular_basis, radial_basis)
                    factor = weights(angular_basis, radial_basis)* &
                        rv(radial_basis)*ad(angular_basis)
                    denominator_derivative(2) = &
                        denominator_derivative(2) + factor
                    numerator_derivative(:, 2) = &
                        numerator_derivative(:, 2) + factor* &
                        control_points(:, angular_basis, radial_basis)
                end do
            end do
            point = numerator/denominator
            do coordinate = 1, 2
                jacobian(:, coordinate) = ( &
                    numerator_derivative(:, coordinate) - &
                    point*denominator_derivative(coordinate))/denominator
            end do
        end subroutine geometry_jacobian

        subroutine inverse_jacobian( &
                jacobian, inverse, determinant, inverse_status)
            real(dp), intent(in) :: jacobian(2, 2)
            real(dp), intent(out) :: inverse(2, 2), determinant
            integer, intent(out) :: inverse_status

            determinant = jacobian(1, 1)*jacobian(2, 2) - &
                jacobian(1, 2)*jacobian(2, 1)
            inverse_status = 1
            inverse = 0.0_dp
            if (determinant <= 128.0_dp*epsilon(1.0_dp)) return
            inverse = reshape([ &
                jacobian(2, 2), -jacobian(2, 1), &
                -jacobian(1, 2), jacobian(1, 1)], [2, 2])/determinant
            inverse_status = 0
        end subroutine inverse_jacobian

        subroutine append_local_matrix( &
                dofs, local, row_indices, column_indices, values, count)
            integer, intent(in) :: dofs(:)
            real(dp), intent(in) :: local(:, :)
            integer, intent(inout) :: row_indices(:), column_indices(:), count
            real(dp), intent(inout) :: values(:)
            integer :: local_column, local_row

            do local_column = 1, size(dofs)
                do local_row = 1, size(dofs)
                    count = count + 1
                    row_indices(count) = dofs(local_row)
                    column_indices(count) = dofs(local_column)
                    values(count) = local(local_row, local_column)
                end do
            end do
        end subroutine append_local_matrix

    end subroutine assemble_bspline_polar_h1_operator_csc

    subroutine restrict_bspline_polar_operator_csc( &
            extraction, tensor_operator, polar_operator, status)
        !! Compute E A E^T for a tensor operator A and polar basis extraction E.
        real(dp), intent(in) :: extraction(:, :)
        type(csc_t), intent(in) :: tensor_operator
        type(csc_t), intent(out) :: polar_operator
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        type(csc_t) :: extraction_csc, extraction_transpose, work
        integer :: column, entry, row

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Polar Galerkin restriction requires compatible dimensions")
        if (size(extraction, 2) /= tensor_operator%nrow .or. &
            tensor_operator%nrow /= tensor_operator%ncol) return
        allocate( &
            rows(count(extraction /= 0.0_dp)), &
            columns(count(extraction /= 0.0_dp)), &
            values(count(extraction /= 0.0_dp)))
        entry = 0
        do column = 1, size(extraction, 2)
            do row = 1, size(extraction, 1)
                if (extraction(row, column) == 0.0_dp) cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = extraction(row, column)
            end do
        end do
        call csc_from_triplet( &
            size(extraction, 1), size(extraction, 2), rows, columns, values, &
            extraction_csc, status)
        if (status%code /= 0) return
        call csc_matmul(extraction_csc, tensor_operator, work, status)
        if (status%code /= 0) return
        call csc_transpose(extraction_csc, extraction_transpose, status)
        if (status%code /= 0) return
        call csc_matmul(work, extraction_transpose, polar_operator, status)
    end subroutine restrict_bspline_polar_operator_csc

end module fortfem_assembly_bspline_polar_2d
