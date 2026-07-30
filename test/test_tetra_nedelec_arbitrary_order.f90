program test_tetra_nedelec_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_nedelec_first_kind, &
        evaluate_tetra_nedelec_first_order, &
        initialize_tetra_nedelec_first_kind, &
        interpolate_reference_tetra_nedelec, tetra_duffy_quadrature, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t, &
        triangle_duffy_quadrature
    use fortfem_generated_tetra_nedelec_coefficients, only: &
        load_tetra_nedelec_coefficients
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    type(tetra_nedelec_first_kind_t) :: basis, copied_basis
    real(dp), allocatable :: curls(:, :), dofs(:), moment_matrix(:, :)
    real(dp), allocatable :: generated_coefficients(:, :), values(:, :)
    real(dp) :: expected(3), p_errors(5), point(3), reconstructed(3)
    integer :: dof_count, order, status
    logical :: all_passed

    all_passed = .true.
    point = [0.23_dp, 0.19_dp, 0.17_dp]
    do order = 1, 5
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        dof_count = tetra_nedelec_dof_count(basis)
        call record_condition(status == 0 .and. dof_count == &
            order * (order + 2) * (order + 3) / 2, &
            "Tetrahedral first-kind Nedelec basis has the exact dimension")
        allocate( &
            values(3, dof_count), curls(3, dof_count), dofs(dof_count), &
            moment_matrix(dof_count, dof_count))
        if (order <= 4) then
            call load_tetra_nedelec_coefficients( &
                order, generated_coefficients, status)
            call record_condition(status == 0 .and. &
                size(generated_coefficients, 1) == dof_count .and. &
                size(generated_coefficients, 2) == dof_count, &
                "Generated tetrahedral basis coefficients have the exact shape")
        end if

        call build_moment_matrix(basis, order, moment_matrix)
        call record_condition(maxval(abs(moment_matrix - identity( &
            dof_count))) < 2.0e-8_dp, &
            "Tetrahedral edge, face, and cell moments form a Kronecker basis")

        call polynomial_gradient_dofs(order, order, dofs)
        call evaluate_tetra_nedelec_first_kind( &
            basis, point, values, curls, status)
        reconstructed = matmul(values, dofs)
        expected = [real(order, dp) * point(1)**(order - 1), 0.0_dp, 0.0_dp]
        call record_condition(maxval(abs(reconstructed - expected)) < &
            2.0e-8_dp, &
            "Tetrahedral Nedelec interpolation reproduces a polynomial gradient")
        call record_condition(maxval(abs(matmul(curls, dofs))) < 2.0e-8_dp, &
            "Tetrahedral Nedelec interpolation preserves zero gradient curl")
        call check_curls_by_finite_difference(basis, dof_count, point)
        call polynomial_gradient_interpolation_error( &
            basis, 4, p_errors(order))
        if (order == 1) call check_lowest_order_compatibility( &
            basis, values, curls)

        copied_basis = basis
        call initialize_tetra_nedelec_first_kind(1, basis, status)
        call evaluate_tetra_nedelec_first_kind( &
            copied_basis, point, values, curls, status)
        call record_condition(status == 0, &
            "Assigned tetrahedral Nedelec basis retains an independent copy")

        deallocate(values, curls, dofs, moment_matrix)
        if (allocated(generated_coefficients)) &
            deallocate(generated_coefficients)
    end do
    call record_condition(all(p_errors(2:4) < p_errors(1:3)), &
        "Tetrahedral Nedelec interpolation improves at every higher order")
    call record_condition(p_errors(4) < 2.0e-11_dp, &
        "Order-four tetrahedral Nedelec interpolation reproduces a cubic field")

    call initialize_tetra_nedelec_first_kind(0, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral first-kind Nedelec basis rejects order zero")
    call initialize_tetra_nedelec_first_kind(6, basis, status)
    call record_condition(status /= 0, &
        "Tetrahedral first-kind Nedelec basis rejects unsupported order six")

    call check_summary("Arbitrary-order tetrahedral first-kind Nedelec basis")
    if (.not. all_passed) error stop 1

contains

    subroutine build_moment_matrix(basis, order, matrix)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        integer, intent(in) :: order
        real(dp), intent(out) :: matrix(:, :)

        real(dp), allocatable :: curls(:, :), edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: tetra_weights(:), triangle_weights(:)
        real(dp), allocatable :: values(:, :), x(:), y(:), z(:)
        real(dp) :: point(3), tangent(3), tangents(3, 2)
        integer :: basis_dof, component, edge, exponent, face
        integer :: moment, node, status, total_degree
        integer :: x_degree, y_degree, z_degree

        allocate( &
            values(3, size(matrix, 2)), curls(3, size(matrix, 2)), &
            edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        matrix = 0.0_dp
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    call evaluate_tetra_nedelec_first_kind( &
                        basis, point, values, curls, status)
                    do basis_dof = 1, size(matrix, 2)
                        matrix(moment, basis_dof) = &
                            matrix(moment, basis_dof) + &
                            edge_weights(node) * &
                            shifted_legendre(exponent, edge_nodes(node)) * &
                            dot_product(values(:, basis_dof), tangent)
                    end do
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order + 2, x, y, triangle_weights, status)
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            call evaluate_tetra_nedelec_first_kind( &
                                basis, point, values, curls, status)
                            do basis_dof = 1, size(matrix, 2)
                                matrix(moment, basis_dof) = &
                                    matrix(moment, basis_dof) + &
                                    triangle_weights(node) * &
                                    x(node)**x_degree * y(node)**y_degree * &
                                    dot_product( &
                                    values(:, basis_dof), &
                                    tangents(:, component))
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2 * order + 2, x, y, z, tetra_weights, status)
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            point = [x(node), y(node), z(node)]
                            call evaluate_tetra_nedelec_first_kind( &
                                basis, point, values, curls, status)
                            matrix(moment, :) = matrix(moment, :) + &
                                tetra_weights(node) * x(node)**x_degree * &
                                y(node)**y_degree * z(node)**z_degree * &
                                values(component, :)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine build_moment_matrix

    subroutine polynomial_gradient_dofs(order, field_degree, dofs)
        integer, intent(in) :: order, field_degree
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: tetra_weights(:), triangle_weights(:)
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: field(3), point(3), tangent(3), tangents(3, 2)
        integer :: component, edge, exponent, face, moment, node, status
        integer :: total_degree, x_degree, y_degree, z_degree

        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    field = gradient_field(field_degree, point)
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(field, tangent)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order + 2, x, y, triangle_weights, status)
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            field = gradient_field(field_degree, point)
                            dofs(moment) = dofs(moment) + &
                                triangle_weights(node) * &
                                x(node)**x_degree * y(node)**y_degree * &
                                dot_product(field, tangents(:, component))
                        end do
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2 * order + 2, x, y, z, tetra_weights, status)
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            point = [x(node), y(node), z(node)]
                            field = gradient_field(field_degree, point)
                            dofs(moment) = dofs(moment) + &
                                tetra_weights(node) * x(node)**x_degree * &
                                y(node)**y_degree * z(node)**z_degree * &
                                field(component)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine polynomial_gradient_dofs

    subroutine polynomial_gradient_interpolation_error( &
            basis, field_degree, error)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        integer, intent(in) :: field_degree
        real(dp), intent(out) :: error

        real(dp), allocatable :: curls(:, :), dofs(:), values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: difference(3), location(3)
        integer :: node, status

        allocate( &
            dofs(tetra_nedelec_dof_count(basis)), &
            values(3, tetra_nedelec_dof_count(basis)), &
            curls(3, tetra_nedelec_dof_count(basis)))
        call interpolate_reference_tetra_nedelec( &
            basis, cubic_gradient, dofs, status)
        call tetra_duffy_quadrature(12, x, y, z, weights, status)
        error = 0.0_dp
        do node = 1, size(weights)
            location = [x(node), y(node), z(node)]
            call evaluate_tetra_nedelec_first_kind( &
                basis, location, values, curls, status)
            difference = matmul(values, dofs) - &
                gradient_field(field_degree, location)
            error = error + weights(node) * dot_product(difference, difference)
        end do
        error = sqrt(error)
    end subroutine polynomial_gradient_interpolation_error

    pure subroutine cubic_gradient(location, field)
        real(dp), intent(in) :: location(3)
        real(dp), intent(out) :: field(3)

        field = gradient_field(4, location)
    end subroutine cubic_gradient

    subroutine check_curls_by_finite_difference( &
            basis, dof_count, point)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        integer, intent(in) :: dof_count
        real(dp), intent(in) :: point(3)

        real(dp), parameter :: step = 2.0e-6_dp
        real(dp), allocatable :: curls(:, :), scratch_curls(:, :)
        real(dp), allocatable :: values_minus(:, :), values_plus(:, :)
        real(dp) :: curl_error, maximum_error, numerical(3), shifted(3)
        integer :: basis_dof, direction, status

        allocate( &
            curls(3, dof_count), scratch_curls(3, dof_count), &
            values_minus(3, dof_count), values_plus(3, dof_count))
        call evaluate_tetra_nedelec_first_kind( &
            basis, point, values_plus, curls, status)
        maximum_error = 0.0_dp
        do basis_dof = 1, dof_count
            numerical = 0.0_dp
            do direction = 1, 3
                shifted = point
                shifted(direction) = shifted(direction) + step
                call evaluate_tetra_nedelec_first_kind( &
                    basis, shifted, values_plus, scratch_curls, status)
                shifted(direction) = shifted(direction) - 2.0_dp * step
                call evaluate_tetra_nedelec_first_kind( &
                    basis, shifted, values_minus, scratch_curls, status)
                select case (direction)
                case (1)
                    numerical(2) = numerical(2) - &
                        (values_plus(3, basis_dof) - &
                        values_minus(3, basis_dof)) / (2.0_dp * step)
                    numerical(3) = numerical(3) + &
                        (values_plus(2, basis_dof) - &
                        values_minus(2, basis_dof)) / (2.0_dp * step)
                case (2)
                    numerical(1) = numerical(1) + &
                        (values_plus(3, basis_dof) - &
                        values_minus(3, basis_dof)) / (2.0_dp * step)
                    numerical(3) = numerical(3) - &
                        (values_plus(1, basis_dof) - &
                        values_minus(1, basis_dof)) / (2.0_dp * step)
                case (3)
                    numerical(1) = numerical(1) - &
                        (values_plus(2, basis_dof) - &
                        values_minus(2, basis_dof)) / (2.0_dp * step)
                    numerical(2) = numerical(2) + &
                        (values_plus(1, basis_dof) - &
                        values_minus(1, basis_dof)) / (2.0_dp * step)
                end select
            end do
            curl_error = maxval(abs(numerical - curls(:, basis_dof)))/ &
                max(1.0_dp, maxval(abs(curls(:, basis_dof))))
            maximum_error = max(maximum_error, curl_error)
        end do
        call record_condition(maximum_error < 3.0e-7_dp, &
            "Tetrahedral basis curls agree with finite differences")
    end subroutine check_curls_by_finite_difference

    subroutine check_lowest_order_compatibility(basis, values, curls)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(out) :: values(:, :), curls(:, :)

        real(dp) :: legacy_curls(3, 6), legacy_values(3, 6)
        integer :: status

        call evaluate_tetra_nedelec_first_kind( &
            basis, point, values, curls, status)
        call evaluate_tetra_nedelec_first_order( &
            point, legacy_values, legacy_curls, status)
        call record_condition(maxval(abs(values - legacy_values)) < &
            2.0e-13_dp .and. maxval(abs(curls - legacy_curls)) < 2.0e-13_dp, &
            "Arbitrary-order tetrahedral basis preserves lowest-order values")
    end subroutine check_lowest_order_compatibility

    pure function gradient_field(order, point) result(field)
        integer, intent(in) :: order
        real(dp), intent(in) :: point(3)
        real(dp) :: field(3)

        field = [real(order, dp) * point(1)**(order - 1), 0.0_dp, 0.0_dp]
    end function gradient_field

    pure function identity(size_) result(matrix)
        integer, intent(in) :: size_
        real(dp) :: matrix(size_, size_)

        integer :: diagonal

        matrix = 0.0_dp
        do diagonal = 1, size_
            matrix(diagonal, diagonal) = 1.0_dp
        end do
    end function identity

    pure subroutine reference_edge(edge, parameter, point, tangent)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(3), tangent(3)

        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), first, second

        call reference_topology(vertices, edge_vertices)
        first = edge_vertices(1, edge)
        second = edge_vertices(2, edge)
        point = (1.0_dp - parameter) * vertices(:, first) + &
            parameter * vertices(:, second)
        tangent = vertices(:, second) - vertices(:, first)
    end subroutine reference_edge

    pure subroutine reference_face(face, u, v, point, tangents)
        integer, intent(in) :: face
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: point(3), tangents(3, 2)

        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), face_vertices(3, 4)
        integer :: first, second, third

        call reference_topology(vertices, edge_vertices)
        face_vertices(:, 1) = [1, 2, 3]
        face_vertices(:, 2) = [1, 2, 4]
        face_vertices(:, 3) = [1, 3, 4]
        face_vertices(:, 4) = [2, 3, 4]
        first = face_vertices(1, face)
        second = face_vertices(2, face)
        third = face_vertices(3, face)
        tangents(:, 1) = vertices(:, second) - vertices(:, first)
        tangents(:, 2) = vertices(:, third) - vertices(:, first)
        point = vertices(:, first) + u * tangents(:, 1) + &
            v * tangents(:, 2)
    end subroutine reference_face

    pure subroutine reference_topology(vertices, edge_vertices)
        real(dp), intent(out) :: vertices(3, 4)
        integer, intent(out) :: edge_vertices(2, 6)

        vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        edge_vertices(:, 1) = [1, 2]
        edge_vertices(:, 2) = [1, 3]
        edge_vertices(:, 3) = [1, 4]
        edge_vertices(:, 4) = [2, 3]
        edge_vertices(:, 5) = [2, 4]
        edge_vertices(:, 6) = [3, 4]
    end subroutine reference_topology

    pure function shifted_legendre(degree, parameter) result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: parameter
        real(dp) :: value

        real(dp) :: current, previous, coordinate
        integer :: polynomial_degree

        coordinate = 2.0_dp * parameter - 1.0_dp
        if (degree == 0) then
            value = 1.0_dp
            return
        end if
        previous = 1.0_dp
        current = coordinate
        do polynomial_degree = 1, degree - 1
            value = ( &
                real(2 * polynomial_degree + 1, dp) * coordinate * current - &
                real(polynomial_degree, dp) * previous) / &
                real(polynomial_degree + 1, dp)
            previous = current
            current = value
        end do
        value = current
    end function shifted_legendre

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_arbitrary_order
