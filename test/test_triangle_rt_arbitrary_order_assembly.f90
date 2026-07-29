program test_triangle_rt_arbitrary_order_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_triangle_rt_div_mass_element, triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    real(dp), allocatable :: dofs(:), matrix(:, :), matrix_times_dofs(:)
    real(dp) :: exact_energy, reference_vertices(2, 3)
    real(dp) :: stretched_vertices(2, 3)
    integer :: degree, status
    logical :: all_passed

    all_passed = .true.
    reference_vertices(:, 1) = [0.0_dp, 0.0_dp]
    reference_vertices(:, 2) = [1.0_dp, 0.0_dp]
    reference_vertices(:, 3) = [0.0_dp, 1.0_dp]
    stretched_vertices(:, 1) = [0.0_dp, 0.0_dp]
    stretched_vertices(:, 2) = [2.0_dp, 0.0_dp]
    stretched_vertices(:, 3) = [0.0_dp, 3.0_dp]

    do degree = 0, 3
        allocate(dofs((degree + 1) * (degree + 3)))
        call polynomial_rt_dofs(degree, dofs)
        call assemble_triangle_rt_div_mass_element( &
            reference_vertices, degree, 2 * degree + 2, matrix, status)
        allocate(matrix_times_dofs(size(dofs)))
        matrix_times_dofs = matmul(matrix, dofs)
        exact_energy = exact_polynomial_energy(degree)
        call record_condition(status == 0 .and. &
            size(matrix, 1) == (degree + 1) * (degree + 3), &
            "Arbitrary-degree RT element matrix has the exact dimension")
        call record_condition(maxval(abs(matrix - transpose(matrix))) < &
            2.0e-12_dp, "Arbitrary-degree RT element matrix is symmetric")
        call record_condition(abs( &
            dot_product(dofs, matrix_times_dofs) - exact_energy) < &
            3.0e-8_dp, &
            "RT element integrates exact polynomial div-plus-mass energy")
        deallocate(dofs, matrix, matrix_times_dofs)
    end do

    allocate(dofs(3))
    dofs = [0.0_dp, 3.0_dp, -3.0_dp]
    call assemble_triangle_rt_div_mass_element( &
        stretched_vertices, 0, 2, matrix, status)
    allocate(matrix_times_dofs(3))
    matrix_times_dofs = matmul(matrix, dofs)
    call record_condition(abs( &
        dot_product(dofs, matrix_times_dofs) - 3.0_dp) < 1.0e-12_dp, &
        "Contravariant Piola assembly preserves constant physical flux energy")
    deallocate(dofs, matrix, matrix_times_dofs)

    call assemble_triangle_rt_div_mass_element( &
        reference_vertices, -1, 2, matrix, status)
    call record_condition(status /= 0, &
        "Arbitrary-degree RT element assembly rejects negative degree")

    call check_summary("Arbitrary-order triangle RT element assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine polynomial_rt_dofs(degree, dofs)
        integer, intent(in) :: degree
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), field(2), normal(2)
        integer :: component, edge, exponent, moment, node, quadrature_status
        integer :: total_degree, x_degree, y_degree

        allocate(edge_nodes(degree + 3), edge_weights(degree + 3))
        call gauss_legendre_ab( &
            degree + 3, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    field = [edge_point(1)**(degree + 1), &
                        edge_point(1)**degree * edge_point(2)]
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(field, normal)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * degree + 2, xi, eta, triangle_weights, quadrature_status)
        do component = 1, 2
            do total_degree = 0, degree - 1
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment = moment + 1
                    if (component == 1) then
                        dofs(moment) = sum(triangle_weights * &
                            xi**(degree + 1 + x_degree) * eta**y_degree)
                    else
                        dofs(moment) = sum(triangle_weights * &
                            xi**(degree + x_degree) * eta**(y_degree + 1))
                    end if
                end do
            end do
        end do
    end subroutine polynomial_rt_dofs

    pure function exact_polynomial_energy(degree) result(value)
        integer, intent(in) :: degree
        real(dp) :: value

        real(dp) :: divergence_energy, first_mass, second_mass

        divergence_energy = real((degree + 2)**2, dp) / &
            real((2 * degree + 1) * (2 * degree + 2), dp)
        first_mass = 1.0_dp / &
            real((2 * degree + 3) * (2 * degree + 4), dp)
        second_mass = 2.0_dp / real( &
            (2 * degree + 1) * (2 * degree + 2) * &
            (2 * degree + 3) * (2 * degree + 4), dp)
        value = divergence_energy + first_mass + second_mass
    end function exact_polynomial_energy

    pure subroutine reference_edge(edge, parameter, point, normal)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(2), normal(2)

        select case (edge)
        case (1)
            point = [parameter, 0.0_dp]
            normal = [0.0_dp, -1.0_dp]
        case (2)
            point = [1.0_dp - parameter, parameter]
            normal = [1.0_dp, 1.0_dp]
        case (3)
            point = [0.0_dp, 1.0_dp - parameter]
            normal = [-1.0_dp, 0.0_dp]
        end select
    end subroutine reference_edge

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

end program test_triangle_rt_arbitrary_order_assembly
