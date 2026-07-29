program test_triangle_rt_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_triangle_raviart_thomas, &
        initialize_triangle_raviart_thomas, triangle_duffy_quadrature, &
        triangle_rt_basis_t, triangle_rt_dof_count
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    type(triangle_rt_basis_t) :: basis, copied_basis
    real(dp), allocatable :: divergences(:), dofs(:), values(:, :)
    real(dp) :: expected(2), expected_divergence, point(2)
    real(dp) :: reconstructed(2), reconstructed_divergence
    integer :: basis_dof, degree, dof_count, status
    logical :: all_passed

    all_passed = .true.
    point = [0.27_dp, 0.19_dp]
    do degree = 0, 3
        call initialize_triangle_raviart_thomas(degree, basis, status)
        dof_count = triangle_rt_dof_count(basis)
        call record_condition(status == 0 .and. &
            dof_count == (degree + 1) * (degree + 3), &
            "Arbitrary-order Raviart-Thomas basis has the exact dimension")
        allocate(values(2, dof_count), divergences(dof_count), dofs(dof_count))

        do basis_dof = 1, dof_count
            call check_basis_moments( &
                basis, degree, basis_dof, dof_count)
        end do

        call polynomial_rt_dofs(degree, dofs)
        call evaluate_triangle_raviart_thomas( &
            basis, point(1), point(2), values, divergences, status)
        reconstructed = matmul(values, dofs)
        reconstructed_divergence = dot_product(divergences, dofs)
        expected = [point(1)**(degree + 1), &
            point(1)**degree * point(2)]
        expected_divergence = real(degree + 2, dp) * point(1)**degree
        call record_condition(status == 0 .and. &
            maxval(abs(reconstructed - expected)) < 3.0e-10_dp, &
            "Raviart-Thomas interpolation reproduces a polynomial field")
        call record_condition(abs( &
            reconstructed_divergence - expected_divergence) < 5.0e-10_dp, &
            "Raviart-Thomas interpolation reproduces polynomial divergence")

        copied_basis = basis
        call initialize_triangle_raviart_thomas(0, basis, status)
        call evaluate_triangle_raviart_thomas( &
            copied_basis, point(1), point(2), values, divergences, status)
        call record_condition(status == 0, &
            "Assigned Raviart-Thomas basis retains an independent deep copy")
        deallocate(values, divergences, dofs)
    end do

    call initialize_triangle_raviart_thomas(-1, basis, status)
    call record_condition(status /= 0, &
        "Raviart-Thomas basis rejects a negative degree")

    call check_summary("Arbitrary-order triangle Raviart-Thomas basis")
    if (.not. all_passed) error stop 1

contains

    subroutine check_basis_moments(basis, degree, basis_dof, dof_count)
        type(triangle_rt_basis_t), intent(in) :: basis
        integer, intent(in) :: degree, basis_dof, dof_count

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: local_divergences(:), local_values(:, :)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), moment, normal(2), polynomial
        integer :: component, edge, exponent, moment_dof, node, status
        integer :: total_degree, x_degree, y_degree

        allocate(edge_nodes(degree + 3), edge_weights(degree + 3))
        allocate(local_values(2, dof_count), local_divergences(dof_count))
        call gauss_legendre_ab( &
            degree + 3, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        moment_dof = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment_dof = moment_dof + 1
                moment = 0.0_dp
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    call evaluate_triangle_raviart_thomas( &
                        basis, edge_point(1), edge_point(2), local_values, &
                        local_divergences, status)
                    polynomial = shifted_legendre(exponent, edge_nodes(node))
                    moment = moment + edge_weights(node) * polynomial * &
                        dot_product(local_values(:, basis_dof), normal)
                end do
                call record_condition(abs(moment - kronecker_delta( &
                    moment_dof, basis_dof)) < 4.0e-10_dp, &
                    "Raviart-Thomas edge flux moments form a Kronecker basis")
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * degree + 2, xi, eta, triangle_weights, status)
        do component = 1, 2
            do total_degree = 0, degree - 1
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment_dof = moment_dof + 1
                    moment = 0.0_dp
                    do node = 1, size(xi)
                        call evaluate_triangle_raviart_thomas( &
                            basis, xi(node), eta(node), local_values, &
                            local_divergences, status)
                        moment = moment + triangle_weights(node) * &
                            local_values(component, basis_dof) * &
                            xi(node)**x_degree * eta(node)**y_degree
                    end do
                    call record_condition(abs(moment - kronecker_delta( &
                        moment_dof, basis_dof)) < 4.0e-10_dp, &
                        "Raviart-Thomas cell moments form a Kronecker basis")
                end do
            end do
        end do
    end subroutine check_basis_moments

    subroutine polynomial_rt_dofs(degree, dofs)
        integer, intent(in) :: degree
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), field(2), normal(2)
        integer :: component, edge, exponent, moment_dof, node, status
        integer :: total_degree, x_degree, y_degree

        allocate(edge_nodes(degree + 3), edge_weights(degree + 3))
        call gauss_legendre_ab( &
            degree + 3, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment_dof = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment_dof = moment_dof + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    field = [edge_point(1)**(degree + 1), &
                        edge_point(1)**degree * edge_point(2)]
                    dofs(moment_dof) = dofs(moment_dof) + &
                        edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(field, normal)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * degree + 2, xi, eta, triangle_weights, status)
        do component = 1, 2
            do total_degree = 0, degree - 1
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment_dof = moment_dof + 1
                    if (component == 1) then
                        dofs(moment_dof) = sum(triangle_weights * &
                            xi**(degree + 1 + x_degree) * eta**y_degree)
                    else
                        dofs(moment_dof) = sum(triangle_weights * &
                            xi**(degree + x_degree) * eta**(y_degree + 1))
                    end if
                end do
            end do
        end do
    end subroutine polynomial_rt_dofs

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

    pure function kronecker_delta(first, second) result(value)
        integer, intent(in) :: first, second
        real(dp) :: value

        if (first == second) then
            value = 1.0_dp
        else
            value = 0.0_dp
        end if
    end function kronecker_delta

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_rt_arbitrary_order
