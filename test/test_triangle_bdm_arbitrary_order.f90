program test_triangle_bdm_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_triangle_bdm, initialize_triangle_bdm, &
        triangle_bdm_basis_t, triangle_bdm_dof_count, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    type(triangle_bdm_basis_t) :: basis, copied_basis
    real(dp), allocatable :: divergences(:), dofs(:), values(:, :)
    real(dp) :: expected(2), expected_divergence, point(2)
    real(dp) :: reconstructed(2), reconstructed_divergence
    integer :: basis_dof, degree, dof_count, status
    logical :: all_passed

    all_passed = .true.
    point = [0.27_dp, 0.19_dp]
    do degree = 1, 4
        call initialize_triangle_bdm(degree, basis, status)
        dof_count = triangle_bdm_dof_count(basis)
        call record_condition(status == 0 .and. &
            dof_count == (degree + 1) * (degree + 2), &
            "BDM basis has the exact vector polynomial dimension")
        allocate(values(2, dof_count), divergences(dof_count), dofs(dof_count))

        do basis_dof = 1, dof_count
            call check_basis_moments( &
                basis, degree, basis_dof, dof_count)
        end do
        call polynomial_dofs(degree, dofs)
        call evaluate_triangle_bdm( &
            basis, point(1), point(2), values, divergences, status)
        reconstructed = matmul(values, dofs)
        reconstructed_divergence = dot_product(divergences, dofs)
        expected = [point(1)**degree + point(2), &
            point(2)**degree - point(1)]
        expected_divergence = real(degree, dp) * ( &
            point(1)**(degree - 1) + point(2)**(degree - 1))
        call record_condition(maxval(abs(reconstructed - expected)) < &
            2.0e-9_dp, "BDM interpolation reproduces [P_k]^2 fields")
        call record_condition(abs( &
            reconstructed_divergence - expected_divergence) < 3.0e-9_dp, &
            "BDM interpolation reproduces exact polynomial divergence")

        copied_basis = basis
        call initialize_triangle_bdm(1, basis, status)
        call evaluate_triangle_bdm( &
            copied_basis, point(1), point(2), values, divergences, status)
        call record_condition(status == 0, &
            "Assigned BDM basis retains an independent deep copy")
        deallocate(values, divergences, dofs)
    end do

    call initialize_triangle_bdm(0, basis, status)
    call record_condition(status /= 0, "BDM basis rejects degree zero")

    call check_summary("Triangle BDM arbitrary-order basis")
    if (.not. all_passed) error stop 1

contains

    subroutine check_basis_moments(basis, degree, basis_dof, dof_count)
        use fortfem_feec, only: &
            evaluate_triangle_raviart_thomas, &
            initialize_triangle_raviart_thomas, triangle_rt_basis_t, &
            triangle_rt_dof_count
        type(triangle_bdm_basis_t), intent(in) :: basis
        integer, intent(in) :: degree, basis_dof, dof_count

        type(triangle_rt_basis_t) :: rt_basis
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: local_divergences(:), local_values(:, :)
        real(dp), allocatable :: rt_divergences(:), rt_values(:, :)
        real(dp), allocatable :: eta(:), triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), moment, normal(2), polynomial
        real(dp) :: rotated_test(2)
        integer :: edge, exponent, moment_dof, node, rt_dof, rt_dof_count
        integer :: status

        allocate(edge_nodes(degree + 2), edge_weights(degree + 2))
        allocate(local_values(2, dof_count), local_divergences(dof_count))
        call gauss_legendre_ab( &
            degree + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        moment_dof = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment_dof = moment_dof + 1
                moment = 0.0_dp
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    call evaluate_triangle_bdm( &
                        basis, edge_point(1), edge_point(2), local_values, &
                        local_divergences, status)
                    polynomial = shifted_legendre(exponent, edge_nodes(node))
                    moment = moment + edge_weights(node) * polynomial * &
                        dot_product(local_values(:, basis_dof), normal)
                end do
                call record_condition(abs(moment - kronecker_delta( &
                    moment_dof, basis_dof)) < 2.0e-8_dp, &
                    "BDM normal edge moments form a Kronecker basis")
            end do
        end do

        if (degree < 2) return
        call initialize_triangle_raviart_thomas( &
            degree - 2, rt_basis, status)
        rt_dof_count = triangle_rt_dof_count(rt_basis)
        allocate(rt_values(2, rt_dof_count), rt_divergences(rt_dof_count))
        call triangle_duffy_quadrature( &
            2 * degree, xi, eta, triangle_weights, status)
        do rt_dof = 1, rt_dof_count
            moment_dof = moment_dof + 1
            moment = 0.0_dp
            do node = 1, size(xi)
                call evaluate_triangle_bdm( &
                    basis, xi(node), eta(node), local_values, &
                    local_divergences, status)
                call evaluate_triangle_raviart_thomas( &
                    rt_basis, xi(node), eta(node), rt_values, rt_divergences, &
                    status)
                call rotate_clockwise(rt_values(:, rt_dof), rotated_test)
                moment = moment + triangle_weights(node) * &
                    dot_product(local_values(:, basis_dof), rotated_test)
            end do
            call record_condition(abs(moment - kronecker_delta( &
                moment_dof, basis_dof)) < 3.0e-8_dp, &
                "BDM interior rotated-RT moments form a Kronecker basis")
        end do
    end subroutine check_basis_moments

    subroutine polynomial_dofs(degree, dofs)
        use fortfem_feec, only: &
            evaluate_triangle_raviart_thomas, &
            initialize_triangle_raviart_thomas, triangle_rt_basis_t, &
            triangle_rt_dof_count
        integer, intent(in) :: degree
        real(dp), intent(out) :: dofs(:)

        type(triangle_rt_basis_t) :: rt_basis
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: eta(:), rt_divergences(:), rt_values(:, :)
        real(dp), allocatable :: triangle_weights(:), xi(:)
        real(dp) :: edge_point(2), field(2), normal(2), polynomial
        real(dp) :: rotated_test(2)
        integer :: edge, exponent, moment, node, rt_dof, rt_dof_count, status

        allocate(edge_nodes(degree + 2), edge_weights(degree + 2))
        call gauss_legendre_ab( &
            degree + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment = 0
        do edge = 1, 3
            do exponent = 0, degree
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), edge_point, normal)
                    call polynomial_field(degree, edge_point, field)
                    polynomial = shifted_legendre(exponent, edge_nodes(node))
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        polynomial * dot_product(field, normal)
                end do
            end do
        end do

        if (degree < 2) return
        call initialize_triangle_raviart_thomas( &
            degree - 2, rt_basis, status)
        rt_dof_count = triangle_rt_dof_count(rt_basis)
        allocate(rt_values(2, rt_dof_count), rt_divergences(rt_dof_count))
        call triangle_duffy_quadrature( &
            2 * degree, xi, eta, triangle_weights, status)
        do rt_dof = 1, rt_dof_count
            moment = moment + 1
            do node = 1, size(xi)
                call polynomial_field( &
                    degree, [xi(node), eta(node)], field)
                call evaluate_triangle_raviart_thomas( &
                    rt_basis, xi(node), eta(node), rt_values, rt_divergences, &
                    status)
                call rotate_clockwise(rt_values(:, rt_dof), rotated_test)
                dofs(moment) = dofs(moment) + triangle_weights(node) * &
                    dot_product(field, rotated_test)
            end do
        end do
    end subroutine polynomial_dofs

    pure subroutine polynomial_field(degree, point, value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: value(2)

        value = [point(1)**degree + point(2), &
            point(2)**degree - point(1)]
    end subroutine polynomial_field

    pure subroutine rotate_clockwise(input, output)
        real(dp), intent(in) :: input(2)
        real(dp), intent(out) :: output(2)

        output = [input(2), -input(1)]
    end subroutine rotate_clockwise

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

end program test_triangle_bdm_arbitrary_order
