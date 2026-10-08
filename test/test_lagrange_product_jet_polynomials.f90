program test_lagrange_product_jet_polynomials
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_triangle_lagrange_arbitrary_order, only: &
        triangle_lagrange_basis_t, initialize_triangle_lagrange_basis, &
        triangle_lagrange_nodes, evaluate_triangle_lagrange_basis
    use fortfem_tetra_lagrange_arbitrary_order, only: tetra_lagrange_t, &
        initialize_tetra_lagrange, tetra_lagrange_nodes, &
        evaluate_tetra_lagrange, evaluate_tetra_lagrange_jvp
    implicit none
    type(triangle_lagrange_basis_t) :: triangle
    type(tetra_lagrange_t) :: tetra
    real(dp), allocatable :: nodes(:, :), values(:), gradients(:, :)
    real(dp), allocatable :: coefficients(:), values_dot(:), gradients_dot(:, :)
    real(dp) :: point(3), direction(3), exact, gradient(3), hessian(3, 3)
    real(dp) :: dummy_gradient(3), dummy_hessian(3, 3), sample
    real(dp) :: largest_scaled_error = 0.0_dp
    integer :: degree, status, n, i, probe
    direction = [0.3_dp, -0.2_dp, 0.4_dp]
    do degree = 0, 7
        call initialize_triangle_lagrange_basis(degree, triangle, status)
        if (status /= 0) error stop 'triangle initialization'
        call triangle_lagrange_nodes(triangle, nodes)
        n = size(nodes, 2)
        allocate(values(n), gradients(2, n), coefficients(n))
        do i = 1, n
            point = [nodes(:, i), 0.0_dp]
            call polynomial(point, degree, 2, coefficients(i), &
                dummy_gradient, dummy_hessian)
        end do
        do probe = 1, 3
            point = probe_point(probe)
            call evaluate_triangle_lagrange_basis(triangle, point(1), &
                point(2), values, gradients, status)
            if (status /= 0) error stop 'triangle evaluation'
            call polynomial(point, degree, 2, exact, gradient, hessian)
            call close(dot_product(coefficients, values), exact, 'triangle value')
            do i = 1, 2
                call close(dot_product(coefficients, gradients(i, :)), &
                    gradient(i), 'triangle gradient')
            end do
            call close(sum(values), 1.0_dp, 'triangle partition')
        end do
        deallocate(nodes, values, gradients, coefficients)
        call initialize_tetra_lagrange(degree, tetra, status)
        if (status /= 0) error stop 'tetra initialization'
        call tetra_lagrange_nodes(tetra, nodes)
        n = size(nodes, 2)
        allocate(values(n), gradients(3, n), coefficients(n), &
            values_dot(n), gradients_dot(3, n))
        do i = 1, n
            call polynomial(nodes(:, i), degree, 3, coefficients(i), &
                dummy_gradient, dummy_hessian)
        end do
        do probe = 1, 3
            point = probe_point(probe)
            call evaluate_tetra_lagrange(tetra, point, values, gradients, status)
            if (status /= 0) error stop 'tetra evaluation'
            call evaluate_tetra_lagrange_jvp(tetra, point, direction, &
                values_dot, gradients_dot, status)
            if (status /= 0) error stop 'tetra JVP'
            call polynomial(point, degree, 3, exact, gradient, hessian)
            call close(dot_product(coefficients, values), exact, 'tetra value')
            call close(dot_product(coefficients, values_dot), &
                dot_product(gradient, direction), 'tetra value JVP')
            do i = 1, 3
                call close(dot_product(coefficients, gradients(i, :)), &
                    gradient(i), 'tetra gradient')
                sample = dot_product(hessian(i, :), direction)
                call close(dot_product(coefficients, gradients_dot(i, :)), &
                    sample, 'tetra gradient JVP')
            end do
            call close(sum(values), 1.0_dp, 'tetra partition')
        end do
        deallocate(nodes, values, gradients, coefficients, values_dot, gradients_dot)
    end do
    print *, 'PASS: degrees 0..7 polynomial value/gradient/JVP, interior/face/vertex'
    print *, 'Maximum scaled polynomial error:', largest_scaled_error
    call check_summary('Lagrange product-jet polynomial reproduction')
contains
    function probe_point(probe) result(point)
        integer, intent(in) :: probe
        real(dp) :: point(3)
        select case (probe)
        case (1)
            point = [0.17_dp, 0.23_dp, 0.19_dp]
        case (2)
            point = [0.0_dp, 0.37_dp, 0.0_dp]
        case default
            point = [1.0_dp, 0.0_dp, 0.0_dp]
        end select
    end function probe_point

    subroutine polynomial(point, degree, dimension, value, gradient, hessian)
        real(dp), intent(in) :: point(3)
        integer, intent(in) :: degree, dimension
        real(dp), intent(out) :: value, gradient(3), hessian(3, 3)
        real(dp), parameter :: scale(3) = [1.0_dp, 2.0_dp, -3.0_dp]
        integer :: axis
        value = 1.0_dp
        gradient = 0.0_dp
        hessian = 0.0_dp
        if (degree == 0) return
        value = 0.0_dp
        do axis = 1, dimension
            value = value + scale(axis)*point(axis)**degree
            gradient(axis) = scale(axis)*degree*point(axis)**(degree - 1)
            if (degree >= 2) then
                hessian(axis, axis) = scale(axis)*degree*(degree - 1)* &
                    point(axis)**(degree - 2)
            end if
        end do
        if (degree < 2) return
        value = value + point(1)*point(2)**(degree - 1)
        gradient(1) = gradient(1) + point(2)**(degree - 1)
        gradient(2) = gradient(2) + &
            (degree - 1)*point(1)*point(2)**(degree - 2)
        hessian(1, 2) = (degree - 1)*point(2)**(degree - 2)
        hessian(2, 1) = hessian(1, 2)
        if (degree >= 3) then
            hessian(2, 2) = hessian(2, 2) + &
                (degree - 1)*(degree - 2)*point(1)*point(2)**(degree - 3)
        end if
        if (dimension /= 3) return
        value = value + point(1)*point(3)**(degree - 1) + &
            point(2)*point(3)**(degree - 1)
        gradient(1) = gradient(1) + point(3)**(degree - 1)
        gradient(2) = gradient(2) + point(3)**(degree - 1)
        gradient(3) = gradient(3) + &
            (degree - 1)*(point(1) + point(2))*point(3)**(degree - 2)
        hessian(1, 3) = (degree - 1)*point(3)**(degree - 2)
        hessian(3, 1) = hessian(1, 3)
        hessian(2, 3) = (degree - 1)*point(3)**(degree - 2)
        hessian(3, 2) = hessian(2, 3)
        if (degree >= 3) then
            hessian(3, 3) = hessian(3, 3) + (degree - 1)*(degree - 2)* &
                (point(1) + point(2))*point(3)**(degree - 3)
        end if
    end subroutine polynomial

    subroutine close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        largest_scaled_error = max(largest_scaled_error, &
            abs(actual - expected)/max(1.0_dp, abs(expected)))
        call check_condition(abs(actual - expected) <= &
            1.0e-10_dp*max(1.0_dp, abs(expected)), label)
        if (abs(actual - expected) <= 1.0e-10_dp*max(1.0_dp, abs(expected))) return
        print *, label, actual, expected
        error stop 'polynomial oracle'
    end subroutine close
end program test_lagrange_product_jet_polynomials
