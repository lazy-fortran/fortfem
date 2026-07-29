program test_triangle_lagrange_arbitrary_order
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_triangle_lagrange_basis, &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_lagrange_basis_t) :: basis, copied_basis
    real(dp), allocatable :: gradients(:, :), nodes(:, :), values(:)
    real(dp), allocatable :: nodal_coefficients(:)
    real(dp) :: expected_gradient(2), expected_value, point(2)
    integer :: degree, dof, status, x_degree, y_degree
    logical :: all_passed

    all_passed = .true.
    point = [0.23_dp, 0.31_dp]
    do degree = 0, 4
        call initialize_triangle_lagrange_basis(degree, basis, status)
        call record_condition(status == 0, &
            "Arbitrary-order triangle Lagrange basis initializes")
        call triangle_lagrange_nodes(basis, nodes)
        allocate(values(size(nodes, 2)), gradients(2, size(nodes, 2)))
        allocate(nodal_coefficients(size(nodes, 2)))

        do dof = 1, size(nodes, 2)
            call evaluate_triangle_lagrange_basis( &
                basis, nodes(1, dof), nodes(2, dof), values, gradients, status)
            call record_condition(status == 0 .and. &
                maxval(abs(values - unit_vector(size(values), dof))) < &
                2.0e-12_dp, &
                "Generated triangle basis has Kronecker nodal values")
        end do

        call evaluate_triangle_lagrange_basis( &
            basis, point(1), point(2), values, gradients, status)
        call record_condition(abs(sum(values) - 1.0_dp) < 2.0e-13_dp .and. &
            maxval(abs(sum(gradients, dim=2))) < 5.0e-13_dp, &
            "Triangle basis satisfies partition of unity and zero gradient sum")

        do x_degree = 0, degree
            do y_degree = 0, degree - x_degree
                nodal_coefficients = nodes(1, :)**x_degree * &
                    nodes(2, :)**y_degree
                expected_value = point(1)**x_degree * point(2)**y_degree
                expected_gradient = monomial_gradient( &
                    point, x_degree, y_degree)
                call record_condition(abs(dot_product( &
                    nodal_coefficients, values) - expected_value) < &
                    2.0e-12_dp, &
                    "Triangle interpolation reproduces a polynomial value")
                call record_condition(maxval(abs(matmul( &
                    gradients, nodal_coefficients) - expected_gradient)) < &
                    8.0e-12_dp, &
                    "Triangle interpolation reproduces a polynomial gradient")
            end do
        end do

        copied_basis = basis
        call initialize_triangle_lagrange_basis(0, basis, status)
        call evaluate_triangle_lagrange_basis( &
            copied_basis, point(1), point(2), values, gradients, status)
        call record_condition(status == 0 .and. &
            abs(sum(values) - 1.0_dp) < 2.0e-13_dp, &
            "Assigned triangle basis retains an independent deep copy")

        deallocate(values, gradients, nodal_coefficients, nodes)
    end do

    call initialize_triangle_lagrange_basis(-1, basis, status)
    call record_condition(status /= 0, &
        "Triangle Lagrange basis rejects a negative degree")

    call check_summary("Arbitrary-order triangle Lagrange basis")
    if (.not. all_passed) error stop 1

contains

    pure function unit_vector(length, active) result(vector)
        integer, intent(in) :: length, active
        real(dp) :: vector(length)

        vector = 0.0_dp
        vector(active) = 1.0_dp
    end function unit_vector

    pure function monomial_gradient(point, x_degree, y_degree) result(gradient)
        real(dp), intent(in) :: point(2)
        integer, intent(in) :: x_degree, y_degree
        real(dp) :: gradient(2)

        gradient = 0.0_dp
        if (x_degree > 0) then
            gradient(1) = real(x_degree, dp) * point(1)**(x_degree - 1) * &
                point(2)**y_degree
        end if
        if (y_degree > 0) then
            gradient(2) = real(y_degree, dp) * point(1)**x_degree * &
                point(2)**(y_degree - 1)
        end if
    end function monomial_gradient

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_lagrange_arbitrary_order
