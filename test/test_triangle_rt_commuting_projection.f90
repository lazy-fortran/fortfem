program test_triangle_rt_commuting_projection
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_triangle_rt_interpolant, initialize_triangle_raviart_thomas, &
        interpolate_triangle_rt, triangle_duffy_quadrature, &
        triangle_rt_basis_t
    use fortfem_kinds, only: dp
    implicit none

    type(triangle_rt_basis_t) :: basis
    real(dp), allocatable :: dofs(:), eta(:), weights(:), xi(:)
    real(dp) :: determinant, divergence_error_moment, exact_divergence
    real(dp) :: interpolated_divergence, interpolated_value(2)
    real(dp) :: jacobian(2, 2), physical_point(2), vertices(2, 3)
    integer :: degree, node, status, total_degree, x_degree, y_degree
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.2_dp, -0.1_dp]
    vertices(:, 2) = [1.7_dp, 0.3_dp]
    vertices(:, 3) = [-0.2_dp, 1.4_dp]
    jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
    jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
    determinant = jacobian(1, 1) * jacobian(2, 2) - &
        jacobian(1, 2) * jacobian(2, 1)
    call triangle_duffy_quadrature(12, xi, eta, weights, status)

    do degree = 0, 3
        call initialize_triangle_raviart_thomas(degree, basis, status)
        call interpolate_triangle_rt( &
            vertices, degree, 12, manufactured_field, dofs, status)
        call record_condition(status == 0 .and. &
            size(dofs) == (degree + 1) * (degree + 3), &
            "Physical RT interpolation has the exact trimmed dimension")
        do total_degree = 0, degree
            do x_degree = 0, total_degree
                y_degree = total_degree - x_degree
                divergence_error_moment = 0.0_dp
                do node = 1, size(weights)
                    physical_point = vertices(:, 1) + &
                        matmul(jacobian, [xi(node), eta(node)])
                    exact_divergence = manufactured_divergence( &
                        physical_point(1), physical_point(2))
                    call evaluate_triangle_rt_interpolant( &
                        vertices, basis, dofs, xi(node), eta(node), &
                        interpolated_value, interpolated_divergence, status)
                    divergence_error_moment = divergence_error_moment + &
                        determinant * weights(node) * &
                        (interpolated_divergence - exact_divergence) * &
                        xi(node)**x_degree * eta(node)**y_degree
                end do
                call record_condition(abs(divergence_error_moment) < &
                    3.0e-9_dp, &
                    "RT divergence projection is orthogonal to physical P_k")
            end do
        end do
        deallocate(dofs)
    end do

    call interpolate_triangle_rt( &
        vertices, -1, 2, manufactured_field, dofs, status)
    call record_condition(status /= 0, &
        "Physical RT interpolation rejects negative degree")

    call check_summary("Triangle RT commuting projection")
    if (.not. all_passed) error stop 1

contains

    pure subroutine manufactured_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value(1) = x**5 + x * y**3
        value(2) = y**5 + x**3 * y
    end subroutine manufactured_field

    pure function manufactured_divergence(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp) :: value

        value = 5.0_dp * x**4 + y**3 + 5.0_dp * y**4 + x**3
    end function manufactured_divergence

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_rt_commuting_projection
