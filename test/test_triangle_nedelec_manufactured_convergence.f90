program test_triangle_nedelec_manufactured_convergence
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_triangle_nedelec_interpolant, &
        initialize_triangle_nedelec_first_kind, interpolate_triangle_nedelec, &
        triangle_duffy_quadrature, triangle_nedelec_first_kind_t
    use fortfem_core, only: mesh_t, unit_square_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: coarse_curl_error, coarse_l2_error
    real(dp) :: curl_rate, fine_curl_error, fine_l2_error, l2_rate
    integer :: order
    logical :: all_passed

    all_passed = .true.
    do order = 1, 4
        call interpolation_errors(3, order, coarse_l2_error, coarse_curl_error)
        call interpolation_errors(5, order, fine_l2_error, fine_curl_error)
        l2_rate = log(coarse_l2_error / fine_l2_error) / log(2.0_dp)
        curl_rate = log(coarse_curl_error / fine_curl_error) / log(2.0_dp)
        call record_condition(fine_l2_error < coarse_l2_error .and. &
            fine_curl_error < coarse_curl_error, &
            "Manufactured H(curl) interpolation improves under refinement")
        call record_condition(l2_rate > real(order, dp) - 0.35_dp, &
            "Manufactured vector-field error reaches the expected h-order")
        call record_condition(curl_rate > real(order, dp) - 0.45_dp, &
            "Manufactured scalar-curl error reaches the expected h-order")
    end do

    call check_summary("Arbitrary-order Nedelec manufactured convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine interpolation_errors( &
            vertex_count, order, l2_error, curl_error)
        integer, intent(in) :: vertex_count, order
        real(dp), intent(out) :: l2_error, curl_error

        type(mesh_t) :: mesh
        type(triangle_nedelec_first_kind_t) :: basis
        real(dp), allocatable :: dofs(:), eta(:), weights(:), xi(:)
        real(dp) :: determinant, exact_curl, exact_value(2)
        real(dp) :: interpolated_curl, interpolated_value(2)
        real(dp) :: jacobian(2, 2), physical_point(2), vertices(2, 3)
        integer :: point, quadrature_degree, status, triangle

        mesh = unit_square_mesh(vertex_count)
        quadrature_degree = 2 * order + 6
        call initialize_triangle_nedelec_first_kind(order, basis, status)
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, status)
        l2_error = 0.0_dp
        curl_error = 0.0_dp
        do triangle = 1, mesh%data%n_triangles
            vertices = mesh%data%vertices(:, mesh%data%triangles(:, triangle))
            jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
            jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
            determinant = jacobian(1, 1) * jacobian(2, 2) - &
                jacobian(1, 2) * jacobian(2, 1)
            call interpolate_triangle_nedelec( &
                vertices, order, quadrature_degree, manufactured_field, &
                dofs, status)
            do point = 1, size(weights)
                physical_point = vertices(:, 1) + &
                    matmul(jacobian, [xi(point), eta(point)])
                call manufactured_field( &
                    physical_point(1), physical_point(2), exact_value)
                exact_curl = manufactured_curl( &
                    physical_point(1), physical_point(2))
                call evaluate_triangle_nedelec_interpolant( &
                    vertices, basis, dofs, xi(point), eta(point), &
                    interpolated_value, interpolated_curl, status)
                l2_error = l2_error + determinant * weights(point) * &
                    sum((interpolated_value - exact_value)**2)
                curl_error = curl_error + determinant * weights(point) * &
                    (interpolated_curl - exact_curl)**2
            end do
            deallocate(dofs)
        end do
        l2_error = sqrt(l2_error)
        curl_error = sqrt(curl_error)
        call mesh%destroy()
    end subroutine interpolation_errors

    pure subroutine manufactured_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        real(dp), parameter :: pi = acos(-1.0_dp)

        value(1) = 0.0_dp
        value(2) = sin(pi * x) * sin(pi * y)
    end subroutine manufactured_field

    pure function manufactured_curl(x, y) result(value)
        real(dp), intent(in) :: x, y
        real(dp) :: value

        real(dp), parameter :: pi = acos(-1.0_dp)

        value = pi * cos(pi * x) * sin(pi * y)
    end function manufactured_curl

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_triangle_nedelec_manufactured_convergence
