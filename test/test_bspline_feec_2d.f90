program test_bspline_feec_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_bspline_derivative_matrix, build_bspline_feec_2d_operators, &
        build_bspline_feec_3d_operators, evaluate_bspline_basis, &
        evaluate_nurbs_surface_geometry
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: knots_x(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_y(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), allocatable :: basis(:), curl_matrix(:, :)
    real(dp), allocatable :: derivative(:), derivative_matrix(:, :)
    real(dp), allocatable :: gradient_matrix(:, :), coefficients(:)
    real(dp), allocatable :: curl_3d(:, :), divergence_3d(:, :)
    real(dp), allocatable :: gradient_3d(:, :)
    real(dp), allocatable :: control_points(:, :, :), nurbs_weights(:, :)
    real(dp) :: jacobian(3, 2), mapped_point(3)
    real(dp) :: points(5), reproduced, x
    integer :: i, status

    points = [0.0_dp, 0.13_dp, 0.35_dp, 0.82_dp, 1.0_dp]
    do i = 1, size(points)
        x = points(i)
        call evaluate_bspline_basis( &
            knots_x, 2, x, basis, derivative, status)
        call check_condition(status == 0 .and. &
            abs(sum(basis) - 1.0_dp) < 3.0e-14_dp, &
            "B-splines form a partition of unity")
        allocate(coefficients(size(basis)))
        coefficients = greville_abscissae(knots_x, 2)
        reproduced = dot_product(coefficients, basis)
        call check_condition(abs(reproduced - x) < 3.0e-14_dp, &
            "Quadratic B-splines reproduce the coordinate exactly")
        call check_condition(abs(dot_product(coefficients, derivative) - &
            1.0_dp) < 2.0e-13_dp, &
            "Analytical B-spline derivatives reproduce a unit gradient")
        deallocate(basis, derivative, coefficients)
    end do

    call build_bspline_derivative_matrix( &
        knots_x, 2, derivative_matrix, status)
    coefficients = greville_abscissae(knots_x, 2)
    call check_condition(status == 0 .and. maxval(abs( &
        matmul(derivative_matrix, coefficients) - 1.0_dp)) < 3.0e-14_dp, &
        "Spline coefficient derivative differentiates the coordinate exactly")

    call build_bspline_feec_2d_operators( &
        knots_x, knots_y, 2, 2, gradient_matrix, curl_matrix, status)
    call check_condition(status == 0 .and. &
        all(shape(gradient_matrix) == [31, 20]) .and. &
        all(shape(curl_matrix) == [12, 31]), &
        "Tensor-product spline de Rham operators have exact dimensions")
    call check_condition(maxval(abs(matmul( &
        curl_matrix, gradient_matrix))) < 3.0e-14_dp, &
        "Tensor-product spline complex satisfies curl grad equals zero")
    call check_condition(maxval(abs(matmul(gradient_matrix, &
        [(1.0_dp, i=1, size(gradient_matrix, 2))]))) < 3.0e-14_dp, &
        "Isogeometric gradient annihilates constants")

    call build_bspline_feec_3d_operators( &
        knots_x, knots_y, knots_y, 2, 2, 2, gradient_3d, curl_3d, &
        divergence_3d, status)
    call check_condition(status == 0 .and. &
        all(shape(gradient_3d) == [184, 80]) .and. &
        all(shape(curl_3d) == [141, 184]) .and. &
        all(shape(divergence_3d) == [36, 141]), &
        "3D spline de Rham operators have exact FEEC dimensions")
    call check_condition(maxval(abs(matmul(curl_3d, gradient_3d))) < &
        4.0e-14_dp, "3D isogeometric complex satisfies curl grad equals zero")
    call check_condition(maxval(abs(matmul(divergence_3d, curl_3d))) < &
        4.0e-14_dp, "3D isogeometric complex satisfies div curl equals zero")

    call make_affine_control_net( &
        knots_x, knots_y, 2, 2, control_points, nurbs_weights)
    call evaluate_nurbs_surface_geometry( &
        knots_x, knots_y, 2, 2, control_points, nurbs_weights, &
        0.23_dp, 0.61_dp, mapped_point, jacobian, status)
    call check_condition(status == 0 .and. maxval(abs(mapped_point - &
        [2.69_dp, 0.22_dp, 0.12_dp])) < 4.0e-14_dp, &
        "NURBS geometry reproduces an affine physical surface exactly")
    call check_condition(maxval(abs(jacobian - reshape([ &
        3.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 2.0_dp, -1.0_dp], [3, 2]))) < &
        3.0e-13_dp, "NURBS geometry returns the analytical Jacobian")

    call check_summary("Structure-preserving B-spline FEEC")

contains

    pure function greville_abscissae(knots, degree) result(points)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable :: points(:)
        integer :: basis

        allocate(points(size(knots) - degree - 1))
        do basis = 1, size(points)
            points(basis) = sum(knots(basis + 1:basis + degree))/real(degree, dp)
        end do
    end function greville_abscissae

    subroutine make_affine_control_net( &
            knots_x_local, knots_y_local, degree_x, degree_y, points, weights)
        real(dp), intent(in) :: knots_x_local(:), knots_y_local(:)
        integer, intent(in) :: degree_x, degree_y
        real(dp), allocatable, intent(out) :: points(:, :, :), weights(:, :)

        real(dp), allocatable :: x_points(:), y_points(:)
        integer :: ix, iy

        x_points = greville_abscissae(knots_x_local, degree_x)
        y_points = greville_abscissae(knots_y_local, degree_y)
        allocate(points(3, size(x_points), size(y_points)))
        allocate(weights(size(x_points), size(y_points)))
        weights = 1.0_dp
        do iy = 1, size(y_points)
            do ix = 1, size(x_points)
                points(:, ix, iy) = [ &
                    2.0_dp + 3.0_dp*x_points(ix), &
                    -1.0_dp + 2.0_dp*y_points(iy), &
                    0.5_dp + x_points(ix) - y_points(iy)]
            end do
        end do
    end subroutine make_affine_control_net

end program test_bspline_feec_2d
