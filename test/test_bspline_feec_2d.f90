program test_bspline_feec_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_bspline_derivative_matrix, build_bspline_feec_2d_operators, &
        build_bspline_feec_3d_operators, evaluate_bspline_basis, &
        evaluate_nurbs_surface_geometry, map_isogeometric_h1_gradient, &
        evaluate_nurbs_volume_geometry, map_isogeometric_hcurl, &
        map_isogeometric_hdiv, map_isogeometric_l2
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
    real(dp), allocatable :: volume_controls(:, :, :, :)
    real(dp), allocatable :: volume_weights(:, :, :)
    real(dp) :: jacobian(3, 2), mapped_point(3)
    real(dp), allocatable :: mapped_hcurl(:, :), mapped_hdiv(:, :)
    real(dp), allocatable :: mapped_gradient(:, :), mapped_l2(:)
    real(dp) :: affine_jacobian(2, 2), determinant
    real(dp) :: reference_vector(2, 1)
    real(dp) :: volume_jacobian(3, 3), volume_point(3)
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

    affine_jacobian = reshape([2.0_dp, 0.25_dp, 0.5_dp, 1.5_dp], [2, 2])
    reference_vector(:, 1) = [1.0_dp, -0.5_dp]
    call map_isogeometric_h1_gradient( &
        affine_jacobian, reference_vector, mapped_gradient, determinant, status)
    call map_isogeometric_hcurl( &
        affine_jacobian, reference_vector, mapped_hcurl, determinant, status)
    call map_isogeometric_hdiv( &
        affine_jacobian, reference_vector, mapped_hdiv, determinant, status)
    call map_isogeometric_l2( &
        affine_jacobian, [2.0_dp], mapped_l2, determinant, status)
    call check_condition(status == 0 .and. &
        abs(determinant - 2.875_dp) < 2.0e-14_dp .and. &
        maxval(abs(mapped_gradient - mapped_hcurl)) < 2.0e-14_dp, &
        "H1 gradients and Hcurl fields use the covariant Piola map")
    call check_condition(abs( &
        determinant*dot_product(mapped_hcurl(:, 1), mapped_hdiv(:, 1)) - &
        dot_product(reference_vector(:, 1), reference_vector(:, 1))) < &
        3.0e-14_dp, "Covariant and contravariant Piola maps preserve pairing")
    call check_condition(abs(mapped_l2(1) - 2.0_dp/2.875_dp) < 2.0e-14_dp, &
        "L2 density uses the determinant-preserving pullback")

    call make_affine_volume_control_net( &
        knots_x, knots_y, knots_y, 2, volume_controls, volume_weights)
    call evaluate_nurbs_volume_geometry( &
        knots_x, knots_y, knots_y, 2, 2, 2, volume_controls, volume_weights, &
        0.2_dp, 0.3_dp, 0.4_dp, volume_point, volume_jacobian, status)
    call check_condition(status == 0 .and. maxval(abs(volume_point - &
        [1.4_dp, -0.1_dp, 2.1_dp])) < 5.0e-14_dp, &
        "NURBS volume reproduces an affine physical map exactly")
    call check_condition(maxval(abs(volume_jacobian - reshape([ &
        2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 4.0_dp], [3, 3]))) < 4.0e-13_dp, &
        "NURBS volume returns the analytical three-dimensional Jacobian")

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

    subroutine make_affine_volume_control_net( &
            knots_x_local, knots_y_local, knots_z_local, degree, points, weights)
        real(dp), intent(in) :: knots_x_local(:), knots_y_local(:)
        real(dp), intent(in) :: knots_z_local(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: points(:, :, :, :)
        real(dp), allocatable, intent(out) :: weights(:, :, :)

        real(dp), allocatable :: x_points(:), y_points(:), z_points(:)
        integer :: ix, iy, iz

        x_points = greville_abscissae(knots_x_local, degree)
        y_points = greville_abscissae(knots_y_local, degree)
        z_points = greville_abscissae(knots_z_local, degree)
        allocate(points(3, size(x_points), size(y_points), size(z_points)))
        allocate(weights(size(x_points), size(y_points), size(z_points)))
        weights = 1.0_dp
        do iz = 1, size(z_points)
            do iy = 1, size(y_points)
                do ix = 1, size(x_points)
                    points(:, ix, iy, iz) = [ &
                        1.0_dp + 2.0_dp*x_points(ix), &
                        -1.0_dp + 3.0_dp*y_points(iy), &
                        0.5_dp + 4.0_dp*z_points(iz)]
                end do
            end do
        end do
    end subroutine make_affine_volume_control_net

end program test_bspline_feec_2d
