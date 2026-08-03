program test_physical_surface_geometry
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        sample_physical_surface_geometry, &
        sample_physical_surface_geometry_jvp, &
        sample_physical_surface_geometry_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: tolerance = 3.0e-11_dp
    real(dp), parameter :: finite_difference_tolerance = 4.0e-9_dp
    real(dp), parameter :: step = 1.0e-6_dp
    real(dp) :: reference_coordinates(2, 3), map_points(3, 3)
    real(dp) :: map_tangents(3, 2, 3), physical_points(3, 3)
    real(dp) :: normals(3, 3), surface_jacobian(3)
    real(dp) :: map_points_dot(3, 3), map_tangents_dot(3, 2, 3)
    real(dp) :: physical_points_dot(3, 3), normals_dot(3, 3)
    real(dp) :: surface_jacobian_dot(3)
    real(dp) :: physical_points_plus(3, 3), physical_points_minus(3, 3)
    real(dp) :: normals_plus(3, 3), normals_minus(3, 3)
    real(dp) :: surface_jacobian_plus(3), surface_jacobian_minus(3)
    real(dp) :: physical_points_bar(3, 3), normals_bar(3, 3)
    real(dp) :: surface_jacobian_bar(3), reference_coordinates_bar(2, 3)
    real(dp) :: map_points_bar(3, 3), map_tangents_bar(3, 2, 3)
    real(dp) :: expected_cross(3), expected_normal(3), expected_jacobian
    real(dp) :: forward_pairing, reverse_pairing
    real(dp) :: degenerate_tangents(3, 2, 3)
    integer :: quadrature
    logical :: all_passed
    type(fortsparse_status_t) :: status

    all_passed = .true.
    reference_coordinates = reshape([ &
        0.1_dp, 0.2_dp, 0.3_dp, &
        0.4_dp, 0.5_dp, 0.6_dp], [2, 3])
    map_points = reshape([ &
        1.0_dp, 2.0_dp, 0.2_dp, &
        1.3_dp, 2.2_dp, 0.4_dp, &
        1.7_dp, 2.5_dp, 0.8_dp], [3, 3])
    map_tangents(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.2_dp, &
        0.1_dp, 1.0_dp, 0.3_dp], [3, 2])
    map_tangents(:, :, 2) = reshape([0.8_dp, 0.2_dp, 0.1_dp, &
        -0.2_dp, 1.1_dp, 0.4_dp], [3, 2])
    map_tangents(:, :, 3) = reshape([1.2_dp, -0.1_dp, 0.3_dp, &
        0.2_dp, 0.9_dp, -0.2_dp], [3, 2])

    call sample_physical_surface_geometry( &
        reference_coordinates, map_points, map_tangents, 1, physical_points, &
        normals, surface_jacobian, status)
    call record_condition(status%code == 0, &
        "surface sampler accepts finite nondegenerate map samples")
    call record_condition(maxval(abs(physical_points - map_points)) < tolerance, &
        "surface sampler returns the physical map coordinates")
    do quadrature = 1, 3
        expected_cross = cross3( &
            map_tangents(:, 1, quadrature), map_tangents(:, 2, quadrature))
        expected_jacobian = sqrt(dot_product(expected_cross, expected_cross))
        expected_normal = expected_cross/expected_jacobian
        call record_condition(abs(surface_jacobian(quadrature) - &
            expected_jacobian) < tolerance, &
            "surface sampler returns the positive cross-product measure")
        call record_condition(maxval(abs(normals(:, quadrature) - &
            expected_normal)) < tolerance, &
            "surface sampler returns an oriented unit normal")
    end do

    map_points_dot = reshape([ &
        0.03_dp, -0.02_dp, 0.04_dp, &
        -0.01_dp, 0.05_dp, 0.02_dp, &
        0.06_dp, -0.03_dp, 0.01_dp], [3, 3])
    map_tangents_dot = reshape([ &
        0.01_dp, -0.03_dp, 0.02_dp, 0.04_dp, 0.02_dp, -0.01_dp, &
        -0.02_dp, 0.03_dp, 0.01_dp, 0.02_dp, -0.04_dp, 0.03_dp, &
        0.03_dp, 0.01_dp, -0.02_dp, -0.01_dp, 0.04_dp, 0.02_dp], [3, 2, 3])
    call sample_physical_surface_geometry_jvp( &
        reference_coordinates, map_points, map_tangents, 1, map_points_dot, &
        map_tangents_dot, physical_points_dot, normals_dot, &
        surface_jacobian_dot, status)
    call record_condition(status%code == 0, &
        "surface sampler fixed-topology JVP succeeds")
    call sample_physical_surface_geometry( &
        reference_coordinates, map_points + step*map_points_dot, &
        map_tangents + step*map_tangents_dot, 1, physical_points_plus, &
        normals_plus, surface_jacobian_plus, status)
    call sample_physical_surface_geometry( &
        reference_coordinates, map_points - step*map_points_dot, &
        map_tangents - step*map_tangents_dot, 1, physical_points_minus, &
        normals_minus, surface_jacobian_minus, status)
    call record_condition(maxval(abs(physical_points_dot - &
        (physical_points_plus - physical_points_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "surface sampler point JVP matches central differences")
    call record_condition(maxval(abs(normals_dot - &
        (normals_plus - normals_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "surface sampler normal JVP matches central differences")
    call record_condition(maxval(abs(surface_jacobian_dot - &
        (surface_jacobian_plus - surface_jacobian_minus)/(2.0_dp*step))) < &
        finite_difference_tolerance, &
        "surface sampler metric JVP matches central differences")

    physical_points_bar = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, 0.1_dp, &
        0.6_dp, 0.2_dp, -0.3_dp], [3, 3])
    normals_bar = reshape([ &
        -0.3_dp, 0.2_dp, 0.4_dp, 0.1_dp, -0.5_dp, 0.2_dp, &
        0.3_dp, 0.6_dp, -0.2_dp], [3, 3])
    surface_jacobian_bar = [0.7_dp, -0.2_dp, 0.4_dp]
    call sample_physical_surface_geometry_vjp( &
        reference_coordinates, map_points, map_tangents, 1, &
        physical_points_bar, normals_bar, surface_jacobian_bar, &
        reference_coordinates_bar, map_points_bar, map_tangents_bar, status)
    call record_condition(status%code == 0, &
        "surface sampler fixed-topology VJP succeeds")
    forward_pairing = sum(physical_points_bar*physical_points_dot) + &
        sum(normals_bar*normals_dot) + &
        sum(surface_jacobian_bar*surface_jacobian_dot)
    reverse_pairing = sum(reference_coordinates_bar*0.0_dp) + &
        sum(map_points_bar*map_points_dot) + &
        sum(map_tangents_bar*map_tangents_dot)
    call record_condition(abs(forward_pairing - reverse_pairing) < tolerance, &
        "surface sampler VJP satisfies the real dot-product identity")

    degenerate_tangents = map_tangents
    degenerate_tangents(:, 2, 2) = degenerate_tangents(:, 1, 2)
    call sample_physical_surface_geometry( &
        reference_coordinates, map_points, degenerate_tangents, 1, &
        physical_points, normals, surface_jacobian, status)
    call record_condition(status%code /= 0, &
        "surface sampler rejects degenerate tangent pairs")
    call sample_physical_surface_geometry( &
        reference_coordinates, map_points, map_tangents, 0, physical_points, &
        normals, surface_jacobian, status)
    call record_condition(status%code /= 0, &
        "surface sampler rejects an invalid orientation sign")

    call check_summary("physical surface geometry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

    pure function cross3(first, second) result(cross_product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: cross_product(3)

        cross_product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross3

end program test_physical_surface_geometry
