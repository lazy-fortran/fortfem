program test_sheet_current_surface_parity
    !! Geometry-independent fitted/resolved sheet-current oracle.
    !!
    !! The three geometry builders below deliberately provide only physical
    !! surface normals and quadrature measures to FortFEM.  The direct sums in
    !! this test are an independent oracle for a cylinder, sphere, and torus;
    !! no geometry-specific implementation is hidden in the library routine.
    use check, only: check_condition, check_summary
    use fortfem_sheet_current_surface_parity, only: &
        compare_sheet_current_surface_representations_direct => &
        compare_sheet_current_surface_representations
    use fortfem_interop, only: &
        compare_sheet_current_surface_representations, &
        compare_sheet_current_surface_representations_jvp, &
        compare_sheet_current_surface_representations_interop => &
        compare_sheet_current_surface_representations, &
        compare_sheet_current_surface_representations_jvp_interop => &
        compare_sheet_current_surface_representations_jvp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: thickness = 0.035_dp
    real(dp), parameter :: half_width = 5.0_dp*thickness
    integer, parameter :: normal_count = 401
    integer, parameter :: sphere_theta = 5, sphere_phi = 8
    integer, parameter :: cylinder_phi = 8, cylinder_z = 5
    integer, parameter :: torus_theta = 5, torus_phi = 8
    real(dp), parameter :: axis_a(3) = [0.35_dp, -0.45_dp, 0.8_dp]
    real(dp), parameter :: axis_b(3) = [-0.55_dp, 0.2_dp, 0.65_dp]
    real(dp), parameter :: eps = 2.0e-7_dp
    real(dp) :: distance(normal_count), normal_weights(normal_count)
    real(dp) :: gaussian, spacing, thickness_dot
    real(dp) :: fitted(3), resolved(3), relative_error
    real(dp) :: fitted_dot(3), resolved_dot(3), relative_error_dot
    real(dp) :: fitted_plus(3), fitted_minus(3), resolved_plus(3), resolved_minus(3)
    real(dp) :: relative_plus, relative_minus
    real(dp), allocatable :: normals(:, :), plus_field(:, :), minus_field(:, :)
    real(dp), allocatable :: normals_dot(:, :), plus_dot(:, :), minus_dot(:, :)
    real(dp), allocatable :: surface_weights(:), surface_weights_dot(:)
    real(dp), allocatable :: surface_weights_plus(:), surface_weights_minus(:)
    real(dp), allocatable :: distance_dot(:), normal_weights_dot(:)
    type(fortsparse_status_t) :: status
    integer :: sample

    abstract interface
        subroutine geometry_builder(normals, weights, plus_field, minus_field, &
                normals_dot, weights_dot, plus_dot, minus_dot, with_direction)
            import dp
            real(dp), allocatable, intent(out) :: normals(:, :), plus_field(:, :)
            real(dp), allocatable, intent(out) :: minus_field(:, :)
            real(dp), allocatable, intent(out) :: normals_dot(:, :)
            real(dp), allocatable, intent(out) :: weights(:), weights_dot(:)
            real(dp), allocatable, intent(out) :: plus_dot(:, :)
            real(dp), allocatable, intent(out) :: minus_dot(:, :)
            logical, intent(in) :: with_direction
        end subroutine geometry_builder
    end interface

    spacing = 2.0_dp*half_width/real(normal_count - 1, dp)
    do sample = 1, normal_count
        distance(sample) = -half_width + spacing*real(sample - 1, dp)
        normal_weights(sample) = spacing
    end do
    normal_weights([1, normal_count]) = 0.5_dp*spacing

    call check_geometry("slab", 4, build_slab)
    call check_geometry("sphere", sphere_theta*sphere_phi, build_sphere)
    call check_geometry("cylinder", cylinder_phi*cylinder_z, build_cylinder)
    call check_geometry("torus", torus_theta*torus_phi, build_torus, .true.)
    call check_summary("sheet-current surface parity")

contains

    subroutine check_geometry(name, point_count, builder, derivative_case)
        character(len=*), intent(in) :: name
        integer, intent(in) :: point_count
        procedure(geometry_builder) :: builder
        logical, intent(in), optional :: derivative_case
        logical :: with_direction
        integer :: point, sample
        real(dp) :: reference(3), independent_resolved(3), dot_current(3)
        real(dp) :: fitted_direct(3), resolved_direct(3), relative_direct
        real(dp) :: fitted_interop(3), resolved_interop(3), relative_interop
        real(dp) :: fitted_dot_interop(3), resolved_dot_interop(3)
        real(dp) :: relative_dot_interop
        real(dp) :: reversed(3), oriented(3), invalid_fitted(3), invalid_resolved(3)
        real(dp) :: invalid_error
        real(dp), allocatable :: invalid_weights(:), invalid_normals(:, :)
        real(dp) :: normal_dot_profile(3), direct_norm

        with_direction = .false.
        if (present(derivative_case)) with_direction = derivative_case
        call builder(normals, surface_weights, plus_field, minus_field, &
            normals_dot, surface_weights_dot, plus_dot, minus_dot, with_direction)
        allocate(surface_weights_plus(point_count), surface_weights_minus(point_count))
        surface_weights_plus = surface_weights + eps*surface_weights_dot
        surface_weights_minus = surface_weights - eps*surface_weights_dot

        call compare_sheet_current_surface_representations( &
            plus_field, minus_field, normals, surface_weights, distance, &
            normal_weights, thickness, fitted, resolved, relative_error, status)
        call check_condition(status%code == 0, trim(name)// &
            " surface representation accepts fixed-topology geometry")
        call compare_sheet_current_surface_representations_direct( &
            plus_field, minus_field, normals, surface_weights, distance, &
            normal_weights, thickness, fitted_direct, resolved_direct, &
            relative_direct, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(fitted_direct - fitted)) < 2.0e-13_dp .and. &
            maxval(abs(resolved_direct - resolved)) < 2.0e-13_dp .and. &
            abs(relative_direct - relative_error) < 2.0e-13_dp, trim(name)// &
            " defining-module and umbrella canonical names agree")
        call compare_sheet_current_surface_representations_interop( &
            plus_field, minus_field, normals, surface_weights, distance, &
            normal_weights, thickness, fitted_interop, resolved_interop, &
            relative_interop, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(fitted_interop - fitted)) < 2.0e-13_dp .and. &
            maxval(abs(resolved_interop - resolved)) < 2.0e-13_dp .and. &
            abs(relative_interop - relative_error) < 2.0e-13_dp, &
            trim(name)//" interop facade preserves surface ledger")

        reference = 0.0_dp
        independent_resolved = 0.0_dp
        do point = 1, point_count
            call cross_product(normals(point, :), &
                plus_field(point, :) - minus_field(point, :), dot_current)
            reference = reference + surface_weights(point)*dot_current
            do sample = 1, normal_count
                gaussian = exp(-(distance(sample)/thickness)**2)/ &
                    (sqrt(pi)*thickness)
                independent_resolved = independent_resolved + &
                    surface_weights(point)*normal_weights(sample)*gaussian*dot_current
            end do
        end do
        call check_condition(maxval(abs(fitted - reference)) < 2.0e-13_dp, &
            trim(name)//" fitted integrated current has independent ledger")
        call check_condition(maxval(abs(resolved - independent_resolved)) < 2.0e-13_dp, &
            trim(name)//" regularized current has independent Gaussian ledger")
        direct_norm = max(1.0_dp, sqrt(sum(reference**2)))
        call check_condition(relative_error < 2.0e-12_dp .and. &
            sqrt(sum((resolved - fitted)**2))/direct_norm < 2.0e-12_dp, &
            trim(name)//" preserves integrated current")

        ! Reversing only the supplied normal reverses the Ampere current.
        call compare_sheet_current_surface_representations( &
            plus_field, minus_field, -normals, surface_weights, distance, &
            normal_weights, thickness, reversed, invalid_resolved, invalid_error, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(reversed + fitted)) < 2.0e-13_dp, trim(name)// &
            " reverses integrated current with the normal")

        ! Reversing both orientation and trace ordering leaves K invariant.
        call compare_sheet_current_surface_representations( &
            minus_field, plus_field, -normals, surface_weights, distance, &
            normal_weights, thickness, oriented, invalid_resolved, invalid_error, status)
        call check_condition(status%code == 0 .and. &
            maxval(abs(oriented - fitted)) < 2.0e-13_dp, trim(name)// &
            " preserves current under a full interface orientation reversal")

        ! A changed row count is a topology error, while a zero measure and a
        ! non-unit normal are invalid quadrature data.  The contract must reject
        ! all three before writing a misleading ledger.
        allocate(invalid_weights(size(surface_weights) - 1), &
            invalid_normals(size(normals, 1), size(normals, 2)))
        invalid_weights = surface_weights(1:size(surface_weights) - 1)
        invalid_normals = normals
        call compare_sheet_current_surface_representations( &
            plus_field, minus_field, normals, invalid_weights, distance, normal_weights, &
            thickness, invalid_fitted, invalid_resolved, invalid_error, status)
        call check_condition(status%code /= 0, trim(name)// &
            " rejects a changed surface quadrature topology")
        invalid_weights = surface_weights
        invalid_weights(1) = 0.0_dp
        call compare_sheet_current_surface_representations( &
            plus_field, minus_field, normals, invalid_weights, distance, normal_weights, &
            thickness, invalid_fitted, invalid_resolved, invalid_error, status)
        call check_condition(status%code /= 0, trim(name)// &
            " rejects a non-positive surface measure")
        invalid_normals = normals
        invalid_normals(1, :) = 1.1_dp*invalid_normals(1, :)
        call compare_sheet_current_surface_representations( &
            plus_field, minus_field, invalid_normals, surface_weights, distance, &
            normal_weights, thickness, invalid_fitted, invalid_resolved, invalid_error, status)
        call check_condition(status%code /= 0, trim(name)// &
            " rejects a non-unit interface normal")
        deallocate(invalid_weights, invalid_normals)

        if (with_direction) then
            allocate(distance_dot(normal_count), normal_weights_dot(normal_count))
            distance_dot = 0.01_dp*distance
            normal_weights_dot = 0.015_dp*normal_weights
            thickness_dot = 0.02_dp
            call compare_sheet_current_surface_representations_jvp( &
                plus_field, minus_field, normals, surface_weights, distance, &
                normal_weights, thickness, plus_dot, minus_dot, normals_dot, &
                surface_weights_dot, distance_dot, normal_weights_dot, thickness_dot, &
                fitted_dot, resolved_dot, relative_error_dot, status)
            call compare_sheet_current_surface_representations_jvp_interop( &
                plus_field, minus_field, normals, surface_weights, distance, &
                normal_weights, thickness, plus_dot, minus_dot, normals_dot, &
                surface_weights_dot, distance_dot, normal_weights_dot, thickness_dot, &
                fitted_dot_interop, resolved_dot_interop, relative_dot_interop, status)
            call check_condition(status%code == 0 .and. &
                maxval(abs(fitted_dot_interop - fitted_dot)) < 2.0e-13_dp .and. &
                maxval(abs(resolved_dot_interop - resolved_dot)) < 2.0e-13_dp .and. &
                abs(relative_dot_interop - relative_error_dot) < 2.0e-13_dp, &
                trim(name)//" interop facade preserves surface JVP")
            call compare_sheet_current_surface_representations( &
                plus_field + eps*plus_dot, minus_field + eps*minus_dot, &
                normals + eps*normals_dot, surface_weights_plus, distance + eps*distance_dot, &
                normal_weights + eps*normal_weights_dot, thickness + eps*thickness_dot, &
                fitted_plus, resolved_plus, relative_plus, status)
            call compare_sheet_current_surface_representations( &
                plus_field - eps*plus_dot, minus_field - eps*minus_dot, &
                normals - eps*normals_dot, surface_weights_minus, distance - eps*distance_dot, &
                normal_weights - eps*normal_weights_dot, thickness - eps*thickness_dot, &
                fitted_minus, resolved_minus, relative_minus, status)
            call check_condition(maxval(abs((fitted_plus - fitted_minus)/(2.0_dp*eps) - &
                fitted_dot)) < 2.0e-8_dp .and. &
                maxval(abs((resolved_plus - resolved_minus)/(2.0_dp*eps) - &
                resolved_dot)) < 2.0e-8_dp, trim(name)// &
                " fixed-topology geometry/current JVP matches central differences")
            deallocate(distance_dot, normal_weights_dot)
        end if

        deallocate(normals, plus_field, minus_field, normals_dot, plus_dot, minus_dot)
        deallocate(surface_weights, surface_weights_dot)
        deallocate(surface_weights_plus, surface_weights_minus)
    end subroutine check_geometry

    subroutine build_slab(slab_normals, weights, slab_plus, slab_minus, &
            normals_direction, weights_direction, plus_direction, minus_direction, &
            with_direction)
        real(dp), allocatable, intent(out) :: slab_normals(:, :), slab_plus(:, :)
        real(dp), allocatable, intent(out) :: slab_minus(:, :), normals_direction(:, :)
        real(dp), allocatable, intent(out) :: weights(:), weights_direction(:)
        real(dp), allocatable, intent(out) :: plus_direction(:, :), minus_direction(:, :)
        logical, intent(in) :: with_direction
        integer :: point

        allocate(slab_normals(4, 3), slab_plus(4, 3), slab_minus(4, 3), &
            normals_direction(4, 3), weights(4), weights_direction(4), &
            plus_direction(4, 3), minus_direction(4, 3))
        slab_normals = 0.0_dp
        slab_normals(:, 1) = 1.0_dp
        weights = 0.25_dp
        slab_minus = 0.0_dp
        slab_plus = 0.0_dp
        slab_plus(:, 2) = 1.0_dp
        normals_direction = 0.0_dp
        weights_direction = 0.0_dp
        plus_direction = 0.0_dp
        minus_direction = 0.0_dp
        if (with_direction) then
            do point = 1, 4
                weights_direction(point) = 0.01_dp*real(point, dp)
                plus_direction(point, 3) = 0.005_dp*real(point, dp)
            end do
        end if
    end subroutine build_slab

    subroutine build_sphere(sphere_normals, weights, sphere_plus, sphere_minus, &
            normals_direction, weights_direction, plus_direction, minus_direction, &
            with_direction)
        real(dp), allocatable, intent(out) :: sphere_normals(:, :), sphere_plus(:, :)
        real(dp), allocatable, intent(out) :: sphere_minus(:, :), normals_direction(:, :)
        real(dp), allocatable, intent(out) :: weights(:), weights_direction(:)
        real(dp), allocatable, intent(out) :: plus_direction(:, :), minus_direction(:, :)
        logical, intent(in) :: with_direction
        integer :: theta_index, phi_index, point
        real(dp) :: theta, phi, dtheta, dphi, normal(3), current(3)
        allocate(sphere_normals(sphere_theta*sphere_phi, 3), &
            sphere_plus(sphere_theta*sphere_phi, 3), &
            sphere_minus(sphere_theta*sphere_phi, 3), &
            normals_direction(sphere_theta*sphere_phi, 3), &
            weights(sphere_theta*sphere_phi), weights_direction(sphere_theta*sphere_phi), &
            plus_direction(sphere_theta*sphere_phi, 3), &
            minus_direction(sphere_theta*sphere_phi, 3))
        dtheta = pi/real(sphere_theta, dp)
        dphi = 2.0_dp*pi/real(sphere_phi, dp)
        point = 0
        do theta_index = 1, sphere_theta
            theta = (real(theta_index, dp) - 0.5_dp)*dtheta
            do phi_index = 1, sphere_phi
                phi = (real(phi_index, dp) - 0.5_dp)*dphi
                point = point + 1
                normal = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
                sphere_normals(point, :) = normal
                weights(point) = sin(theta)*dtheta*dphi
                call tangent_current(normal, current)
                sphere_plus(point, :) = cross(current, normal)
                sphere_minus(point, :) = 0.0_dp
                normals_direction(point, :) = cross(normal, [0.04_dp, -0.02_dp, 0.03_dp])
                weights_direction(point) = 0.02_dp*weights(point)*sin(phi)
                plus_direction(point, :) = [0.01_dp*sin(phi), -0.02_dp*cos(phi), &
                    0.015_dp*sin(theta)]
                minus_direction(point, :) = [0.005_dp*cos(phi), 0.003_dp*sin(phi), &
                    -0.004_dp*cos(theta)]
            end do
        end do
    end subroutine build_sphere

    subroutine build_cylinder(cylinder_normals, weights, cylinder_plus, cylinder_minus, &
            normals_direction, weights_direction, plus_direction, minus_direction, &
            with_direction)
        real(dp), allocatable, intent(out) :: cylinder_normals(:, :), cylinder_plus(:, :)
        real(dp), allocatable, intent(out) :: cylinder_minus(:, :), normals_direction(:, :)
        real(dp), allocatable, intent(out) :: weights(:), weights_direction(:)
        real(dp), allocatable, intent(out) :: plus_direction(:, :), minus_direction(:, :)
        logical, intent(in) :: with_direction
        integer :: z_index, phi_index, point
        real(dp) :: phi, z, dz, dphi, normal(3), current(3)
        allocate(cylinder_normals(cylinder_phi*cylinder_z, 3), &
            cylinder_plus(cylinder_phi*cylinder_z, 3), &
            cylinder_minus(cylinder_phi*cylinder_z, 3), &
            normals_direction(cylinder_phi*cylinder_z, 3), &
            weights(cylinder_phi*cylinder_z), weights_direction(cylinder_phi*cylinder_z), &
            plus_direction(cylinder_phi*cylinder_z, 3), &
            minus_direction(cylinder_phi*cylinder_z, 3))
        dz = 2.0_dp/real(cylinder_z, dp)
        dphi = 2.0_dp*pi/real(cylinder_phi, dp)
        point = 0
        do z_index = 1, cylinder_z
            z = -1.0_dp + (real(z_index, dp) - 0.5_dp)*dz
            do phi_index = 1, cylinder_phi
                phi = (real(phi_index, dp) - 0.5_dp)*dphi
                point = point + 1
                normal = [cos(phi), sin(phi), 0.0_dp]
                cylinder_normals(point, :) = normal
                weights(point) = dz*dphi
                call tangent_current(normal, current)
                cylinder_plus(point, :) = cross(current, normal)
                cylinder_minus(point, :) = 0.0_dp
                normals_direction(point, :) = cross(normal, [0.02_dp, 0.01_dp, -0.03_dp])
                weights_direction(point) = -0.017_dp*weights(point)*cos(phi)
                plus_direction(point, :) = [0.015_dp*z, 0.01_dp*sin(phi), -0.02_dp*cos(phi)]
                minus_direction(point, :) = [0.004_dp*sin(phi), -0.005_dp*cos(phi), &
                    0.003_dp*z]
            end do
        end do
    end subroutine build_cylinder

    subroutine build_torus(torus_normals, weights, torus_plus, torus_minus, &
            normals_direction, weights_direction, plus_direction, minus_direction, &
            with_direction)
        real(dp), allocatable, intent(out) :: torus_normals(:, :), torus_plus(:, :)
        real(dp), allocatable, intent(out) :: torus_minus(:, :), normals_direction(:, :)
        real(dp), allocatable, intent(out) :: weights(:), weights_direction(:)
        real(dp), allocatable, intent(out) :: plus_direction(:, :), minus_direction(:, :)
        logical, intent(in) :: with_direction
        integer :: theta_index, phi_index, point
        real(dp) :: theta, phi, dtheta, dphi, major, minor, normal(3), current(3)
        allocate(torus_normals(torus_theta*torus_phi, 3), &
            torus_plus(torus_theta*torus_phi, 3), &
            torus_minus(torus_theta*torus_phi, 3), &
            normals_direction(torus_theta*torus_phi, 3), &
            weights(torus_theta*torus_phi), weights_direction(torus_theta*torus_phi), &
            plus_direction(torus_theta*torus_phi, 3), &
            minus_direction(torus_theta*torus_phi, 3))
        major = 2.4_dp
        minor = 0.65_dp
        dtheta = 2.0_dp*pi/real(torus_theta, dp)
        dphi = 2.0_dp*pi/real(torus_phi, dp)
        point = 0
        do theta_index = 1, torus_theta
            theta = (real(theta_index, dp) - 0.5_dp)*dtheta
            do phi_index = 1, torus_phi
                phi = (real(phi_index, dp) - 0.5_dp)*dphi
                point = point + 1
                normal = [cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)]
                torus_normals(point, :) = normal
                weights(point) = minor*(major + minor*cos(theta))*dtheta*dphi
                call tangent_current(normal, current)
                torus_plus(point, :) = cross(current, normal)
                torus_minus(point, :) = 0.0_dp
                normals_direction(point, :) = cross(normal, [0.03_dp, 0.02_dp, -0.01_dp])
                weights_direction(point) = 0.012_dp*weights(point)*sin(theta + phi)
                plus_direction(point, :) = [0.013_dp*sin(theta), -0.008_dp*cos(phi), &
                    0.011_dp*sin(phi)]
                minus_direction(point, :) = [0.003_dp*cos(theta), 0.002_dp*sin(phi), &
                    -0.006_dp*cos(phi)]
            end do
        end do
    end subroutine build_torus

    subroutine tangent_current(normal, current)
        real(dp), intent(in) :: normal(3)
        real(dp), intent(out) :: current(3)
        call cross_product(normal, axis_a, current)
        current = current*dot_product(normal, axis_b)
    end subroutine tangent_current

    pure function cross(first, second) result(result)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: result(3)
        result = [first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross

    pure subroutine cross_product(first, second, result)
        real(dp), intent(in) :: first(3), second(3)
        real(dp), intent(out) :: result(3)
        result = cross(first, second)
    end subroutine cross_product

end program test_sheet_current_surface_parity
