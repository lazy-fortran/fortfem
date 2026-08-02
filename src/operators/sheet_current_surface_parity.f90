module fortfem_sheet_current_surface_parity
    !! Geometry-neutral fitted/resolved surface-current parity contract.
    !!
    !! A caller owns a fixed surface quadrature: the normal and positive
    !! surface measure may come from a fitted mesh, a cut surface, DG/HDG
    !! skeleton, or an IGA patch.  The contract builds the Ampere trace
    !! current from plus/minus fields, integrates it on that surface, and
    !! compares it with the same current represented by a normalized Gaussian
    !! layer in a caller-owned signed-distance quadrature.  No cylinder,
    !! sphere, or torus geometry is selected here; those geometries provide
    !! the physical normals and measures to the same ledger.
    !!
    !! The JVP holds the quadrature topology and orientation fixed.  Surface
    !! normals and measures may vary smoothly, but a quadrature-row or
    !! orientation change requires a new contract evaluation.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_regularized_surface_current_layer, only: &
        evaluate_regularized_surface_current_layer, &
        evaluate_regularized_surface_current_layer_jvp
    use fortfem_surface_current, only: &
        assemble_interface_surface_current, &
        assemble_interface_surface_current_jvp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: compare_sheet_current_surface_representations
    public :: compare_sheet_current_surface_representations_jvp

contains

    subroutine compare_sheet_current_surface_representations( &
            plus_field, minus_field, normals, surface_weights, signed_distance, &
            normal_weights, thickness, fitted_integrated, regularized_integrated, &
            relative_error, status)
        !! Compare a fitted Ampere trace with its resolved Gaussian layer.
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        real(dp), intent(in) :: signed_distance(:), normal_weights(:), thickness
        real(dp), intent(out) :: fitted_integrated(:), regularized_integrated(:)
        real(dp), intent(out) :: relative_error
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: surface_current(:, :)
        real(dp), allocatable :: profile(:, :)
        integer :: surface_point

        fitted_integrated = 0.0_dp
        regularized_integrated = 0.0_dp
        relative_error = 0.0_dp
        call validate_inputs( &
            plus_field, minus_field, normals, surface_weights, signed_distance, &
            normal_weights, thickness, fitted_integrated, regularized_integrated, status)
        if (status%code /= FORTSPARSE_OK) return

        allocate(surface_current(size(surface_weights), 3), &
            profile(size(signed_distance), 3))
        call assemble_interface_surface_current( &
            plus_field, minus_field, normals, surface_weights, surface_current, &
            fitted_integrated, status)
        if (status%code /= FORTSPARSE_OK) return

        do surface_point = 1, size(surface_weights)
            call evaluate_regularized_surface_current_layer( &
                signed_distance, &
                spread(surface_current(surface_point, :), 1, size(signed_distance)), &
                thickness, profile, status)
            if (status%code /= FORTSPARSE_OK) return
            regularized_integrated = regularized_integrated + &
                surface_weights(surface_point)*sum_weighted_profile( &
                normal_weights, profile)
        end do
        call evaluate_relative_error( &
            fitted_integrated, regularized_integrated, relative_error, status)
    end subroutine compare_sheet_current_surface_representations

    subroutine compare_sheet_current_surface_representations_jvp( &
            plus_field, minus_field, normals, surface_weights, signed_distance, &
            normal_weights, thickness, plus_dot, minus_dot, normals_dot, &
            surface_weights_dot, signed_distance_dot, normal_weights_dot, thickness_dot, &
            fitted_integrated_dot, regularized_integrated_dot, relative_error_dot, status)
        !! Fixed-topology JVP of both integrated current representations.
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :)
        real(dp), intent(in) :: normals(:, :), surface_weights(:)
        real(dp), intent(in) :: signed_distance(:), normal_weights(:), thickness
        real(dp), intent(in) :: plus_dot(:, :), minus_dot(:, :), normals_dot(:, :)
        real(dp), intent(in) :: surface_weights_dot(:), signed_distance_dot(:)
        real(dp), intent(in) :: normal_weights_dot(:), thickness_dot
        real(dp), intent(out) :: fitted_integrated_dot(:)
        real(dp), intent(out) :: regularized_integrated_dot(:)
        real(dp), intent(out) :: relative_error_dot
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: surface_current(:, :), surface_current_dot(:, :)
        real(dp), allocatable :: profile(:, :), profile_dot(:, :)
        real(dp), allocatable :: profile_current(:, :), profile_current_dot(:, :)
        real(dp) :: fitted_integrated(3), regularized_integrated(3)
        real(dp) :: relative_error
        integer :: surface_point

        fitted_integrated_dot = 0.0_dp
        regularized_integrated_dot = 0.0_dp
        relative_error_dot = 0.0_dp
        call validate_inputs( &
            plus_field, minus_field, normals, surface_weights, signed_distance, &
            normal_weights, thickness, fitted_integrated_dot, &
            regularized_integrated_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(plus_dot, 1) /= size(surface_weights) .or. &
            size(plus_dot, 2) /= 3 .or. size(minus_dot, 1) /= size(surface_weights) .or. &
            size(minus_dot, 2) /= 3 .or. size(normals_dot, 1) /= size(surface_weights) .or. &
            size(normals_dot, 2) /= 3 .or. size(surface_weights_dot) /= size(surface_weights) .or. &
            size(signed_distance_dot) /= size(signed_distance) .or. &
            size(normal_weights_dot) /= size(normal_weights)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current parity JVP received incompatible increments")
            return
        end if
        if (.not. all(ieee_is_finite(plus_dot)) .or. &
            .not. all(ieee_is_finite(minus_dot)) .or. &
            .not. all(ieee_is_finite(normals_dot)) .or. &
            .not. all(ieee_is_finite(surface_weights_dot)) .or. &
            .not. all(ieee_is_finite(signed_distance_dot)) .or. &
            .not. all(ieee_is_finite(normal_weights_dot)) .or. &
            .not. ieee_is_finite(thickness_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current parity JVP received non-finite increments")
            return
        end if

        allocate(surface_current(size(surface_weights), 3), &
            surface_current_dot(size(surface_weights), 3), &
            profile_current(size(signed_distance), 3), &
            profile_current_dot(size(signed_distance), 3), &
            profile(size(signed_distance), 3), profile_dot(size(signed_distance), 3))
        call assemble_interface_surface_current( &
            plus_field, minus_field, normals, surface_weights, surface_current, &
            fitted_integrated, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_interface_surface_current_jvp( &
            plus_field, minus_field, normals, surface_weights, plus_dot, minus_dot, &
            normals_dot, surface_weights_dot, surface_current_dot, &
            fitted_integrated_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        do surface_point = 1, size(surface_weights)
            profile_current = spread(surface_current(surface_point, :), 1, &
                size(signed_distance))
            profile_current_dot = spread(surface_current_dot(surface_point, :), 1, &
                size(signed_distance))
            call evaluate_regularized_surface_current_layer( &
                signed_distance, profile_current, thickness, profile, status)
            if (status%code /= FORTSPARSE_OK) return
            call evaluate_regularized_surface_current_layer_jvp( &
                signed_distance, profile_current, thickness, signed_distance_dot, &
                profile_current_dot, thickness_dot, profile_dot, status)
            if (status%code /= FORTSPARSE_OK) return
            regularized_integrated_dot = regularized_integrated_dot + &
                surface_weights_dot(surface_point)*sum_weighted_profile( &
                normal_weights, profile) + &
                surface_weights(surface_point)*sum_weighted_profile_dot( &
                normal_weights, normal_weights_dot, profile, profile_dot)
        end do

        call compare_sheet_current_surface_representations( &
            plus_field, minus_field, normals, surface_weights, signed_distance, &
            normal_weights, thickness, fitted_integrated, regularized_integrated, &
            relative_error, status)
        if (status%code /= FORTSPARSE_OK) return
        call evaluate_relative_error_jvp( &
            fitted_integrated, regularized_integrated, fitted_integrated_dot, &
            regularized_integrated_dot, relative_error_dot, status)
    end subroutine compare_sheet_current_surface_representations_jvp

    pure function sum_weighted_profile(weights, profile) result(result)
        real(dp), intent(in) :: weights(:), profile(:, :)
        real(dp) :: result(size(profile, 2))
        integer :: sample, component

        result = 0.0_dp
        do sample = 1, size(weights)
            do component = 1, size(profile, 2)
                result(component) = result(component) + &
                    weights(sample)*profile(sample, component)
            end do
        end do
    end function sum_weighted_profile

    pure function sum_weighted_profile_dot( &
            weights, weights_dot, profile, profile_dot) result(result)
        real(dp), intent(in) :: weights(:), weights_dot(:)
        real(dp), intent(in) :: profile(:, :), profile_dot(:, :)
        real(dp) :: result(size(profile, 2))
        integer :: sample, component

        result = 0.0_dp
        do sample = 1, size(weights)
            do component = 1, size(profile, 2)
                result(component) = result(component) + &
                    weights_dot(sample)*profile(sample, component) + &
                    weights(sample)*profile_dot(sample, component)
            end do
        end do
    end function sum_weighted_profile_dot

    subroutine validate_inputs( &
            plus_field, minus_field, normals, surface_weights, signed_distance, &
            normal_weights, thickness, fitted_integrated, regularized_integrated, status)
        real(dp), intent(in) :: plus_field(:, :), minus_field(:, :), normals(:, :)
        real(dp), intent(in) :: surface_weights(:), signed_distance(:), normal_weights(:)
        real(dp), intent(in) :: thickness
        real(dp), intent(in) :: fitted_integrated(:), regularized_integrated(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: surface_point

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "surface-current parity received incompatible arrays")
        if (size(plus_field, 1) < 1 .or. size(plus_field, 2) /= 3 .or. &
            any(shape(minus_field) /= shape(plus_field)) .or. &
            any(shape(normals) /= shape(plus_field)) .or. &
            size(surface_weights) /= size(plus_field, 1) .or. &
            size(signed_distance) < 2 .or. size(normal_weights) /= size(signed_distance) .or. &
            size(fitted_integrated) /= 3 .or. size(regularized_integrated) /= 3) return
        if (.not. ieee_is_finite(thickness) .or. thickness <= 0.0_dp) return
        if (.not. all(ieee_is_finite(plus_field)) .or. &
            .not. all(ieee_is_finite(minus_field)) .or. &
            .not. all(ieee_is_finite(normals)) .or. &
            .not. all(ieee_is_finite(surface_weights)) .or. &
            .not. all(ieee_is_finite(signed_distance)) .or. &
            .not. all(ieee_is_finite(normal_weights)) .or. &
            any(surface_weights <= 0.0_dp) .or. any(normal_weights <= 0.0_dp)) return
        do surface_point = 1, size(surface_weights)
            if (abs(dot_product(normals(surface_point, :), &
                normals(surface_point, :)) - 1.0_dp) > unit_tolerance) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_inputs

    subroutine evaluate_relative_error( &
            fitted_integrated, regularized_integrated, relative_error, status)
        real(dp), intent(in) :: fitted_integrated(:), regularized_integrated(:)
        real(dp), intent(out) :: relative_error
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: scale

        scale = max(1.0_dp, sqrt(sum(fitted_integrated**2)))
        relative_error = sqrt(sum((regularized_integrated - fitted_integrated)**2))/scale
        if (.not. ieee_is_finite(relative_error)) then
            relative_error = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current parity produced a non-finite error")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_relative_error

    subroutine evaluate_relative_error_jvp( &
            fitted_integrated, regularized_integrated, fitted_integrated_dot, &
            regularized_integrated_dot, relative_error_dot, status)
        real(dp), intent(in) :: fitted_integrated(:), regularized_integrated(:)
        real(dp), intent(in) :: fitted_integrated_dot(:), regularized_integrated_dot(:)
        real(dp), intent(out) :: relative_error_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp) :: difference(3), difference_dot(3), scale, scale_dot, norm_difference

        difference = regularized_integrated - fitted_integrated
        difference_dot = regularized_integrated_dot - fitted_integrated_dot
        norm_difference = sqrt(sum(difference**2))
        scale = max(1.0_dp, sqrt(sum(fitted_integrated**2)))
        scale_dot = 0.0_dp
        if (scale > 1.0_dp) scale_dot = dot_product(fitted_integrated, &
            fitted_integrated_dot)/scale
        relative_error_dot = 0.0_dp
        if (norm_difference > 0.0_dp) relative_error_dot = &
            dot_product(difference, difference_dot)/(norm_difference*scale) - &
            norm_difference*scale_dot/scale**2
        if (.not. ieee_is_finite(relative_error_dot)) then
            relative_error_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "surface-current parity JVP produced a non-finite error")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_relative_error_jvp

end module fortfem_sheet_current_surface_parity
