module fortfem_regularized_surface_current_layer
    !! Differentiable Gaussian approximation of K delta_Gamma.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_generated_regularized_surface_current, only: &
        generated_regularized_surface_current
    use fortfem_generated_regularized_surface_current_jvp, only: &
        generated_regularized_surface_current_jvp
    use fortfem_generated_regularized_surface_current_vjp, only: &
        generated_regularized_surface_current_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    real(dp), parameter :: inverse_sqrt_pi = 1.0_dp/sqrt(acos(-1.0_dp))

    public :: evaluate_regularized_surface_current_layer
    public :: evaluate_regularized_surface_current_layer_jvp
    public :: evaluate_regularized_surface_current_layer_vjp
    public :: evaluate_regularized_surface_current_integral

contains

    pure subroutine evaluate_regularized_surface_current_layer( &
            signed_distance, sheet_current, thickness, volume_current, status)
        real(dp), intent(in) :: signed_distance(:), sheet_current(:, :), thickness
        real(dp), intent(out) :: volume_current(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, component

        volume_current = 0.0_dp
        call validate_primal( &
            signed_distance, sheet_current, thickness, volume_current, status)
        if (status%code /= FORTSPARSE_OK) return

        do sample = 1, size(signed_distance)
            do component = 1, size(sheet_current, 2)
                call generated_regularized_surface_current( &
                    signed_distance(sample), sheet_current(sample, component), &
                    thickness, inverse_sqrt_pi, &
                    volume_current(sample, component))
            end do
        end do
        if (.not. all(ieee_is_finite(volume_current))) then
            volume_current = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current layer produced non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_regularized_surface_current_layer

    pure subroutine evaluate_regularized_surface_current_layer_jvp( &
            signed_distance, sheet_current, thickness, signed_distance_dot, &
            sheet_current_dot, thickness_dot, volume_current_dot, status)
        real(dp), intent(in) :: signed_distance(:), sheet_current(:, :), thickness
        real(dp), intent(in) :: signed_distance_dot(:), sheet_current_dot(:, :)
        real(dp), intent(in) :: thickness_dot
        real(dp), intent(out) :: volume_current_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, component

        volume_current_dot = 0.0_dp
        call validate_primal( &
            signed_distance, sheet_current, thickness, volume_current_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(signed_distance_dot) /= size(signed_distance)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current JVP has incompatible distances")
            return
        end if
        if (any(shape(sheet_current_dot) /= shape(sheet_current))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current JVP has incompatible currents")
            return
        end if
        if (.not. all(ieee_is_finite(signed_distance_dot)) .or. &
            .not. all(ieee_is_finite(sheet_current_dot)) .or. &
            .not. ieee_is_finite(thickness_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current JVP received non-finite directions")
            return
        end if

        do sample = 1, size(signed_distance)
            do component = 1, size(sheet_current, 2)
                call generated_regularized_surface_current_jvp( &
                    signed_distance(sample), sheet_current(sample, component), &
                    thickness, inverse_sqrt_pi, signed_distance_dot(sample), &
                    sheet_current_dot(sample, component), thickness_dot, &
                    volume_current_dot(sample, component))
            end do
        end do
        if (.not. all(ieee_is_finite(volume_current_dot))) then
            volume_current_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current JVP produced non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_regularized_surface_current_layer_jvp

    pure subroutine evaluate_regularized_surface_current_layer_vjp( &
            signed_distance, sheet_current, thickness, volume_current_bar, &
            signed_distance_bar, sheet_current_bar, thickness_bar, status)
        real(dp), intent(in) :: signed_distance(:), sheet_current(:, :), thickness
        real(dp), intent(in) :: volume_current_bar(:, :)
        real(dp), intent(out) :: signed_distance_bar(:), sheet_current_bar(:, :)
        real(dp), intent(out) :: thickness_bar
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, component
        real(dp) :: distance_contribution, thickness_contribution

        signed_distance_bar = 0.0_dp
        sheet_current_bar = 0.0_dp
        thickness_bar = 0.0_dp
        call validate_primal( &
            signed_distance, sheet_current, thickness, volume_current_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(signed_distance_bar) /= size(signed_distance) .or. &
            any(shape(sheet_current_bar) /= shape(sheet_current))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current VJP has incompatible cotangents")
            return
        end if
        if (.not. all(ieee_is_finite(volume_current_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current VJP received non-finite cotangents")
            return
        end if

        do sample = 1, size(signed_distance)
            do component = 1, size(sheet_current, 2)
                call generated_regularized_surface_current_vjp( &
                    signed_distance(sample), sheet_current(sample, component), &
                    thickness, inverse_sqrt_pi, &
                    volume_current_bar(sample, component), distance_contribution, &
                    sheet_current_bar(sample, component), thickness_contribution)
                signed_distance_bar(sample) = signed_distance_bar(sample) + &
                    distance_contribution
                thickness_bar = thickness_bar + thickness_contribution
            end do
        end do
        if (.not. all(ieee_is_finite(signed_distance_bar)) .or. &
            .not. all(ieee_is_finite(sheet_current_bar)) .or. &
            .not. ieee_is_finite(thickness_bar)) then
            signed_distance_bar = 0.0_dp
            sheet_current_bar = 0.0_dp
            thickness_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current VJP produced non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_regularized_surface_current_layer_vjp

    pure subroutine evaluate_regularized_surface_current_integral( &
            signed_distance, normal_weights, sheet_current, thickness, &
            normalization, integrated_current, status)
        real(dp), intent(in) :: signed_distance(:), normal_weights(:)
        real(dp), intent(in) :: sheet_current(:, :), thickness
        real(dp), intent(out) :: normalization, integrated_current(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, component
        real(dp) :: profile

        normalization = 0.0_dp
        integrated_current = 0.0_dp
        call validate_primal( &
            signed_distance, sheet_current, thickness, sheet_current, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(normal_weights) /= size(signed_distance) .or. &
            size(integrated_current) /= size(sheet_current, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current diagnostic has incompatible arrays")
            return
        end if
        if (.not. all(ieee_is_finite(normal_weights)) .or. &
            any(normal_weights <= 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current diagnostic needs positive weights")
            return
        end if

        do sample = 1, size(signed_distance)
            call generated_regularized_surface_current( &
                signed_distance(sample), 1.0_dp, thickness, inverse_sqrt_pi, profile)
            normalization = normalization + normal_weights(sample)*profile
            do component = 1, size(sheet_current, 2)
                integrated_current(component) = integrated_current(component) + &
                    normal_weights(sample)*profile*sheet_current(sample, component)
            end do
        end do
        if (.not. ieee_is_finite(normalization) .or. &
            .not. all(ieee_is_finite(integrated_current))) then
            normalization = 0.0_dp
            integrated_current = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "regularized surface-current diagnostic produced non-finite values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_regularized_surface_current_integral

    pure subroutine validate_primal( &
            signed_distance, sheet_current, thickness, volume_current, status)
        real(dp), intent(in) :: signed_distance(:), sheet_current(:, :), thickness
        real(dp), intent(in) :: volume_current(:, :)
        type(fortsparse_status_t), intent(out) :: status

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "regularized surface-current layer received incompatible arrays")
        if (size(signed_distance) < 1 .or. size(sheet_current, 2) < 1) return
        if (size(sheet_current, 1) /= size(signed_distance)) return
        if (any(shape(volume_current) /= shape(sheet_current))) return
        if (.not. ieee_is_finite(thickness)) return
        if (thickness <= 0.0_dp) return
        if (.not. all(ieee_is_finite(signed_distance))) return
        if (.not. all(ieee_is_finite(sheet_current))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_primal

end module fortfem_regularized_surface_current_layer
