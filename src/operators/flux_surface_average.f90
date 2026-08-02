module fortfem_flux_surface_average
    !! Fixed-topology weighted reduction on a sampled surface.
    !!
    !! For samples s(c,q), quadrature weights w(q), and optional positive
    !! geometric factors g(q), this module evaluates
    !!
    !!   D       = sum_q w(q) g(q),
    !!   average = sum_q w(q) g(q) s(:,q) / D.
    !!
    !! The component-major sample layout is shared by scalar (one component)
    !! and vector-valued diagnostics.  Geometry and topology are caller-owned;
    !! the optional factor is useful for a surface Jacobian, volume factor, or
    !! any other positive measure already evaluated by a geometry layer.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_flux_surface_average
    public :: evaluate_flux_surface_average_jvp
    public :: evaluate_flux_surface_average_vjp

contains

    subroutine evaluate_flux_surface_average( &
            samples, weights, average, denominator, status, measure_factors)
        real(dp), intent(in) :: samples(:, :), weights(:)
        real(dp), intent(out) :: average(:), denominator
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: measure_factors(:)
        real(dp), allocatable :: effective_weights(:)
        integer :: component, sample

        average = 0.0_dp
        denominator = 0.0_dp
        if (size(average) /= size(samples, 1)) then
            call invalid_status(status, "flux-surface average output shape mismatch")
            return
        end if
        if (present(measure_factors)) then
            call validate_base_inputs(samples, weights, measure_factors, status)
        else
            call validate_base_inputs_without_factors(samples, weights, status)
        end if
        if (status%code /= FORTSPARSE_OK) return

        allocate(effective_weights(size(weights)))
        if (present(measure_factors)) then
            effective_weights = weights*measure_factors
        else
            effective_weights = weights
        end if
        if (.not. all(ieee_is_finite(effective_weights))) then
            call invalid_status(status, "flux-surface effective weights are non-finite")
            deallocate(effective_weights)
            return
        end if
        denominator = sum(effective_weights)
        if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_dp) then
            call invalid_status(status, "flux-surface denominator is invalid")
            denominator = 0.0_dp
            deallocate(effective_weights)
            return
        end if
        do component = 1, size(samples, 1)
            do sample = 1, size(samples, 2)
                average(component) = average(component) + effective_weights(sample)* &
                    samples(component, sample)
            end do
            average(component) = average(component)/denominator
        end do
        if (.not. all(ieee_is_finite(average))) then
            call invalid_status(status, "flux-surface average is non-finite")
            average = 0.0_dp
            denominator = 0.0_dp
            deallocate(effective_weights)
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
        deallocate(effective_weights)
    end subroutine evaluate_flux_surface_average

    subroutine evaluate_flux_surface_average_jvp( &
            samples, weights, samples_dot, weights_dot, average_dot, &
            denominator_dot, status, measure_factors, measure_factors_dot)
        real(dp), intent(in) :: samples(:, :), weights(:)
        real(dp), intent(in) :: samples_dot(:, :), weights_dot(:)
        real(dp), intent(out) :: average_dot(:), denominator_dot
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: measure_factors(:), measure_factors_dot(:)
        real(dp), allocatable :: effective_weights(:), effective_weights_dot(:)
        real(dp), allocatable :: average(:)
        real(dp) :: denominator, numerator_dot
        integer :: component, sample

        average_dot = 0.0_dp
        denominator_dot = 0.0_dp
        if (size(average_dot) /= size(samples, 1)) then
            call invalid_status(status, "flux-surface JVP output shape mismatch")
            return
        end if
        if (present(measure_factors_dot)) then
            if (.not. present(measure_factors)) then
                call invalid_status(status, "flux-surface factor direction has no base")
                return
            end if
            if (size(measure_factors_dot) /= size(weights)) then
                call invalid_status(status, &
                    "flux-surface factor direction shape mismatch")
                return
            end if
            if (.not. all(ieee_is_finite(measure_factors_dot))) then
                call invalid_status(status, &
                    "flux-surface factor direction is non-finite")
                return
            end if
        end if
        if (present(measure_factors)) then
            call validate_base_inputs(samples, weights, measure_factors, status)
        else
            call validate_base_inputs_without_factors(samples, weights, status)
        end if
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_direction_shapes( &
            samples, weights, samples_dot, weights_dot)) then
            call invalid_status(status, "flux-surface JVP direction shape mismatch")
            return
        end if
        if (.not. all(ieee_is_finite(samples_dot)) .or. &
            .not. all(ieee_is_finite(weights_dot))) then
            call invalid_status(status, "flux-surface JVP direction is non-finite")
            return
        end if

        allocate(effective_weights(size(weights)))
        allocate(effective_weights_dot(size(weights)))
        allocate(average(size(samples, 1)))
        if (present(measure_factors)) then
            effective_weights = weights*measure_factors
            effective_weights_dot = weights_dot*measure_factors
            if (present(measure_factors_dot)) then
                effective_weights_dot = effective_weights_dot + &
                    weights*measure_factors_dot
            end if
        else
            effective_weights = weights
            effective_weights_dot = weights_dot
        end if
        if (.not. all(ieee_is_finite(effective_weights)) .or. &
            .not. all(ieee_is_finite(effective_weights_dot))) then
            call invalid_status(status, "flux-surface JVP weights are non-finite")
            deallocate(effective_weights, effective_weights_dot, average)
            return
        end if
        denominator_dot = sum(effective_weights_dot)
        call evaluate_average_from_effective_weights( &
            samples, effective_weights, average, denominator, status)
        if (status%code /= FORTSPARSE_OK) then
            denominator_dot = 0.0_dp
            deallocate(effective_weights, effective_weights_dot, average)
            return
        end if
        do component = 1, size(samples, 1)
            numerator_dot = 0.0_dp
            do sample = 1, size(samples, 2)
                numerator_dot = numerator_dot + effective_weights(sample)* &
                    samples_dot(component, sample) + effective_weights_dot(sample)* &
                    samples(component, sample)
            end do
            average_dot(component) = (numerator_dot - average(component)* &
                denominator_dot)/denominator
        end do
        if (.not. ieee_is_finite(denominator_dot) .or. &
            .not. all(ieee_is_finite(average_dot))) then
            average_dot = 0.0_dp
            denominator_dot = 0.0_dp
            call invalid_status(status, "flux-surface JVP is non-finite")
            deallocate(effective_weights, effective_weights_dot, average)
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
        deallocate(effective_weights, effective_weights_dot, average)
    end subroutine evaluate_flux_surface_average_jvp

    subroutine evaluate_flux_surface_average_vjp( &
            samples, weights, average_bar, denominator_bar, samples_bar, &
            weights_bar, status, measure_factors, measure_factors_bar)
        real(dp), intent(in) :: samples(:, :), weights(:), average_bar(:)
        real(dp), intent(in) :: denominator_bar
        real(dp), intent(out) :: samples_bar(:, :), weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: measure_factors(:)
        real(dp), intent(out), optional :: measure_factors_bar(:)
        real(dp), allocatable :: effective_weights(:), average(:)
        real(dp) :: denominator, effective_weight_bar
        integer :: component, sample

        samples_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (present(measure_factors_bar)) measure_factors_bar = 0.0_dp
        if (present(measure_factors_bar)) then
            if (.not. present(measure_factors)) then
                call invalid_status(status, "flux-surface factor cotangent has no base")
                return
            end if
            if (size(measure_factors_bar) /= size(weights)) then
                call invalid_status(status, &
                    "flux-surface factor cotangent shape mismatch")
                return
            end if
        end if
        if (size(samples_bar, 1) /= size(samples, 1) .or. &
            size(samples_bar, 2) /= size(samples, 2) .or. &
            size(weights_bar) /= size(weights) .or. &
            size(average_bar) /= size(samples, 1)) then
            call invalid_status(status, "flux-surface VJP output shape mismatch")
            return
        end if
        if (.not. ieee_is_finite(denominator_bar) .or. &
            .not. all(ieee_is_finite(average_bar))) then
            call invalid_status(status, "flux-surface VJP cotangent is non-finite")
            return
        end if
        if (present(measure_factors)) then
            call validate_base_inputs(samples, weights, measure_factors, status)
        else
            call validate_base_inputs_without_factors(samples, weights, status)
        end if
        if (status%code /= FORTSPARSE_OK) return

        allocate(effective_weights(size(weights)))
        allocate(average(size(samples, 1)))
        if (present(measure_factors)) then
            effective_weights = weights*measure_factors
        else
            effective_weights = weights
        end if
        call evaluate_average_from_effective_weights( &
            samples, effective_weights, average, denominator, status)
        if (status%code /= FORTSPARSE_OK) then
            deallocate(effective_weights, average)
            return
        end if
        do sample = 1, size(weights)
            effective_weight_bar = denominator_bar
            do component = 1, size(samples, 1)
                samples_bar(component, sample) = average_bar(component)* &
                    effective_weights(sample)/denominator
                effective_weight_bar = effective_weight_bar + average_bar(component)* &
                    (samples(component, sample) - average(component))/denominator
            end do
            if (present(measure_factors)) then
                weights_bar(sample) = measure_factors(sample)*effective_weight_bar
            else
                weights_bar(sample) = effective_weight_bar
            end if
            if (present(measure_factors_bar)) then
                measure_factors_bar(sample) = weights(sample)*effective_weight_bar
            end if
        end do
        if (.not. all(ieee_is_finite(samples_bar)) .or. &
            .not. all(ieee_is_finite(weights_bar))) then
            samples_bar = 0.0_dp
            weights_bar = 0.0_dp
            if (present(measure_factors_bar)) measure_factors_bar = 0.0_dp
            call invalid_status(status, "flux-surface VJP is non-finite")
            deallocate(effective_weights, average)
            return
        end if
        if (present(measure_factors_bar)) then
            if (.not. all(ieee_is_finite(measure_factors_bar))) then
                samples_bar = 0.0_dp
                weights_bar = 0.0_dp
                measure_factors_bar = 0.0_dp
                call invalid_status(status, "flux-surface factor VJP is non-finite")
                deallocate(effective_weights, average)
                return
            end if
        end if
        call status_set(status, FORTSPARSE_OK, "")
        deallocate(effective_weights, average)
    end subroutine evaluate_flux_surface_average_vjp

    subroutine validate_base_inputs(samples, weights, measure_factors, status)
        real(dp), intent(in) :: samples(:, :), weights(:), measure_factors(:)
        type(fortsparse_status_t), intent(out) :: status

        if (size(samples, 1) <= 0 .or. size(samples, 2) <= 0 .or. &
            size(weights) <= 0 .or. size(samples, 2) /= size(weights)) then
            call invalid_status(status, &
                "flux-surface samples and weights are empty or incompatible")
            return
        end if
        if (size(measure_factors) /= size(weights)) then
            call invalid_status(status, "flux-surface geometric factors shape mismatch")
            return
        end if
        if (.not. all(ieee_is_finite(samples)) .or. &
            .not. all(ieee_is_finite(weights)) .or. &
            .not. all(ieee_is_finite(measure_factors)) .or. &
            any(weights <= 0.0_dp) .or. any(measure_factors <= 0.0_dp)) then
            call invalid_status(status, &
                "flux-surface inputs must be finite and positive in measure")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_base_inputs

    subroutine validate_base_inputs_without_factors(samples, weights, status)
        real(dp), intent(in) :: samples(:, :), weights(:)
        type(fortsparse_status_t), intent(out) :: status

        if (size(samples, 1) <= 0 .or. size(samples, 2) <= 0 .or. &
            size(weights) <= 0 .or. size(samples, 2) /= size(weights)) then
            call invalid_status(status, &
                "flux-surface samples and weights are empty or incompatible")
            return
        end if
        if (.not. all(ieee_is_finite(samples)) .or. &
            .not. all(ieee_is_finite(weights)) .or. any(weights <= 0.0_dp)) then
            call invalid_status(status, &
                "flux-surface inputs must be finite and positive in measure")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_base_inputs_without_factors

    subroutine evaluate_average_from_effective_weights( &
            samples, effective_weights, average, denominator, status)
        real(dp), intent(in) :: samples(:, :), effective_weights(:)
        real(dp), intent(out) :: average(:), denominator
        type(fortsparse_status_t), intent(out) :: status
        integer :: component, sample

        average = 0.0_dp
        if (size(average) /= size(samples, 1)) then
            call invalid_status(status, "flux-surface effective average shape mismatch")
            denominator = 0.0_dp
            return
        end if
        if (size(effective_weights) /= size(samples, 2)) then
            call invalid_status(status, "flux-surface effective weight shape mismatch")
            denominator = 0.0_dp
            return
        end if
        denominator = sum(effective_weights)
        if (denominator <= 0.0_dp) then
            call invalid_status(status, "flux-surface effective denominator is invalid")
            denominator = 0.0_dp
            return
        end if
        if (.not. ieee_is_finite(denominator)) then
            call invalid_status(status, "flux-surface effective denominator is invalid")
            denominator = 0.0_dp
            return
        end if
        do component = 1, size(samples, 1)
            do sample = 1, size(samples, 2)
                average(component) = average(component) + effective_weights(sample)* &
                    samples(component, sample)
            end do
            average(component) = average(component)/denominator
        end do
        if (.not. all(ieee_is_finite(average))) then
            call invalid_status(status, "flux-surface effective average is non-finite")
            average = 0.0_dp
            denominator = 0.0_dp
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_average_from_effective_weights

    logical function valid_direction_shapes( &
            samples, weights, samples_dot, weights_dot) result(valid)
        real(dp), intent(in) :: samples(:, :), weights(:)
        real(dp), intent(in) :: samples_dot(:, :), weights_dot(:)

        valid = size(samples_dot, 1) == size(samples, 1) .and. &
            size(samples_dot, 2) == size(samples, 2) .and. &
            size(weights_dot) == size(weights)
    end function valid_direction_shapes

    subroutine invalid_status(status, message)
        type(fortsparse_status_t), intent(out) :: status
        character(len=*), intent(in) :: message

        call status_set(status, FORTSPARSE_INVALID_MATRIX, message)
    end subroutine invalid_status

end module fortfem_flux_surface_average
