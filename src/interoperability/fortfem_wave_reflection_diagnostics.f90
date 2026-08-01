module fortfem_wave_reflection_diagnostics
    !! Weighted complex error and reflection metrics for open-boundary oracles.
    !!
    !! The routines operate on caller-owned physical samples.  They do not
    !! assume a mesh, coordinate system, boundary condition, or wave model.
    !! A reflection coefficient is the weighted norm of total-minus-incident
    !! divided by the weighted incident norm.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: evaluate_weighted_complex_error
    public :: evaluate_weighted_complex_error_jvp
    public :: evaluate_weighted_complex_error_vjp
    public :: evaluate_weighted_reflection_coefficient
    public :: evaluate_weighted_reflection_coefficient_jvp
    public :: evaluate_weighted_reflection_coefficient_vjp

contains

    subroutine evaluate_weighted_complex_error( &
            reference, candidate, weights, absolute_error, relative_error, &
            status)
        complex(dp), intent(in) :: reference(:, :), candidate(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: absolute_error, relative_error
        integer, intent(out) :: status

        real(dp) :: error_squared, reference_squared
        integer :: component, sample

        absolute_error = huge(1.0_dp)
        relative_error = huge(1.0_dp)
        status = 1
        if (.not. compatible_shapes(reference, candidate, weights)) return
        if (.not. finite_complex_array(reference)) return
        if (.not. finite_complex_array(candidate)) return
        if (.not. finite_real_array(weights)) return
        if (any(weights <= 0.0_dp)) return
        error_squared = 0.0_dp
        reference_squared = 0.0_dp
        do sample = 1, size(weights)
            do component = 1, size(reference, 1)
                error_squared = error_squared + weights(sample)*abs( &
                    candidate(component, sample) - reference(component, sample))**2
                reference_squared = reference_squared + weights(sample)*abs( &
                    reference(component, sample))**2
            end do
        end do
        absolute_error = sqrt(error_squared)
        relative_error = absolute_error/max( &
            sqrt(reference_squared), tiny(1.0_dp))
        status = 0
    end subroutine evaluate_weighted_complex_error

    subroutine evaluate_weighted_complex_error_jvp( &
            reference, candidate, weights, reference_dot, candidate_dot, &
            weights_dot, absolute_error_dot, relative_error_dot, status)
        complex(dp), intent(in) :: reference(:, :), candidate(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: reference_dot(:, :), candidate_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: absolute_error_dot, relative_error_dot
        integer, intent(out) :: status

        real(dp) :: absolute_error, relative_error, error_squared_dot
        real(dp) :: reference_norm, reference_norm_dot
        complex(dp) :: difference, difference_dot
        integer :: component, sample

        absolute_error_dot = 0.0_dp
        relative_error_dot = 0.0_dp
        status = 1
        call validate_error_inputs( &
            reference, candidate, weights, reference_dot, candidate_dot, &
            weights_dot, status)
        if (status /= 0) return
        call evaluate_weighted_complex_error( &
            reference, candidate, weights, absolute_error, relative_error, status)
        if (status /= 0) return
        reference_norm = weighted_norm(reference, weights)
        if (absolute_error <= tiny(1.0_dp) .or. &
            reference_norm <= tiny(1.0_dp)) return
        error_squared_dot = 0.0_dp
        reference_norm_dot = 0.0_dp
        do sample = 1, size(weights)
            do component = 1, size(reference, 1)
                difference = candidate(component, sample) - &
                    reference(component, sample)
                difference_dot = candidate_dot(component, sample) - &
                    reference_dot(component, sample)
                error_squared_dot = error_squared_dot + weights_dot(sample)* &
                    abs(difference)**2 + 2.0_dp*weights(sample)*real( &
                    conjg(difference)*difference_dot, dp)
                reference_norm_dot = reference_norm_dot + 0.5_dp/ &
                    reference_norm*(weights_dot(sample)*abs( &
                    reference(component, sample))**2 + 2.0_dp*weights(sample)* &
                    real(conjg(reference(component, sample))* &
                    reference_dot(component, sample), dp))
            end do
        end do
        absolute_error_dot = error_squared_dot/(2.0_dp*absolute_error)
        relative_error_dot = absolute_error_dot/reference_norm - &
            absolute_error*reference_norm_dot/reference_norm**2
        status = 0
    end subroutine evaluate_weighted_complex_error_jvp

    subroutine evaluate_weighted_complex_error_vjp( &
            reference, candidate, weights, absolute_error_bar, &
            relative_error_bar, reference_bar, candidate_bar, weights_bar, &
            status)
        complex(dp), intent(in) :: reference(:, :), candidate(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: absolute_error_bar, relative_error_bar
        complex(dp), intent(out) :: reference_bar(:, :), candidate_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        integer, intent(out) :: status

        real(dp) :: absolute_error, relative_error, reference_norm
        real(dp) :: absolute_partial, reference_partial
        complex(dp) :: difference
        integer :: component, sample

        reference_bar = cmplx(0.0_dp, 0.0_dp, dp)
        candidate_bar = cmplx(0.0_dp, 0.0_dp, dp)
        weights_bar = 0.0_dp
        status = 1
        if (.not. ieee_is_finite(absolute_error_bar)) return
        if (.not. ieee_is_finite(relative_error_bar)) return
        call evaluate_weighted_complex_error( &
            reference, candidate, weights, absolute_error, relative_error, status)
        if (status /= 0) return
        if (any(shape(reference_bar) /= shape(reference))) return
        if (any(shape(candidate_bar) /= shape(candidate))) return
        if (size(weights_bar) /= size(weights)) return
        reference_norm = weighted_norm(reference, weights)
        if (absolute_error <= tiny(1.0_dp) .or. &
            reference_norm <= tiny(1.0_dp)) then
            status = 1
            return
        end if
        absolute_partial = absolute_error_bar + relative_error_bar/reference_norm
        reference_partial = -relative_error_bar*absolute_error/reference_norm**2
        do sample = 1, size(weights)
            do component = 1, size(reference, 1)
                difference = candidate(component, sample) - &
                    reference(component, sample)
                candidate_bar(component, sample) = absolute_partial*weights(sample)* &
                    difference/absolute_error
                reference_bar(component, sample) = -candidate_bar( &
                    component, sample) + reference_partial*weights(sample)* &
                    reference(component, sample)/reference_norm
                weights_bar(sample) = weights_bar(sample) + &
                    absolute_partial*abs(difference)**2/(2.0_dp*absolute_error) + &
                    reference_partial*abs(reference(component, sample))**2/( &
                    2.0_dp*reference_norm)
            end do
        end do
        status = 0
    end subroutine evaluate_weighted_complex_error_vjp

    subroutine evaluate_weighted_reflection_coefficient( &
            incident, total, weights, coefficient, status)
        complex(dp), intent(in) :: incident(:, :), total(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: coefficient
        integer, intent(out) :: status

        real(dp) :: absolute_error

        coefficient = huge(1.0_dp)
        call evaluate_weighted_complex_error( &
            incident, total, weights, absolute_error, coefficient, status)
    end subroutine evaluate_weighted_reflection_coefficient

    subroutine evaluate_weighted_reflection_coefficient_jvp( &
            incident, total, weights, incident_dot, total_dot, weights_dot, &
            coefficient_dot, status)
        complex(dp), intent(in) :: incident(:, :), total(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: incident_dot(:, :), total_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: coefficient_dot
        integer, intent(out) :: status

        real(dp) :: absolute_error_dot

        coefficient_dot = 0.0_dp
        call evaluate_weighted_complex_error_jvp( &
            incident, total, weights, incident_dot, total_dot, weights_dot, &
            absolute_error_dot, coefficient_dot, status)
    end subroutine evaluate_weighted_reflection_coefficient_jvp

    subroutine evaluate_weighted_reflection_coefficient_vjp( &
            incident, total, weights, coefficient_bar, incident_bar, total_bar, &
            weights_bar, status)
        complex(dp), intent(in) :: incident(:, :), total(:, :)
        real(dp), intent(in) :: weights(:), coefficient_bar
        complex(dp), intent(out) :: incident_bar(:, :), total_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        integer, intent(out) :: status

        call evaluate_weighted_complex_error_vjp( &
            incident, total, weights, 0.0_dp, coefficient_bar, incident_bar, &
            total_bar, weights_bar, status)
    end subroutine evaluate_weighted_reflection_coefficient_vjp

    subroutine validate_error_inputs( &
            reference, candidate, weights, reference_dot, candidate_dot, &
            weights_dot, status)
        complex(dp), intent(in) :: reference(:, :), candidate(:, :)
        real(dp), intent(in) :: weights(:), weights_dot(:)
        complex(dp), intent(in) :: reference_dot(:, :), candidate_dot(:, :)
        integer, intent(out) :: status

        status = 1
        if (.not. compatible_shapes(reference, candidate, weights)) return
        if (any(shape(reference_dot) /= shape(reference))) return
        if (any(shape(candidate_dot) /= shape(candidate))) return
        if (size(weights_dot) /= size(weights)) return
        if (.not. finite_complex_array(reference)) return
        if (.not. finite_complex_array(candidate)) return
        if (.not. finite_complex_array(reference_dot)) return
        if (.not. finite_complex_array(candidate_dot)) return
        if (.not. finite_real_array(weights)) return
        if (.not. finite_real_array(weights_dot)) return
        if (any(weights <= 0.0_dp)) return
        status = 0
    end subroutine validate_error_inputs

    logical function compatible_shapes(reference, candidate, weights) result(valid)
        complex(dp), intent(in) :: reference(:, :), candidate(:, :)
        real(dp), intent(in) :: weights(:)

        valid = .false.
        if (size(reference, 1) < 1 .or. size(reference, 2) < 1) return
        if (any(shape(candidate) /= shape(reference))) return
        if (size(weights) /= size(reference, 2)) return
        valid = .true.
    end function compatible_shapes

    real(dp) function weighted_norm(values, weights) result(norm)
        complex(dp), intent(in) :: values(:, :)
        real(dp), intent(in) :: weights(:)

        integer :: component, sample

        norm = 0.0_dp
        do sample = 1, size(weights)
            do component = 1, size(values, 1)
                norm = norm + weights(sample)*abs(values(component, sample))**2
            end do
        end do
        norm = sqrt(norm)
    end function weighted_norm

    logical function finite_complex_array(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_array

    logical function finite_real_array(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_real_array

end module fortfem_wave_reflection_diagnostics
