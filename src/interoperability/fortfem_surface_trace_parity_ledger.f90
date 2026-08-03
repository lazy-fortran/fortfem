module fortfem_surface_trace_parity_ledger
    !! Work/reciprocity ledger for two fixed-topology surface representations.
    !!
    !! A reference trace/response pair and a candidate trace/response pair are
    !! sampled at the same physical points.  The ledger reports the weighted
    !! reference norm, candidate trace error, relative error, and the normalized
    !! reciprocal work defect
    !!
    !!   | Re <reference_trace, candidate_response>_W
    !!       - Re <candidate_trace, reference_response>_W |.
    !!
    !! This is deliberately a neutral FEM/BEM/DtN/PML and fitted/cut/DG/IGA
    !! comparison primitive.  Meshes, kernels, constitutive laws, and physical
    !! normalization remain caller-owned.  The existing typed boundary
    !! contracts provide the fixed-topology/equation/units checks; this module
    !! adds only the paired work ledger that the multi-backend value parity
    !! report does not contain.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_boundary_operator_contract, only: &
        boundary_operator_contract_t, validate_boundary_operator_contract
    implicit none
    private

    public :: evaluate_surface_trace_parity_ledger
    public :: evaluate_surface_trace_parity_ledger_jvp
    public :: evaluate_surface_trace_parity_ledger_vjp

contains

    subroutine evaluate_surface_trace_parity_ledger( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance, reference_norm, &
            absolute_error, relative_error, reciprocity_defect, within_tolerance, status)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        real(dp), intent(out) :: reference_norm, absolute_error, relative_error
        real(dp), intent(out) :: reciprocity_defect
        logical, intent(out) :: within_tolerance
        integer, intent(out) :: status

        real(dp) :: reference_squared, error_squared
        real(dp) :: pairing_one, pairing_two, work_scale

        reference_norm = 0.0_dp
        absolute_error = 0.0_dp
        relative_error = 0.0_dp
        reciprocity_defect = huge(1.0_dp)
        within_tolerance = .false.
        status = 1
        if (.not. valid_inputs( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance)) return

        reference_squared = weighted_squared_norm(reference_trace, weights)
        error_squared = weighted_squared_norm(candidate_trace-reference_trace, weights)
        reference_norm = sqrt(reference_squared)
        absolute_error = sqrt(error_squared)
        relative_error = absolute_error/max(reference_norm, tiny(1.0_dp))
        call work_pairings(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, pairing_one, pairing_two, work_scale)
        reciprocity_defect = abs(pairing_one-pairing_two)/work_scale
        within_tolerance = absolute_error <= absolute_tolerance .or. &
            relative_error <= relative_tolerance
        status = 0
    end subroutine evaluate_surface_trace_parity_ledger

    subroutine evaluate_surface_trace_parity_ledger_jvp( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance, &
            reference_trace_dot, reference_response_dot, candidate_trace_dot, &
            candidate_response_dot, weights_dot, reference_norm_dot, absolute_error_dot, &
            relative_error_dot, reciprocity_defect_dot, status)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        complex(dp), intent(in) :: reference_trace_dot(:, :), reference_response_dot(:, :)
        complex(dp), intent(in) :: candidate_trace_dot(:, :), candidate_response_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: reference_norm_dot, absolute_error_dot
        real(dp), intent(out) :: relative_error_dot, reciprocity_defect_dot
        integer, intent(out) :: status

        real(dp) :: reference_norm, absolute_error, relative_error, reciprocity_defect
        logical :: within_tolerance
        real(dp) :: reference_squared_dot, error_squared_dot
        real(dp) :: pairing_one, pairing_two, pairing_one_dot, pairing_two_dot
        real(dp) :: work_scale, work_scale_dot, defect, defect_dot
        complex(dp) :: reference_difference, reference_difference_dot
        integer :: component, sample

        reference_norm_dot = 0.0_dp
        absolute_error_dot = 0.0_dp
        relative_error_dot = 0.0_dp
        reciprocity_defect_dot = 0.0_dp
        call evaluate_surface_trace_parity_ledger( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance, reference_norm, &
            absolute_error, relative_error, reciprocity_defect, within_tolerance, status)
        if (status /= 0) return
        if (.not. valid_jvp_inputs(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, reference_trace_dot, reference_response_dot, &
            candidate_trace_dot, candidate_response_dot, weights_dot)) then
            status = 1
            return
        end if
        if (reference_norm <= tiny(1.0_dp) .or. absolute_error <= tiny(1.0_dp)) then
            status = 1
            return
        end if

        reference_squared_dot = 0.0_dp
        error_squared_dot = 0.0_dp
        do sample = 1, size(weights)
            do component = 1, size(reference_trace, 1)
                reference_squared_dot = reference_squared_dot + &
                    weights_dot(sample)*abs(reference_trace(component, sample))**2 + &
                    2.0_dp*weights(sample)*real(conjg(reference_trace(component, sample))* &
                    reference_trace_dot(component, sample), dp)
                reference_difference = candidate_trace(component, sample) - &
                    reference_trace(component, sample)
                reference_difference_dot = candidate_trace_dot(component, sample) - &
                    reference_trace_dot(component, sample)
                error_squared_dot = error_squared_dot + &
                    weights_dot(sample)*abs(reference_difference)**2 + &
                    2.0_dp*weights(sample)*real(conjg(reference_difference)* &
                    reference_difference_dot, dp)
            end do
        end do
        reference_norm_dot = reference_squared_dot/(2.0_dp*reference_norm)
        absolute_error_dot = error_squared_dot/(2.0_dp*absolute_error)
        relative_error_dot = (absolute_error_dot*reference_norm - absolute_error* &
            reference_norm_dot)/reference_norm**2

        call work_pairings_jvp(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, reference_trace_dot, reference_response_dot, &
            candidate_trace_dot, candidate_response_dot, weights_dot, pairing_one, &
            pairing_two, pairing_one_dot, pairing_two_dot, work_scale, work_scale_dot)
        defect = pairing_one-pairing_two
        defect_dot = pairing_one_dot-pairing_two_dot
        if (abs(defect) > tiny(1.0_dp)) then
            reciprocity_defect_dot = sign(1.0_dp, defect)* &
                (defect_dot*work_scale-abs(defect)*work_scale_dot)/work_scale**2
        else
            reciprocity_defect_dot = 0.0_dp
        end if
        if (.not. ieee_is_finite(reference_norm_dot) .or. &
            .not. ieee_is_finite(absolute_error_dot) .or. &
            .not. ieee_is_finite(relative_error_dot) .or. &
            .not. ieee_is_finite(reciprocity_defect_dot)) then
            reference_norm_dot = 0.0_dp
            absolute_error_dot = 0.0_dp
            relative_error_dot = 0.0_dp
            reciprocity_defect_dot = 0.0_dp
            status = 1
            return
        end if
        status = 0
    end subroutine evaluate_surface_trace_parity_ledger_jvp

    subroutine evaluate_surface_trace_parity_ledger_vjp( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance, reference_norm_bar, &
            absolute_error_bar, relative_error_bar, reciprocity_defect_bar, &
            reference_trace_bar, reference_response_bar, candidate_trace_bar, &
            candidate_response_bar, weights_bar, status)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        real(dp), intent(in) :: reference_norm_bar, absolute_error_bar
        real(dp), intent(in) :: relative_error_bar, reciprocity_defect_bar
        complex(dp), intent(out) :: reference_trace_bar(:, :), reference_response_bar(:, :)
        complex(dp), intent(out) :: candidate_trace_bar(:, :), candidate_response_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        integer, intent(out) :: status

        real(dp) :: reference_norm, absolute_error, relative_error, reciprocity_defect
        logical :: within_tolerance
        real(dp) :: pairing_one, pairing_two, work_scale
        real(dp) :: reference_bar_effective, error_bar_effective
        real(dp) :: pairing_bar, work_scale_bar, defect
        complex(dp) :: difference
        integer :: component, sample

        reference_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        reference_response_bar = cmplx(0.0_dp, 0.0_dp, dp)
        candidate_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
        candidate_response_bar = cmplx(0.0_dp, 0.0_dp, dp)
        weights_bar = 0.0_dp
        call evaluate_surface_trace_parity_ledger( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance, reference_norm, &
            absolute_error, relative_error, reciprocity_defect, within_tolerance, status)
        if (status /= 0) return
        if (.not. valid_vjp_inputs(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, reference_norm_bar, absolute_error_bar, &
            relative_error_bar, reciprocity_defect_bar, reference_trace_bar, &
            reference_response_bar, candidate_trace_bar, candidate_response_bar, weights_bar)) then
            status = 1
            return
        end if
        if (reference_norm <= tiny(1.0_dp) .or. absolute_error <= tiny(1.0_dp)) then
            status = 1
            return
        end if

        reference_bar_effective = reference_norm_bar - relative_error_bar* &
            absolute_error/reference_norm**2
        error_bar_effective = absolute_error_bar + relative_error_bar/reference_norm
        do sample = 1, size(weights)
            do component = 1, size(reference_trace, 1)
                reference_trace_bar(component, sample) = reference_bar_effective* &
                    weights(sample)*reference_trace(component, sample)/reference_norm
                weights_bar(sample) = weights_bar(sample) + reference_bar_effective* &
                    abs(reference_trace(component, sample))**2/(2.0_dp*reference_norm)
                difference = candidate_trace(component, sample) - &
                    reference_trace(component, sample)
                candidate_trace_bar(component, sample) = error_bar_effective* &
                    weights(sample)*difference/absolute_error
                reference_trace_bar(component, sample) = reference_trace_bar( &
                    component, sample) - error_bar_effective*weights(sample)* &
                    difference/absolute_error
                weights_bar(sample) = weights_bar(sample) + error_bar_effective* &
                    abs(difference)**2/(2.0_dp*absolute_error)
            end do
        end do

        call work_pairings(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, pairing_one, pairing_two, work_scale)
        defect = pairing_one-pairing_two
        if (abs(defect) > tiny(1.0_dp) .and. abs(reciprocity_defect_bar) > 0.0_dp) then
            pairing_bar = sign(1.0_dp, defect)*reciprocity_defect_bar/work_scale
            work_scale_bar = -abs(defect)*reciprocity_defect_bar/work_scale**2
            if (abs(pairing_one) >= max(1.0_dp, abs(pairing_two))) then
                pairing_bar = pairing_bar + sign(1.0_dp, pairing_one)*work_scale_bar
            else if (abs(pairing_two) >= 1.0_dp) then
                pairing_bar = pairing_bar - sign(1.0_dp, pairing_two)*work_scale_bar
            end if
            do sample = 1, size(weights)
                do component = 1, size(reference_trace, 1)
                    reference_trace_bar(component, sample) = &
                        reference_trace_bar(component, sample) + pairing_bar*weights(sample)* &
                        conjg(candidate_response(component, sample))
                    candidate_response_bar(component, sample) = &
                        candidate_response_bar(component, sample) + pairing_bar*weights(sample)* &
                        conjg(reference_trace(component, sample))
                    candidate_trace_bar(component, sample) = &
                        candidate_trace_bar(component, sample) - pairing_bar*weights(sample)* &
                        conjg(reference_response(component, sample))
                    reference_response_bar(component, sample) = &
                        reference_response_bar(component, sample) - pairing_bar*weights(sample)* &
                        conjg(candidate_trace(component, sample))
                    weights_bar(sample) = weights_bar(sample) + pairing_bar*real( &
                        reference_trace(component, sample)*candidate_response(component, sample) - &
                        candidate_trace(component, sample)*reference_response(component, sample), dp)
                end do
            end do
        end if
        if (.not. finite_complex_2d(reference_trace_bar) .or. &
            .not. finite_complex_2d(reference_response_bar) .or. &
            .not. finite_complex_2d(candidate_trace_bar) .or. &
            .not. finite_complex_2d(candidate_response_bar) .or. &
            .not. finite_real(weights_bar)) then
            reference_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
            reference_response_bar = cmplx(0.0_dp, 0.0_dp, dp)
            candidate_trace_bar = cmplx(0.0_dp, 0.0_dp, dp)
            candidate_response_bar = cmplx(0.0_dp, 0.0_dp, dp)
            weights_bar = 0.0_dp
            status = 1
            return
        end if
        status = 0
    end subroutine evaluate_surface_trace_parity_ledger_vjp

    logical function valid_inputs( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, contracts, absolute_tolerance, relative_tolerance) result(valid)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        type(boundary_operator_contract_t), intent(in) :: contracts(:)
        real(dp), intent(in) :: absolute_tolerance, relative_tolerance
        integer :: contract_status

        valid = .false.
        if (size(reference_trace, 1) < 1 .or. size(reference_trace, 2) < 1 .or. &
            any(shape(reference_response) /= shape(reference_trace)) .or. &
            any(shape(candidate_trace) /= shape(reference_trace)) .or. &
            any(shape(candidate_response) /= shape(reference_trace)) .or. &
            size(weights) /= size(reference_trace, 2) .or. size(contracts) /= 2) return
        if (.not. finite_complex_2d(reference_trace) .or. &
            .not. finite_complex_2d(reference_response) .or. &
            .not. finite_complex_2d(candidate_trace) .or. &
            .not. finite_complex_2d(candidate_response) .or. &
            .not. finite_real(weights) .or. any(weights <= 0.0_dp)) return
        if (.not. ieee_is_finite(absolute_tolerance) .or. &
            .not. ieee_is_finite(relative_tolerance) .or. &
            absolute_tolerance < 0.0_dp .or. relative_tolerance < 0.0_dp) return
        if (.not. validate_boundary_operator_contract(contracts(1), contract_status)) return
        if (.not. validate_boundary_operator_contract(contracts(2), contract_status)) return
        if (contracts(1)%backend_kind == contracts(2)%backend_kind .or. &
            trim(contracts(1)%topology_id) /= trim(contracts(2)%topology_id) .or. &
            trim(contracts(1)%equation) /= trim(contracts(2)%equation) .or. &
            trim(contracts(1)%space) /= trim(contracts(2)%space) .or. &
            trim(contracts(1)%units) /= trim(contracts(2)%units) .or. &
            trim(contracts(1)%normalization) /= trim(contracts(2)%normalization)) return
        valid = .true.
    end function valid_inputs

    logical function valid_jvp_inputs( &
            reference_trace, reference_response, candidate_trace, candidate_response, &
            weights, reference_trace_dot, reference_response_dot, candidate_trace_dot, &
            candidate_response_dot, weights_dot) result(valid)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: reference_trace_dot(:, :), reference_response_dot(:, :)
        complex(dp), intent(in) :: candidate_trace_dot(:, :), candidate_response_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)

        valid = .not. any(shape(reference_trace_dot) /= shape(reference_trace)) .and. &
            .not. any(shape(reference_response_dot) /= shape(reference_response)) .and. &
            .not. any(shape(candidate_trace_dot) /= shape(candidate_trace)) .and. &
            .not. any(shape(candidate_response_dot) /= shape(candidate_response)) .and. &
            size(weights_dot) == size(weights) .and. finite_complex_2d(reference_trace_dot) .and. &
            finite_complex_2d(reference_response_dot) .and. finite_complex_2d(candidate_trace_dot) .and. &
            finite_complex_2d(candidate_response_dot) .and. finite_real(weights_dot)
    end function valid_jvp_inputs

    logical function valid_vjp_inputs( &
            reference_trace, reference_response, candidate_trace, candidate_response, weights, &
            reference_norm_bar, absolute_error_bar, relative_error_bar, reciprocity_defect_bar, &
            reference_trace_bar, reference_response_bar, candidate_trace_bar, &
            candidate_response_bar, weights_bar) result(valid)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:), reference_norm_bar, absolute_error_bar
        real(dp), intent(in) :: relative_error_bar, reciprocity_defect_bar
        complex(dp), intent(in) :: reference_trace_bar(:, :), reference_response_bar(:, :)
        complex(dp), intent(in) :: candidate_trace_bar(:, :), candidate_response_bar(:, :)
        real(dp), intent(in) :: weights_bar(:)

        valid = ieee_is_finite(reference_norm_bar) .and. ieee_is_finite(absolute_error_bar) .and. &
            ieee_is_finite(relative_error_bar) .and. ieee_is_finite(reciprocity_defect_bar) .and. &
            all(shape(reference_trace_bar) == shape(reference_trace)) .and. &
            all(shape(reference_response_bar) == shape(reference_response)) .and. &
            all(shape(candidate_trace_bar) == shape(candidate_trace)) .and. &
            all(shape(candidate_response_bar) == shape(candidate_response)) .and. &
            size(weights_bar) == size(weights)
    end function valid_vjp_inputs

    subroutine work_pairings(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, pairing_one, pairing_two, work_scale)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: pairing_one, pairing_two, work_scale

        real(dp) :: weighted_reference, weighted_candidate

        weighted_reference = 0.0_dp
        weighted_candidate = 0.0_dp
        pairing_one = real(sum(spread(weights, 1, size(reference_trace, 1))* &
            reference_trace*candidate_response), dp)
        pairing_two = real(sum(spread(weights, 1, size(reference_trace, 1))* &
            candidate_trace*reference_response), dp)
        work_scale = max(1.0_dp, abs(pairing_one), abs(pairing_two))
    end subroutine work_pairings

    subroutine work_pairings_jvp(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, reference_trace_dot, reference_response_dot, &
            candidate_trace_dot, candidate_response_dot, weights_dot, pairing_one, &
            pairing_two, pairing_one_dot, pairing_two_dot, work_scale, work_scale_dot)
        complex(dp), intent(in) :: reference_trace(:, :), reference_response(:, :)
        complex(dp), intent(in) :: candidate_trace(:, :), candidate_response(:, :)
        real(dp), intent(in) :: weights(:)
        complex(dp), intent(in) :: reference_trace_dot(:, :), reference_response_dot(:, :)
        complex(dp), intent(in) :: candidate_trace_dot(:, :), candidate_response_dot(:, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: pairing_one, pairing_two, pairing_one_dot
        real(dp), intent(out) :: pairing_two_dot, work_scale, work_scale_dot

        complex(dp) :: first_dot, second_dot
        call work_pairings(reference_trace, reference_response, candidate_trace, &
            candidate_response, weights, pairing_one, pairing_two, work_scale)
        first_dot = sum(spread(weights_dot, 1, size(reference_trace, 1))* &
            reference_trace*candidate_response + spread(weights, 1, size(reference_trace, 1))* &
            (reference_trace_dot*candidate_response + reference_trace*candidate_response_dot))
        second_dot = sum(spread(weights_dot, 1, size(reference_trace, 1))* &
            candidate_trace*reference_response + spread(weights, 1, size(reference_trace, 1))* &
            (candidate_trace_dot*reference_response + candidate_trace*reference_response_dot))
        pairing_one_dot = real(first_dot, dp)
        pairing_two_dot = real(second_dot, dp)
        work_scale_dot = 0.0_dp
        if (abs(pairing_one) >= max(1.0_dp, abs(pairing_two))) then
            work_scale_dot = sign(1.0_dp, pairing_one)*pairing_one_dot
        else if (abs(pairing_two) >= 1.0_dp) then
            work_scale_dot = sign(1.0_dp, pairing_two)*pairing_two_dot
        end if
    end subroutine work_pairings_jvp

    real(dp) function weighted_squared_norm(values, weights) result(norm_squared)
        complex(dp), intent(in) :: values(:, :)
        real(dp), intent(in) :: weights(:)

        norm_squared = sum(spread(weights, 1, size(values, 1))*abs(values)**2)
    end function weighted_squared_norm

    logical function finite_complex_2d(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_2d

    logical function finite_real(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_real

end module fortfem_surface_trace_parity_ledger
