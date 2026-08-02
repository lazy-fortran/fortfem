module fortfem_resistive_mhd_branch_history
    !! Closure-neutral continuation branch-history diagnostics.
    !!
    !! A caller supplies a fixed-topology sequence of continuation parameters,
    !! states, residual samples, energy samples, and integer branch labels.
    !! This module validates strict parameter monotonicity, compares two
    !! histories for deterministic branch multiplicity and hysteresis, and
    !! supplies a differentiable path metric.  It does not select a branch,
    !! implement a nonlinear solver, or infer a resistive-MHD closure.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: resistive_mhd_branch_history_t
        !! Caller-owned samples along one strictly monotone path.
        real(dp), allocatable :: parameter(:)
        real(dp), allocatable :: state(:, :)
        real(dp), allocatable :: residual(:, :)
        real(dp), allocatable :: energy(:)
        integer, allocatable :: branch_id(:)
        integer :: monotonic_direction = 1
    end type resistive_mhd_branch_history_t

    type, public :: resistive_mhd_branch_diagnostics_t
        !! Discrete topology labels plus continuous hysteresis magnitudes.
        integer :: branch_multiplicity = 0
        integer :: hysteresis_sample_count = 0
        logical :: hysteresis_detected = .false.
        real(dp) :: max_state_hysteresis = 0.0_dp
        real(dp) :: max_residual_hysteresis = 0.0_dp
        real(dp) :: max_energy_hysteresis = 0.0_dp
    end type resistive_mhd_branch_diagnostics_t

    public :: initialize_resistive_mhd_branch_history
    public :: validate_resistive_mhd_branch_history
    public :: evaluate_resistive_mhd_branch_diagnostics
    public :: compare_resistive_mhd_branch_histories
    public :: evaluate_resistive_mhd_branch_path_metric
    public :: evaluate_resistive_mhd_branch_path_metric_jvp
    public :: evaluate_resistive_mhd_branch_path_metric_vjp

contains

    subroutine initialize_resistive_mhd_branch_history( &
            parameter, state, residual, energy, branch_id, history, status, &
            monotonic_direction)
        real(dp), intent(in) :: parameter(:), state(:, :), residual(:, :), energy(:)
        integer, intent(in) :: branch_id(:)
        type(resistive_mhd_branch_history_t), intent(out) :: history
        type(fortsparse_status_t), intent(out) :: status
        integer, intent(in), optional :: monotonic_direction
        integer :: direction

        history = resistive_mhd_branch_history_t()
        direction = 1
        if (present(monotonic_direction)) direction = monotonic_direction
        if (.not. valid_sample_arrays(parameter, state, residual, energy, branch_id, &
                direction)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch history has incompatible samples")
            return
        end if
        history%parameter = parameter
        history%state = state
        history%residual = residual
        history%energy = energy
        history%branch_id = branch_id
        history%monotonic_direction = direction
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_resistive_mhd_branch_history

    subroutine validate_resistive_mhd_branch_history(history, status)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        type(fortsparse_status_t), intent(out) :: status

        if (.not. allocated_history_arrays(history)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch history is not allocated")
            return
        end if
        if (.not. valid_sample_arrays(history%parameter, history%state, &
                history%residual, history%energy, history%branch_id, &
                history%monotonic_direction)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch history is invalid")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_resistive_mhd_branch_history

    subroutine evaluate_resistive_mhd_branch_diagnostics(history, diagnostics, status)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        type(resistive_mhd_branch_diagnostics_t), intent(out) :: diagnostics
        type(fortsparse_status_t), intent(out) :: status

        diagnostics = resistive_mhd_branch_diagnostics_t()
        call validate_resistive_mhd_branch_history(history, status)
        if (status%code /= FORTSPARSE_OK) return
        diagnostics%branch_multiplicity = count_unique_labels(history%branch_id)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_resistive_mhd_branch_diagnostics

    subroutine compare_resistive_mhd_branch_histories( &
            first, second, tolerance, diagnostics, status)
        type(resistive_mhd_branch_history_t), intent(in) :: first, second
        real(dp), intent(in) :: tolerance
        type(resistive_mhd_branch_diagnostics_t), intent(out) :: diagnostics
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample, label
        real(dp) :: state_difference, residual_difference, energy_difference

        diagnostics = resistive_mhd_branch_diagnostics_t()
        call validate_resistive_mhd_branch_history(first, status)
        if (status%code /= FORTSPARSE_OK) return
        call validate_resistive_mhd_branch_history(second, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp .or. &
                size(first%parameter) /= size(second%parameter) .or. &
                size(first%state, 1) /= size(second%state, 1) .or. &
                first%monotonic_direction /= second%monotonic_direction .or. &
                any(abs(first%parameter - second%parameter) > tolerance)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch histories cannot be compared")
            return
        end if

        diagnostics%branch_multiplicity = count_unique_union(first%branch_id, second%branch_id)
        do sample = 1, size(first%branch_id)
            state_difference = maxval(abs(first%state(:, sample) - second%state(:, sample)))
            residual_difference = maxval(abs(first%residual(:, sample) - &
                second%residual(:, sample)))
            energy_difference = abs(first%energy(sample) - second%energy(sample))
            diagnostics%max_state_hysteresis = max( &
                diagnostics%max_state_hysteresis, state_difference)
            diagnostics%max_residual_hysteresis = max( &
                diagnostics%max_residual_hysteresis, residual_difference)
            diagnostics%max_energy_hysteresis = max( &
                diagnostics%max_energy_hysteresis, energy_difference)
            if (state_difference > tolerance .or. residual_difference > tolerance .or. &
                    energy_difference > tolerance) then
                diagnostics%hysteresis_detected = .true.
                diagnostics%hysteresis_sample_count = &
                    diagnostics%hysteresis_sample_count + 1
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine compare_resistive_mhd_branch_histories

    subroutine evaluate_resistive_mhd_branch_path_metric(history, metric, status)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        real(dp), intent(out) :: metric
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: left_integrand, right_integrand, parameter_step

        metric = 0.0_dp
        call validate_resistive_mhd_branch_history(history, status)
        if (status%code /= FORTSPARSE_OK) return
        do sample = 1, size(history%parameter) - 1
            left_integrand = sample_integrand(history, sample)
            right_integrand = sample_integrand(history, sample + 1)
            parameter_step = history%parameter(sample + 1) - history%parameter(sample)
            metric = metric + 0.5_dp*(left_integrand + right_integrand)*parameter_step
        end do
        if (.not. ieee_is_finite(metric)) then
            metric = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch path metric is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_resistive_mhd_branch_path_metric

    subroutine evaluate_resistive_mhd_branch_path_metric_jvp( &
            history, parameter_dot, state_dot, residual_dot, energy_dot, metric_dot, status)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        real(dp), intent(in) :: parameter_dot(:), state_dot(:, :), residual_dot(:, :)
        real(dp), intent(in) :: energy_dot(:)
        real(dp), intent(out) :: metric_dot
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: left_integrand, right_integrand, left_dot, right_dot
        real(dp) :: parameter_step, parameter_step_dot

        metric_dot = 0.0_dp
        call validate_resistive_mhd_branch_history(history, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_directions(history, parameter_dot, state_dot, residual_dot, &
                energy_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch path metric JVP has incompatible directions")
            return
        end if
        do sample = 1, size(history%parameter) - 1
            left_integrand = sample_integrand(history, sample)
            right_integrand = sample_integrand(history, sample + 1)
            left_dot = sample_integrand_dot(history, parameter_dot, state_dot, &
                residual_dot, energy_dot, sample)
            right_dot = sample_integrand_dot(history, parameter_dot, state_dot, &
                residual_dot, energy_dot, sample + 1)
            parameter_step = history%parameter(sample + 1) - history%parameter(sample)
            parameter_step_dot = parameter_dot(sample + 1) - parameter_dot(sample)
            metric_dot = metric_dot + 0.5_dp*(left_dot + right_dot)*parameter_step + &
                0.5_dp*(left_integrand + right_integrand)*parameter_step_dot
        end do
        if (.not. ieee_is_finite(metric_dot)) then
            metric_dot = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch path metric JVP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_resistive_mhd_branch_path_metric_jvp

    subroutine evaluate_resistive_mhd_branch_path_metric_vjp( &
            history, metric_bar, parameter_bar, state_bar, residual_bar, energy_bar, status)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        real(dp), intent(in) :: metric_bar
        real(dp), intent(out) :: parameter_bar(:), state_bar(:, :), residual_bar(:, :)
        real(dp), intent(out) :: energy_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: sample
        real(dp) :: left_integrand, right_integrand, left_bar, right_bar
        real(dp) :: parameter_step, interval_bar

        parameter_bar = 0.0_dp
        state_bar = 0.0_dp
        residual_bar = 0.0_dp
        energy_bar = 0.0_dp
        call validate_resistive_mhd_branch_history(history, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(metric_bar) .or. size(parameter_bar) /= &
                size(history%parameter) .or. size(state_bar, 1) /= size(history%state, 1) .or. &
                size(state_bar, 2) /= size(history%state, 2) .or. &
                size(residual_bar, 1) /= size(history%residual, 1) .or. &
                size(residual_bar, 2) /= size(history%residual, 2) .or. &
                size(energy_bar) /= size(history%energy)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch path metric VJP has incompatible cotangents")
            return
        end if
        do sample = 1, size(history%parameter) - 1
            left_integrand = sample_integrand(history, sample)
            right_integrand = sample_integrand(history, sample + 1)
            parameter_step = history%parameter(sample + 1) - history%parameter(sample)
            left_bar = 0.5_dp*metric_bar*parameter_step
            right_bar = left_bar
            parameter_bar(sample) = parameter_bar(sample) - &
                0.5_dp*metric_bar*(left_integrand + right_integrand)
            parameter_bar(sample + 1) = parameter_bar(sample + 1) + &
                0.5_dp*metric_bar*(left_integrand + right_integrand)
            state_bar(:, sample) = state_bar(:, sample) + left_bar*history%state(:, sample)
            state_bar(:, sample + 1) = state_bar(:, sample + 1) + &
                right_bar*history%state(:, sample + 1)
            residual_bar(:, sample) = residual_bar(:, sample) + &
                left_bar*history%residual(:, sample)
            residual_bar(:, sample + 1) = residual_bar(:, sample + 1) + &
                right_bar*history%residual(:, sample + 1)
            energy_bar(sample) = energy_bar(sample) + left_bar
            energy_bar(sample + 1) = energy_bar(sample + 1) + right_bar
        end do
        if (.not. all(ieee_is_finite(parameter_bar)) .or. &
                .not. all(ieee_is_finite(state_bar)) .or. &
                .not. all(ieee_is_finite(residual_bar)) .or. &
                .not. all(ieee_is_finite(energy_bar))) then
            parameter_bar = 0.0_dp
            state_bar = 0.0_dp
            residual_bar = 0.0_dp
            energy_bar = 0.0_dp
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "resistive-MHD branch path metric VJP is non-finite")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_resistive_mhd_branch_path_metric_vjp

    integer function count_unique_labels(labels) result(count)
        integer, intent(in) :: labels(:)
        integer :: index, previous
        logical :: seen

        count = 0
        do index = 1, size(labels)
            seen = .false.
            do previous = 1, index - 1
                if (labels(previous) == labels(index)) then
                    seen = .true.
                    exit
                end if
            end do
            if (.not. seen) count = count + 1
        end do
    end function count_unique_labels

    integer function count_unique_union(first, second) result(count)
        integer, intent(in) :: first(:), second(:)
        integer, allocatable :: labels(:)

        labels = [first, second]
        count = count_unique_labels(labels)
    end function count_unique_union

    logical function allocated_history_arrays(history) result(valid)
        type(resistive_mhd_branch_history_t), intent(in) :: history

        valid = allocated(history%parameter) .and. allocated(history%state) .and. &
            allocated(history%residual) .and. allocated(history%energy) .and. &
            allocated(history%branch_id)
    end function allocated_history_arrays

    logical function valid_sample_arrays( &
            parameter, state, residual, energy, branch_id, direction) result(valid)
        real(dp), intent(in) :: parameter(:), state(:, :), residual(:, :), energy(:)
        integer, intent(in) :: branch_id(:), direction
        integer :: sample

        valid = direction == 1 .or. direction == -1
        if (.not. valid) return
        valid = size(parameter) > 0 .and. size(state, 2) == size(parameter) .and. &
            size(residual, 2) == size(parameter) .and. size(state, 1) > 0 .and. &
            size(residual, 1) == size(state, 1) .and. size(energy) == size(parameter) .and. &
            size(branch_id) == size(parameter)
        if (.not. valid) return
        valid = all(ieee_is_finite(parameter)) .and. all(ieee_is_finite(state)) .and. &
            all(ieee_is_finite(residual)) .and. all(ieee_is_finite(energy)) .and. &
            all(branch_id >= 1)
        if (.not. valid) return
        do sample = 1, size(parameter) - 1
            if (direction == 1) then
                valid = parameter(sample + 1) > parameter(sample)
            else
                valid = parameter(sample + 1) < parameter(sample)
            end if
            if (.not. valid) return
        end do
    end function valid_sample_arrays

    logical function valid_directions( &
            history, parameter_dot, state_dot, residual_dot, energy_dot) result(valid)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        real(dp), intent(in) :: parameter_dot(:), state_dot(:, :), residual_dot(:, :)
        real(dp), intent(in) :: energy_dot(:)

        valid = size(parameter_dot) == size(history%parameter) .and. &
            size(state_dot, 1) == size(history%state, 1) .and. &
            size(state_dot, 2) == size(history%state, 2) .and. &
            size(residual_dot, 1) == size(history%residual, 1) .and. &
            size(residual_dot, 2) == size(history%residual, 2) .and. &
            size(energy_dot) == size(history%energy)
        if (.not. valid) return
        valid = all(ieee_is_finite(parameter_dot)) .and. all(ieee_is_finite(state_dot)) .and. &
            all(ieee_is_finite(residual_dot)) .and. all(ieee_is_finite(energy_dot))
    end function valid_directions

    real(dp) function sample_integrand(history, sample) result(value)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        integer, intent(in) :: sample

        value = history%energy(sample) + 0.5_dp* &
            (dot_product(history%state(:, sample), history%state(:, sample)) + &
            dot_product(history%residual(:, sample), history%residual(:, sample)))
    end function sample_integrand

    real(dp) function sample_integrand_dot( &
            history, parameter_dot, state_dot, residual_dot, energy_dot, sample) result(value)
        type(resistive_mhd_branch_history_t), intent(in) :: history
        real(dp), intent(in) :: parameter_dot(:), state_dot(:, :), residual_dot(:, :)
        real(dp), intent(in) :: energy_dot(:)
        integer, intent(in) :: sample

        associate (unused_parameter_dot => parameter_dot)
            value = energy_dot(sample) + &
                dot_product(history%state(:, sample), state_dot(:, sample)) + &
                dot_product(history%residual(:, sample), residual_dot(:, sample))
        end associate
    end function sample_integrand_dot

end module fortfem_resistive_mhd_branch_history
