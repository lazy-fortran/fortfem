module fortfem_oracle_performance_ledger
    !! Differentiable weighted performance ledger for external oracle runs.
    !!
    !! The input records remain the metadata-owned `oracle_manifest_t` objects;
    !! this adapter extracts only continuous benchmark metrics: six phase times,
    !! total time, peak memory, and four declared tolerances.  It reports their
    !! weighted means and exact real JVP/VJP products.  Code identities may
    !! differ, but case/checksum/normalization/hardware provenance, success,
    !! warmup count, and repetition count must agree before metrics are mixed.
    !! No external executable, reader, or license is selected here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortfem_oracle_manifest, only: &
        oracle_manifest_t, validate_oracle_manifest
    implicit none
    private

    public :: evaluate_oracle_performance_ledger
    public :: evaluate_oracle_performance_ledger_jvp
    public :: evaluate_oracle_performance_ledger_vjp

contains

    subroutine evaluate_oracle_performance_ledger( &
            manifests, weights, phase_means, total_mean, memory_mean, &
            tolerance_means, repetition_count, status)
        type(oracle_manifest_t), intent(in) :: manifests(:)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: phase_means(:), total_mean, memory_mean
        real(dp), intent(out) :: tolerance_means(:)
        integer, intent(out) :: repetition_count, status

        real(dp) :: phase_values(6, size(manifests))
        real(dp) :: total_values(size(manifests)), memory_values(size(manifests))
        real(dp) :: tolerance_values(4, size(manifests))

        phase_means = 0.0_dp
        total_mean = 0.0_dp
        memory_mean = 0.0_dp
        tolerance_means = 0.0_dp
        repetition_count = 0
        status = 1
        if (size(phase_means) /= 6 .or. size(tolerance_means) /= 4) return
        if (.not. validate_records(manifests, weights, repetition_count)) return
        call extract_values(manifests, phase_values, total_values, memory_values, &
            tolerance_values)
        call weighted_means(phase_values, total_values, memory_values, &
            tolerance_values, &
            weights, phase_means, total_mean, memory_mean, tolerance_means)
        if (.not. finite_real(phase_means) .or. &
            .not. ieee_is_finite(total_mean) .or. &
            .not. ieee_is_finite(memory_mean) .or. &
            .not. finite_real(tolerance_means)) return
        status = 0
    end subroutine evaluate_oracle_performance_ledger

    subroutine evaluate_oracle_performance_ledger_jvp( &
            manifests, weights, phase_dot, total_dot, memory_dot, tolerance_dot, &
            weights_dot, phase_mean_dot, total_mean_dot, memory_mean_dot, &
            tolerance_mean_dot, status)
        type(oracle_manifest_t), intent(in) :: manifests(:)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: phase_dot(:, :), total_dot(:), memory_dot(:)
        real(dp), intent(in) :: tolerance_dot(:, :), weights_dot(:)
        real(dp), intent(out) :: phase_mean_dot(:), total_mean_dot
        real(dp), intent(out) :: memory_mean_dot, tolerance_mean_dot(:)
        integer, intent(out) :: status

        real(dp) :: phase_values(6, size(manifests))
        real(dp) :: total_values(size(manifests)), memory_values(size(manifests))
        real(dp) :: tolerance_values(4, size(manifests))
        real(dp) :: weight_sum, weight_sum_dot
        real(dp) :: phase_numerator_dot(6), total_numerator_dot
        real(dp) :: memory_numerator_dot, tolerance_numerator_dot(4)
        real(dp) :: phase_mean(6), total_mean, memory_mean, tolerance_mean(4)
        integer :: repetition_count

        phase_mean_dot = 0.0_dp
        total_mean_dot = 0.0_dp
        memory_mean_dot = 0.0_dp
        tolerance_mean_dot = 0.0_dp
        status = 1
        if (size(phase_mean_dot) /= 6 .or. size(tolerance_mean_dot) /= 4) return
        if (.not. validate_records(manifests, weights, repetition_count)) return
        if (.not. valid_jvp_shapes(size(manifests), phase_dot, total_dot, memory_dot, &
            tolerance_dot, weights_dot)) return
        call extract_values(manifests, phase_values, total_values, memory_values, &
            tolerance_values)
        call weighted_means(phase_values, total_values, memory_values, &
            tolerance_values, &
            weights, phase_mean, total_mean, memory_mean, tolerance_mean)
        weight_sum = sum(weights)
        weight_sum_dot = sum(weights_dot)
        phase_numerator_dot = matmul(phase_dot, weights) + &
            matmul(phase_values, weights_dot)
        total_numerator_dot = dot_product(total_dot, weights) + &
            dot_product(total_values, weights_dot)
        memory_numerator_dot = dot_product(memory_dot, weights) + &
            dot_product(memory_values, weights_dot)
        tolerance_numerator_dot = matmul(tolerance_dot, weights) + &
            matmul(tolerance_values, weights_dot)
        phase_mean_dot = (phase_numerator_dot-weight_sum_dot*phase_mean)/weight_sum
        total_mean_dot = (total_numerator_dot-weight_sum_dot*total_mean)/weight_sum
        memory_mean_dot = (memory_numerator_dot-weight_sum_dot*memory_mean)/weight_sum
        tolerance_mean_dot = (tolerance_numerator_dot-weight_sum_dot*tolerance_mean)/ &
            weight_sum
        if (.not. finite_real(phase_mean_dot) .or. &
            .not. ieee_is_finite(total_mean_dot) .or. &
            .not. ieee_is_finite(memory_mean_dot) .or. &
            .not. finite_real(tolerance_mean_dot)) then
            phase_mean_dot = 0.0_dp
            total_mean_dot = 0.0_dp
            memory_mean_dot = 0.0_dp
            tolerance_mean_dot = 0.0_dp
            return
        end if
        status = 0
    end subroutine evaluate_oracle_performance_ledger_jvp

    subroutine evaluate_oracle_performance_ledger_vjp( &
            manifests, weights, phase_bar, total_bar, memory_bar, tolerance_bar, &
            phase_metric_bar, total_metric_bar, memory_metric_bar, tolerance_metric_bar, &
            weights_bar, status)
        type(oracle_manifest_t), intent(in) :: manifests(:)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: phase_bar(:), total_bar, memory_bar, tolerance_bar(:)
        real(dp), intent(out) :: phase_metric_bar(:, :), total_metric_bar(:)
        real(dp), intent(out) :: memory_metric_bar(:), tolerance_metric_bar(:, :)
        real(dp), intent(out) :: weights_bar(:)
        integer, intent(out) :: status

        real(dp) :: phase_values(6, size(manifests))
        real(dp) :: total_values(size(manifests)), memory_values(size(manifests))
        real(dp) :: tolerance_values(4, size(manifests))
        real(dp) :: phase_mean(6), total_mean, memory_mean, tolerance_mean(4)
        real(dp) :: weight_sum
        integer :: repetition_count, run

        phase_metric_bar = 0.0_dp
        total_metric_bar = 0.0_dp
        memory_metric_bar = 0.0_dp
        tolerance_metric_bar = 0.0_dp
        weights_bar = 0.0_dp
        status = 1
        if (size(phase_bar) /= 6 .or. size(tolerance_bar) /= 4 .or. &
            size(phase_metric_bar, 1) /= 6 .or. size(tolerance_metric_bar, 1) /= 4) return
        if (.not. ieee_is_finite(total_bar) .or. .not. ieee_is_finite(memory_bar) .or. &
            .not. finite_real(phase_bar) .or. .not. finite_real(tolerance_bar)) return
        if (.not. validate_records(manifests, weights, repetition_count)) return
        if (size(phase_metric_bar, 2) /= size(manifests) .or. &
            size(total_metric_bar) /= size(manifests) .or. &
            size(memory_metric_bar) /= size(manifests) .or. &
            size(tolerance_metric_bar, 2) /= size(manifests) .or. &
            size(weights_bar) /= size(manifests)) return
        call extract_values(manifests, phase_values, total_values, memory_values, &
            tolerance_values)
        call weighted_means(phase_values, total_values, memory_values, tolerance_values, &
            weights, phase_mean, total_mean, memory_mean, tolerance_mean)
        weight_sum = sum(weights)
        do run = 1, size(manifests)
            phase_metric_bar(:, run) = weights(run)*phase_bar/weight_sum
            total_metric_bar(run) = weights(run)*total_bar/weight_sum
            memory_metric_bar(run) = weights(run)*memory_bar/weight_sum
            tolerance_metric_bar(:, run) = weights(run)*tolerance_bar/weight_sum
            weights_bar(run) = (dot_product(phase_bar, phase_values(:, run)-phase_mean) + &
                total_bar*(total_values(run)-total_mean) + &
                memory_bar*(memory_values(run)-memory_mean) + &
                dot_product(tolerance_bar, tolerance_values(:, run)-tolerance_mean))/weight_sum
        end do
        if (.not. finite_real(phase_metric_bar) .or. .not. finite_real(total_metric_bar) .or. &
            .not. finite_real(memory_metric_bar) .or. .not. finite_real(tolerance_metric_bar) .or. &
            .not. finite_real(weights_bar)) then
            phase_metric_bar = 0.0_dp
            total_metric_bar = 0.0_dp
            memory_metric_bar = 0.0_dp
            tolerance_metric_bar = 0.0_dp
            weights_bar = 0.0_dp
            return
        end if
        status = 0
    end subroutine evaluate_oracle_performance_ledger_vjp

    logical function validate_records(manifests, weights, repetition_count) result(valid)
        type(oracle_manifest_t), intent(in) :: manifests(:)
        real(dp), intent(in) :: weights(:)
        integer, intent(out) :: repetition_count
        integer :: run, status

        valid = .false.
        repetition_count = 0
        if (size(manifests) < 1 .or. size(weights) /= size(manifests)) return
        if (.not. finite_real(weights) .or. any(weights <= 0.0_dp)) return
        do run = 1, size(manifests)
            if (.not. validate_oracle_manifest(manifests(run), status)) return
            if (.not. manifests(run)%success) return
        end do
        repetition_count = manifests(1)%timing%repetition_count
        do run = 2, size(manifests)
            if (manifests(run)%timing%repetition_count /= repetition_count .or. &
                manifests(run)%timing%warmup_count /= manifests(1)%timing%warmup_count) return
            if (.not. same_provenance(manifests(1), manifests(run))) return
        end do
        valid = .true.
    end function validate_records

    logical function same_provenance(reference, candidate) result(same)
        type(oracle_manifest_t), intent(in) :: reference, candidate

        same = trim(reference%case_name) == trim(candidate%case_name) .and. &
            trim(reference%case_revision) == trim(candidate%case_revision) .and. &
            trim(reference%coordinate_system) == trim(candidate%coordinate_system) .and. &
            trim(reference%checksum_algorithm) == trim(candidate%checksum_algorithm) .and. &
            trim(reference%coordinate_checksum) == trim(candidate%coordinate_checksum) .and. &
            trim(reference%sample_checksum) == trim(candidate%sample_checksum) .and. &
            trim(reference%runner_hardware) == trim(candidate%runner_hardware) .and. &
            trim(reference%fortfem_commit) == trim(candidate%fortfem_commit) .and. &
            same_normalization(reference, candidate)
    end function same_provenance

    logical function same_normalization(reference, candidate) result(same)
        type(oracle_manifest_t), intent(in) :: reference, candidate

        same = trim(reference%normalization%normalization_name) == &
            trim(candidate%normalization%normalization_name) .and. &
            trim(reference%normalization%length_unit) == trim(candidate%normalization%length_unit) .and. &
            trim(reference%normalization%time_unit) == trim(candidate%normalization%time_unit) .and. &
            trim(reference%normalization%magnetic_field_unit) == &
            trim(candidate%normalization%magnetic_field_unit) .and. &
            trim(reference%normalization%pressure_unit) == trim(candidate%normalization%pressure_unit) .and. &
            trim(reference%normalization%current_unit) == trim(candidate%normalization%current_unit) .and. &
            all([reference%normalization%length_scale, reference%normalization%time_scale, &
                reference%normalization%magnetic_field_scale, reference%normalization%pressure_scale, &
                reference%normalization%current_scale] == [candidate%normalization%length_scale, &
                candidate%normalization%time_scale, candidate%normalization%magnetic_field_scale, &
                candidate%normalization%pressure_scale, candidate%normalization%current_scale])
    end function same_normalization

    subroutine extract_values(manifests, phase_values, total_values, memory_values, tolerance_values)
        type(oracle_manifest_t), intent(in) :: manifests(:)
        real(dp), intent(out) :: phase_values(:, :), total_values(:), memory_values(:)
        real(dp), intent(out) :: tolerance_values(:, :)
        integer :: run

        do run = 1, size(manifests)
            phase_values(:, run) = [manifests(run)%timing%mesh_seconds, &
                manifests(run)%timing%assembly_seconds, manifests(run)%timing%factorization_seconds, &
                manifests(run)%timing%solve_seconds, manifests(run)%timing%output_seconds, &
                manifests(run)%timing%total_seconds]
            total_values(run) = manifests(run)%timing%total_seconds
            memory_values(run) = real(manifests(run)%timing%peak_memory_bytes, dp)
            tolerance_values(:, run) = [manifests(run)%tolerances%coordinate, &
                manifests(run)%tolerances%absolute, manifests(run)%tolerances%relative, &
                manifests(run)%tolerances%residual]
        end do
    end subroutine extract_values

    subroutine weighted_means(phase_values, total_values, memory_values, tolerance_values, &
            weights, phase_mean, total_mean, memory_mean, tolerance_mean)
        real(dp), intent(in) :: phase_values(:, :), total_values(:), memory_values(:)
        real(dp), intent(in) :: tolerance_values(:, :), weights(:)
        real(dp), intent(out) :: phase_mean(:), total_mean, memory_mean, tolerance_mean(:)
        real(dp) :: weight_sum

        weight_sum = sum(weights)
        phase_mean = matmul(phase_values, weights)/weight_sum
        total_mean = dot_product(total_values, weights)/weight_sum
        memory_mean = dot_product(memory_values, weights)/weight_sum
        tolerance_mean = matmul(tolerance_values, weights)/weight_sum
    end subroutine weighted_means

    logical function valid_jvp_shapes(run_count, phase_dot, total_dot, memory_dot, &
            tolerance_dot, weights_dot) result(valid)
        integer, intent(in) :: run_count
        real(dp), intent(in) :: phase_dot(:, :), total_dot(:), memory_dot(:)
        real(dp), intent(in) :: tolerance_dot(:, :), weights_dot(:)

        valid = size(phase_dot, 1) == 6 .and. size(phase_dot, 2) == run_count .and. &
            size(total_dot) == run_count .and. size(memory_dot) == run_count .and. &
            size(tolerance_dot, 1) == 4 .and. size(tolerance_dot, 2) == run_count .and. &
            size(weights_dot) == run_count .and. finite_real(phase_dot) .and. &
            finite_real(total_dot) .and. finite_real(memory_dot) .and. &
            finite_real(tolerance_dot) .and. finite_real(weights_dot)
    end function valid_jvp_shapes

    logical function finite_real(values) result(valid)
        real(dp), intent(in) :: values(..)

        valid = .true.
        select rank (values)
        rank (1)
            valid = all(ieee_is_finite(values))
        rank (2)
            valid = all(ieee_is_finite(values))
        rank default
            valid = .false.
        end select
    end function finite_real

end module fortfem_oracle_performance_ledger
