program test_oracle_performance_ledger
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use check, only: check_condition, check_summary
    use fortfem_oracle_manifest, only: &
        oracle_manifest_t, oracle_normalization_t, oracle_tolerance_t, &
        oracle_timing_t, initialize_oracle_manifest
    use fortfem_oracle_performance_ledger, only: &
        evaluate_oracle_performance_ledger, &
        evaluate_oracle_performance_ledger_jvp, &
        evaluate_oracle_performance_ledger_vjp
    implicit none

    integer, parameter :: dp = real64, run_count = 2
    real(dp), parameter :: epsilon = 1.0e-6_dp
    type(oracle_manifest_t) :: manifests(run_count), manifests_plus(run_count)
    type(oracle_normalization_t) :: normalization
    type(oracle_tolerance_t) :: tolerances
    type(oracle_timing_t) :: timing
    real(dp) :: weights(run_count), weights_dot(run_count), weights_plus(run_count)
    real(dp) :: phase_dot(6, run_count), memory_dot(run_count)
    real(dp) :: total_dot(run_count), tolerance_dot(4, run_count)
    real(dp) :: phase_mean(6), phase_mean_dot(6), phase_mean_plus(6)
    real(dp) :: total_mean, total_mean_dot, total_mean_plus
    real(dp) :: memory_mean, memory_mean_dot, memory_mean_plus
    real(dp) :: tolerance_mean(4), tolerance_mean_dot(4), tolerance_mean_plus(4)
    real(dp) :: phase_bar(6), total_bar, memory_bar, tolerance_bar(4)
    real(dp) :: phase_metric_bar(6, run_count), total_metric_bar(run_count)
    real(dp) :: memory_metric_bar(run_count), tolerance_metric_bar(4, run_count)
    real(dp) :: weights_bar(run_count), lhs, rhs
    real(dp) :: phase_raw(6, run_count), total_raw(run_count)
    real(dp) :: memory_raw(run_count), tolerance_raw(4, run_count)
    real(dp) :: weight_sum, expected_phase(6), expected_total, expected_memory
    real(dp) :: expected_tolerance(4)
    real(dp) :: expected_memory_dot
    integer :: repetition_count, status, run, phase, tolerance_index

    normalization%normalization_name = "SI-reference"
    normalization%length_unit = "m"
    normalization%time_unit = "s"
    normalization%magnetic_field_unit = "T"
    normalization%pressure_unit = "Pa"
    normalization%current_unit = "A"
    tolerances%coordinate = 1.0e-10_dp
    tolerances%absolute = 2.0e-8_dp
    tolerances%relative = 3.0e-7_dp
    tolerances%residual = 4.0e-9_dp
    timing%warmup_count = 1
    timing%repetition_count = 4
    timing%mesh_seconds = 0.2_dp
    timing%assembly_seconds = 0.4_dp
    timing%factorization_seconds = 0.7_dp
    timing%solve_seconds = 1.1_dp
    timing%output_seconds = 0.3_dp
    timing%total_seconds = 2.7_dp
    timing%peak_memory_bytes = 4096_int64
    call initialize_manifest(manifests(1), "FEM", "1", timing, tolerances, status)
    call check_condition(status == 0, "first oracle performance provenance initializes")
    timing%mesh_seconds = 0.5_dp
    timing%assembly_seconds = 0.3_dp
    timing%factorization_seconds = 0.8_dp
    timing%solve_seconds = 1.4_dp
    timing%output_seconds = 0.4_dp
    timing%total_seconds = 3.4_dp
    timing%peak_memory_bytes = 8192_int64
    tolerances%absolute = 4.0e-8_dp
    tolerances%relative = 5.0e-7_dp
    call initialize_manifest(manifests(2), "BEM", "2", timing, tolerances, status)
    call check_condition(status == 0, "second oracle performance provenance initializes")
    weights = [1.0_dp, 3.0_dp]
    weights_dot = [0.1_dp, -0.2_dp]
    phase_raw(:, 1) = [0.2_dp, 0.4_dp, 0.7_dp, 1.1_dp, 0.3_dp, 2.7_dp]
    phase_raw(:, 2) = [0.5_dp, 0.3_dp, 0.8_dp, 1.4_dp, 0.4_dp, 3.4_dp]
    total_raw = phase_raw(6, :)
    memory_raw = [4096.0_dp, 8192.0_dp]
    tolerance_raw(:, 1) = [1.0e-10_dp, 2.0e-8_dp, 3.0e-7_dp, 4.0e-9_dp]
    tolerance_raw(:, 2) = [1.0e-10_dp, 4.0e-8_dp, 5.0e-7_dp, 4.0e-9_dp]
    phase_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, 0.04_dp, -0.01_dp, 0.05_dp, &
        -0.02_dp, 0.03_dp, 0.01_dp, -0.04_dp, 0.02_dp, -0.03_dp], &
        shape(phase_dot))
    total_dot = phase_dot(6, :)
    memory_dot = [1.0e6_dp, -1.0e6_dp]
    tolerance_dot = reshape([ &
        0.1e-10_dp, 0.2e-8_dp, -0.3e-7_dp, 0.4e-9_dp, &
        -0.2e-10_dp, 0.3e-8_dp, 0.1e-7_dp, -0.2e-9_dp], shape(tolerance_dot))

    call evaluate_oracle_performance_ledger( &
        manifests, weights, phase_mean, total_mean, memory_mean, tolerance_mean, &
        repetition_count, status)
    weight_sum = sum(weights)
    expected_phase = matmul(phase_raw, weights)/weight_sum
    expected_total = sum(weights*total_raw)/weight_sum
    expected_memory = sum(weights*memory_raw)/weight_sum
    expected_tolerance = matmul(tolerance_raw, weights)/weight_sum
    call check_condition(status == 0 .and. repetition_count == 4 .and. &
        maxval(abs(phase_mean-expected_phase)) < 1.0e-14_dp .and. &
        abs(total_mean-expected_total) < 1.0e-14_dp .and. &
        abs(memory_mean-expected_memory) < 1.0e-12_dp .and. &
        maxval(abs(tolerance_mean-expected_tolerance)) < 1.0e-14_dp, &
        "weighted phase, memory, total, and tolerance metrics match oracle")

    call evaluate_oracle_performance_ledger_jvp( &
        manifests, weights, phase_dot, total_dot, memory_dot, tolerance_dot, weights_dot, &
        phase_mean_dot, total_mean_dot, memory_mean_dot, tolerance_mean_dot, status)
    weights_plus = weights + epsilon*weights_dot
    manifests_plus = manifests
    do run = 1, run_count
        manifests_plus(run)%timing%mesh_seconds = manifests(run)%timing%mesh_seconds + &
            epsilon*phase_dot(1, run)
        manifests_plus(run)%timing%assembly_seconds = manifests(run)%timing%assembly_seconds + &
            epsilon*phase_dot(2, run)
        manifests_plus(run)%timing%factorization_seconds = manifests(run)%timing%factorization_seconds + &
            epsilon*phase_dot(3, run)
        manifests_plus(run)%timing%solve_seconds = manifests(run)%timing%solve_seconds + &
            epsilon*phase_dot(4, run)
        manifests_plus(run)%timing%output_seconds = manifests(run)%timing%output_seconds + &
            epsilon*phase_dot(5, run)
        manifests_plus(run)%timing%total_seconds = manifests(run)%timing%total_seconds + &
            epsilon*total_dot(run)
        manifests_plus(run)%timing%peak_memory_bytes = int( &
            nint(real(manifests(run)%timing%peak_memory_bytes, dp) + epsilon*memory_dot(run)), int64)
        manifests_plus(run)%tolerances%coordinate = manifests(run)%tolerances%coordinate + &
            epsilon*tolerance_dot(1, run)
        manifests_plus(run)%tolerances%absolute = manifests(run)%tolerances%absolute + &
            epsilon*tolerance_dot(2, run)
        manifests_plus(run)%tolerances%relative = manifests(run)%tolerances%relative + &
            epsilon*tolerance_dot(3, run)
        manifests_plus(run)%tolerances%residual = manifests(run)%tolerances%residual + &
            epsilon*tolerance_dot(4, run)
    end do
    call evaluate_oracle_performance_ledger( &
        manifests_plus, weights_plus, phase_mean_plus, total_mean_plus, memory_mean_plus, &
        tolerance_mean_plus, repetition_count, status)
    expected_memory_dot = (dot_product(memory_dot, weights) + &
        dot_product(memory_raw, weights_dot) - sum(weights_dot)*memory_mean)/weight_sum
    call check_condition(status == 0 .and. &
        maxval(abs(phase_mean_dot-(phase_mean_plus-phase_mean)/epsilon)) < 1.0e-6_dp .and. &
        abs(total_mean_dot-(total_mean_plus-total_mean)/epsilon) < 1.0e-6_dp .and. &
        abs(memory_mean_dot-expected_memory_dot) < 1.0e-8_dp .and. &
        maxval(abs(tolerance_mean_dot-(tolerance_mean_plus-tolerance_mean)/epsilon)) < 1.0e-6_dp, &
        "weighted metrics JVP matches deterministic perturbation oracle")

    phase_bar = [0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.2_dp]
    total_bar = 0.6_dp
    memory_bar = -0.7_dp
    tolerance_bar = [0.8_dp, -0.9_dp, 0.4_dp, -0.3_dp]
    call evaluate_oracle_performance_ledger_vjp( &
        manifests, weights, phase_bar, total_bar, memory_bar, tolerance_bar, &
        phase_metric_bar, total_metric_bar, memory_metric_bar, tolerance_metric_bar, &
        weights_bar, status)
    lhs = dot_product(phase_bar, phase_mean_dot) + total_bar*total_mean_dot + &
        memory_bar*memory_mean_dot + dot_product(tolerance_bar, tolerance_mean_dot)
    rhs = sum(phase_metric_bar*phase_dot) + dot_product(total_metric_bar, total_dot) + &
        dot_product(memory_metric_bar, memory_dot) + sum(tolerance_metric_bar*tolerance_dot) + &
        dot_product(weights_bar, weights_dot)
    call check_condition(status == 0 .and. abs(lhs-rhs) < 2.0e-7_dp, &
        "weighted metrics VJP satisfies the real ledger adjoint oracle")

    manifests(2)%case_revision = "different-case"
    call evaluate_oracle_performance_ledger( &
        manifests, weights, phase_mean, total_mean, memory_mean, tolerance_mean, &
        repetition_count, status)
    call check_condition(status /= 0, "ledger rejects mismatched case provenance")
    manifests(2)%case_revision = manifests(1)%case_revision
    manifests(2)%timing%repetition_count = 5
    call evaluate_oracle_performance_ledger( &
        manifests, weights, phase_mean, total_mean, memory_mean, tolerance_mean, &
        repetition_count, status)
    call check_condition(status /= 0, "ledger rejects mismatched repetition counts")
    manifests(2)%timing%repetition_count = manifests(1)%timing%repetition_count
    weights(1) = 0.0_dp
    call evaluate_oracle_performance_ledger( &
        manifests, weights, phase_mean, total_mean, memory_mean, tolerance_mean, &
        repetition_count, status)
    call check_condition(status /= 0, "ledger rejects nonpositive benchmark weights")

    call check_summary("oracle performance ledger")

contains

    subroutine initialize_manifest(manifest, code_name, code_version, timing, tolerances, status)
        type(oracle_manifest_t), intent(out) :: manifest
        character(len=*), intent(in) :: code_name, code_version
        type(oracle_timing_t), intent(in) :: timing
        type(oracle_tolerance_t), intent(in) :: tolerances
        integer, intent(out) :: status

        call initialize_oracle_manifest(manifest, code_name, code_version, &
            "revision-"//trim(code_name), "BSD-3-Clause", "poisson-mms", "mms-v1", &
            "cartesian", "sha256:coordinates", "sha256:samples", 2, 12, &
            normalization, tolerances, timing, "runner-ci", &
            "CI-x86_64 / 32 GB / Linux / 1 thread", "fortfem-test", &
            success=.true., status=status)
    end subroutine initialize_manifest

end program test_oracle_performance_ledger
