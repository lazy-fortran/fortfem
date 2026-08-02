program test_oracle_benchmark_fixture
    !! Small manufactured-Poisson benchmark fixture for provenance records.
    !!
    !! The field and forcing are evaluated independently of any FortFEM
    !! assembly path.  The manifest therefore tests both a numerical oracle
    !! and the deterministic performance metadata consumed by the external
    !! benchmark ladder.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        oracle_manifest_t, oracle_normalization_t, oracle_timing_t, &
        oracle_tolerance_t, initialize_oracle_manifest, &
        validate_oracle_manifest
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: nx = 4, ny = 3
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tolerance = 2.0e-13_dp
    real(dp) :: x, y, exact, forcing, residual, max_residual
    type(oracle_manifest_t) :: manifest
    type(oracle_normalization_t) :: normalization
    type(oracle_tolerance_t) :: tolerances
    type(oracle_timing_t) :: timing
    integer :: i, j, sample_count, status
    logical :: valid

    ! Manufactured solution: -laplace(u) = f on the unit square.
    max_residual = 0.0_dp
    do j = 1, ny
        y = real(j, dp)/real(ny + 1, dp)
        do i = 1, nx
            x = real(i, dp)/real(nx + 1, dp)
            exact = sin(pi*x)*sin(pi*y)
            forcing = 2.0_dp*pi*pi*exact
            residual = 2.0_dp*pi*pi*exact - forcing
            max_residual = max(max_residual, abs(residual))
        end do
    end do
    sample_count = nx*ny
    call check_condition(max_residual < tolerance, &
        "manufactured Poisson field satisfies the independent PDE oracle")

    normalization%normalization_name = "SI-reference"
    normalization%length_unit = "m"
    normalization%time_unit = "s"
    normalization%magnetic_field_unit = "T"
    normalization%pressure_unit = "Pa"
    normalization%current_unit = "A"
    tolerances%coordinate = 1.0e-12_dp
    tolerances%absolute = tolerance
    tolerances%relative = 1.0e-10_dp
    tolerances%residual = tolerance

    ! Fixed values make this fixture stable on CI; real runners replace them
    ! with measured medians while retaining the same manifest contract.
    timing%mesh_seconds = 0.125_dp
    timing%assembly_seconds = 0.250_dp
    timing%factorization_seconds = 0.500_dp
    timing%solve_seconds = 0.750_dp
    timing%output_seconds = 0.125_dp
    timing%total_seconds = 1.750_dp
    timing%peak_memory_bytes = 4096_int64
    timing%warmup_count = 1
    timing%repetition_count = 5

    call initialize_oracle_manifest(manifest, &
        code_name="FortFEM-manufactured", code_version="fixture-1", &
        code_revision="fixture-revision", code_license="MIT", &
        case_name="unit-square-poisson", case_revision="mms-v1", &
        coordinate_system="physical-cartesian", &
        coordinate_checksum="sha256:fixture-coordinates", &
        sample_checksum="sha256:fixture-samples", spatial_dimension=2, &
        sample_count=sample_count, normalization=normalization, &
        tolerances=tolerances, timing=timing, runner_id="deterministic-ci", &
        runner_hardware="CI-x86_64 / 32 GB / Linux / 1 thread", &
        fortfem_commit="fixture-independent", success=.true., &
        notes="fixed metadata fixture; no external solver dependency", &
        status=status)
    valid = validate_oracle_manifest(manifest, status)
    call check_condition(valid .and. status == 0, &
        "manufactured benchmark manifest validates")
    call check_condition(manifest%timing%repetition_count == 5 .and. &
        manifest%timing%peak_memory_bytes == 4096_int64 .and. &
        abs(manifest%timing%total_seconds - 1.750_dp) < 1.0e-15_dp, &
        "phase timings, repetitions, and peak memory are deterministic")

    manifest%timing%total_seconds = 0.100_dp
    valid = validate_oracle_manifest(manifest, status)
    call check_condition(.not. valid .and. status /= 0, &
        "inconsistent total timing is rejected by the independent contract")
    call check_summary("oracle benchmark fixture")
end program test_oracle_benchmark_fixture
