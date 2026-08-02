program test_oracle_manifest_ladder
    !! Validate the metadata-only external oracle ladder.
    !!
    !! This fixture deliberately does not launch a third-party code or ship
    !! any of its data.  It gives every named target the same synthetic case
    !! and checks that a runner must provide provenance and phase timings
    !! before an external record can enter the comparison pipeline.
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        oracle_manifest_t, oracle_normalization_t, oracle_timing_t, &
        oracle_tolerance_t, initialize_oracle_manifest, &
        validate_oracle_manifest
    implicit none

    integer, parameter :: dp = real64
    character(len=32), parameter :: targets(*) = [character(len=32) :: &
        "CHEASE", "FreeGS", "VMEC/PARVMEC", "GVEC", "DESC", &
        "SPEC/SIESTA", "GPEC", "MARS-F", "GLISS", "STARWALL", &
        "JOREK", "FreeFEM", "MFEM", "FEniCSx"]
    real(dp), parameter :: phase_total = 1.750_dp
    type(oracle_manifest_t) :: manifest, copied
    type(oracle_normalization_t) :: normalization
    type(oracle_tolerance_t) :: tolerances
    type(oracle_timing_t) :: timing
    integer :: index, status
    logical :: valid

    normalization%normalization_name = "SI-reference"
    normalization%length_unit = "m"
    normalization%time_unit = "s"
    normalization%magnetic_field_unit = "T"
    normalization%pressure_unit = "Pa"
    normalization%current_unit = "A"
    tolerances%coordinate = 1.0e-11_dp
    tolerances%absolute = 2.0e-10_dp
    tolerances%relative = 3.0e-9_dp
    tolerances%residual = 4.0e-8_dp
    timing%mesh_seconds = 0.125_dp
    timing%assembly_seconds = 0.250_dp
    timing%factorization_seconds = 0.500_dp
    timing%solve_seconds = 0.750_dp
    timing%output_seconds = 0.125_dp
    timing%total_seconds = phase_total
    timing%peak_memory_bytes = 4096_int64
    timing%warmup_count = 1
    timing%repetition_count = 5

    call check_condition(size(targets) == 14, &
        "external oracle ladder has the expected target count")
    do index = 1, size(targets)
        call initialize_oracle_manifest(manifest, &
            code_name=targets(index), code_version="fixture-v1", &
            code_revision="runner-supplies-revision", &
            code_license="runner-supplies-license", &
            case_name="unit-square-poisson", case_revision="mms-v1", &
            coordinate_system="physical-cartesian", &
            coordinate_checksum="sha256:fixture-coordinates", &
            sample_checksum="sha256:fixture-samples", spatial_dimension=2, &
            sample_count=12, normalization=normalization, &
            tolerances=tolerances, timing=timing, runner_id="fixture-runner", &
            runner_hardware="CI-x86_64 / 32 GB / Linux / 1 thread", &
            fortfem_commit="fixture-commit", &
            sister_repository_uri="https://example.invalid/fortfem-benchmarks", &
            success=.true., notes="metadata-only external target fixture", &
            status=status)
        valid = validate_oracle_manifest(manifest, status)
        call check_condition(valid .and. status == 0 .and. &
            trim(manifest%code_name) == trim(targets(index)), &
            "target manifest validates: "//trim(targets(index)))
        call check_condition(manifest%timing%total_seconds >= &
            manifest%timing%mesh_seconds + manifest%timing%assembly_seconds + &
            manifest%timing%factorization_seconds + &
            manifest%timing%solve_seconds + manifest%timing%output_seconds - &
            1.0e-15_dp .and. manifest%timing%repetition_count > &
            manifest%timing%warmup_count .and. &
            manifest%timing%peak_memory_bytes > 0_int64, &
            "target performance metadata is complete: "//trim(targets(index)))
    end do

    ! Copy semantics must retain the provenance of the last target.
    copied = manifest
    manifest%code_name = "changed-after-copy"
    call check_condition(trim(copied%code_name) == trim(targets(size(targets))), &
        "target manifest provenance survives deep copy")

    ! A runner cannot publish a record without an immutable revision.
    copied%code_revision = ""
    valid = validate_oracle_manifest(copied, status)
    call check_condition(.not. valid .and. status /= 0, &
        "missing external revision is rejected")
    call check_summary("oracle manifest ladder")
end program test_oracle_manifest_ladder
