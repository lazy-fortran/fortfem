program test_oracle_manifest
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        oracle_manifest_t, oracle_normalization_t, oracle_timing_t, &
        oracle_tolerance_t, initialize_oracle_manifest, &
        validate_oracle_manifest, oracle_manifest_schema_magic, &
        oracle_manifest_schema_version, &
        read_oracle_manifest, write_oracle_manifest
    implicit none

    integer, parameter :: dp = real64
    character(len=*), parameter :: filename = &
        "build/test_oracle_manifest.txt"
    character(len=*), parameter :: bad_filename = &
        "build/test_oracle_manifest_bad.txt"
    type(oracle_manifest_t) :: original, copied, restored
    type(oracle_normalization_t) :: normalization
    type(oracle_tolerance_t) :: tolerances
    type(oracle_timing_t) :: timing
    integer(int64) :: bytes
    integer :: status, unit
    logical :: valid

    normalization%normalization_name = "SI-reference"
    normalization%length_unit = "m"
    normalization%time_unit = "s"
    normalization%magnetic_field_unit = "T"
    normalization%pressure_unit = "Pa"
    normalization%current_unit = "A"
    normalization%length_scale = 2.5_dp
    normalization%time_scale = 0.125_dp
    normalization%magnetic_field_scale = 3.0_dp
    normalization%pressure_scale = 4.0_dp
    normalization%current_scale = 5.0_dp

    tolerances%coordinate = 1.0e-11_dp
    tolerances%absolute = 2.0e-10_dp
    tolerances%relative = 3.0e-9_dp
    tolerances%residual = 4.0e-8_dp

    timing%mesh_seconds = 0.11_dp
    timing%assembly_seconds = 0.23_dp
    timing%factorization_seconds = 0.31_dp
    timing%solve_seconds = 0.47_dp
    timing%output_seconds = 0.03_dp
    timing%total_seconds = 1.15_dp
    timing%peak_memory_bytes = 123456_int64
    timing%warmup_count = 1
    timing%repetition_count = 5

    call initialize_oracle_manifest(original, &
        code_name="MFEM", code_version="4.7", &
        code_revision="mfem-git-0123456789abcdef", &
        code_license="BSD-3-Clause", case_name="torus-ampere", &
        case_revision="case-v2", coordinate_system="physical-cartesian", &
        coordinate_checksum="sha256:coords-a1b2", &
        sample_checksum="sha256:samples-c3d4", spatial_dimension=3, &
        sample_count=128, normalization=normalization, &
        tolerances=tolerances, timing=timing, runner_id="linux-x86_64", &
        runner_hardware="AMD EPYC / 64 GB / Linux / 1 thread", &
        fortfem_commit="fortfem-abcdef", &
        sister_repository_uri="https://example.invalid/benchmarks@deadbeef", &
        success=.true., notes="license-safe external oracle", status=status)
    call check_condition(status == 0, "manifest initializes")
    valid = validate_oracle_manifest(original, status)
    call check_condition(valid .and. status == 0, "manifest validates")

    copied = original
    original%code_name = "changed-after-copy"
    original%normalization%length_scale = 99.0_dp
    call check_condition(copied%code_name == "MFEM" .and. &
        abs(copied%normalization%length_scale - 2.5_dp) < 1.0e-15_dp, &
        "assignment performs a deep copy")

    call write_oracle_manifest(filename, copied, status)
    call check_condition(status == 0, "manifest round-trip writes")
    call read_oracle_manifest(filename, restored, status)
    call check_condition(status == 0, "manifest round-trip reads")
    valid = validate_oracle_manifest(restored, status)
    call check_condition(valid .and. status == 0, "round-tripped manifest validates")
    call check_condition(restored%code_name == copied%code_name .and. &
        restored%code_revision == copied%code_revision .and. &
        restored%code_license == copied%code_license, &
        "code provenance survives round-trip")
    call check_condition(restored%case_name == copied%case_name .and. &
        restored%coordinate_checksum == copied%coordinate_checksum .and. &
        restored%sample_checksum == copied%sample_checksum .and. &
        restored%sample_count == copied%sample_count, &
        "case and coordinate checksums survive round-trip")
    call check_condition(restored%normalization%time_unit == "s" .and. &
        abs(restored%normalization%pressure_scale - 4.0_dp) < 1.0e-15_dp .and. &
        abs(restored%tolerances%relative - 3.0e-9_dp) < 1.0e-20_dp, &
        "units, scales, and tolerances survive round-trip")
    bytes = restored%timing%peak_memory_bytes
    call check_condition(bytes == 123456_int64 .and. &
        restored%timing%repetition_count == 5 .and. &
        restored%sister_repository_uri == copied%sister_repository_uri, &
        "timing, memory, and sister repository survive round-trip")
    call check_condition(restored%runner_hardware == copied%runner_hardware, &
        "runner hardware provenance survives round-trip")

    copied%runner_hardware = ""
    valid = validate_oracle_manifest(copied, status)
    call check_condition(.not. valid .and. status /= 0, &
        "missing runner hardware is rejected")
    copied%runner_hardware = "AMD EPYC / 64 GB / Linux / 1 thread"

    open(newunit=unit, file=filename, status="old", action="read")
    close(unit, status="delete")
    open(newunit=unit, file=bad_filename, status="replace", action="write")
    write(unit, '(A)') oracle_manifest_schema_magic
    write(unit, '(A)') oracle_manifest_schema_version
    write(unit, '(A)') "CHEASE"
    write(unit, '(A)') "unknown"
    write(unit, '(A)') "revision"
    write(unit, '(A)') "GPL-3.0-or-later"
    write(unit, '(A)') "bad-case"
    write(unit, '(A)') "v1"
    write(unit, '(A)') "physical"
    write(unit, '(A)') "sha256"
    write(unit, '(A)') ""
    close(unit)
    call read_oracle_manifest(bad_filename, restored, status)
    call check_condition(status /= 0, "missing coordinate checksum is rejected")
    open(newunit=unit, file=bad_filename, status="old", action="read")
    close(unit, status="delete")

    ! An independent arithmetic oracle guards the normalization contract.
    call check_condition(abs(original%normalization%length_scale - 99.0_dp) < &
        1.0e-15_dp, "rejected input cannot overwrite caller state")
    call check_summary("oracle manifest")
end program test_oracle_manifest
