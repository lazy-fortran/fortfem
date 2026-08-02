module fortfem_oracle_manifest
    !! License-safe provenance contract for external benchmark oracles.
    !!
    !! This module stores metadata only.  It neither reads an external code's
    !! input/output format nor embeds its source or physics.  A caller-owned
    !! adapter can emit one manifest next to sampled fields in the separate
    !! benchmark-data repository.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    implicit none
    private

    character(len=*), parameter, public :: oracle_manifest_schema_magic = &
        "FORTFEM_ORACLE_MANIFEST_TEXT 1"
    character(len=*), parameter, public :: oracle_manifest_schema_version = &
        "fortfem-oracle-manifest-1"

    type, public :: oracle_normalization_t
        character(len=64) :: normalization_name = ""
        character(len=32) :: length_unit = "1"
        character(len=32) :: time_unit = "1"
        character(len=32) :: magnetic_field_unit = "1"
        character(len=32) :: pressure_unit = "1"
        character(len=32) :: current_unit = "1"
        real(dp) :: length_scale = 1.0_dp
        real(dp) :: time_scale = 1.0_dp
        real(dp) :: magnetic_field_scale = 1.0_dp
        real(dp) :: pressure_scale = 1.0_dp
        real(dp) :: current_scale = 1.0_dp
    end type oracle_normalization_t

    type, public :: oracle_tolerance_t
        real(dp) :: coordinate = 0.0_dp
        real(dp) :: absolute = 0.0_dp
        real(dp) :: relative = 0.0_dp
        real(dp) :: residual = 0.0_dp
    end type oracle_tolerance_t

    type, public :: oracle_timing_t
        real(dp) :: mesh_seconds = 0.0_dp
        real(dp) :: assembly_seconds = 0.0_dp
        real(dp) :: factorization_seconds = 0.0_dp
        real(dp) :: solve_seconds = 0.0_dp
        real(dp) :: output_seconds = 0.0_dp
        real(dp) :: total_seconds = 0.0_dp
        integer(int64) :: peak_memory_bytes = 0_int64
        integer :: warmup_count = 0
        integer :: repetition_count = 1
    end type oracle_timing_t

    type, public :: oracle_manifest_t
        character(len=32) :: schema_version = &
            oracle_manifest_schema_version
        character(len=64) :: code_name = ""
        character(len=64) :: code_version = ""
        character(len=128) :: code_revision = ""
        character(len=128) :: code_license = ""
        character(len=128) :: case_name = ""
        character(len=64) :: case_revision = ""
        character(len=64) :: coordinate_system = ""
        character(len=32) :: checksum_algorithm = ""
        character(len=128) :: coordinate_checksum = ""
        character(len=128) :: sample_checksum = ""
        character(len=64) :: runner_id = ""
        character(len=64) :: fortfem_commit = ""
        character(len=256) :: sister_repository_uri = ""
        character(len=256) :: notes = ""
        integer :: spatial_dimension = 0
        integer :: sample_count = 0
        logical :: success = .false.
        type(oracle_normalization_t) :: normalization
        type(oracle_tolerance_t) :: tolerances
        type(oracle_timing_t) :: timing
    end type oracle_manifest_t

    public :: initialize_oracle_manifest
    public :: validate_oracle_manifest
    public :: read_oracle_manifest
    public :: write_oracle_manifest

contains

    subroutine initialize_oracle_manifest(manifest, code_name, code_version, &
            code_revision, code_license, case_name, case_revision, &
            coordinate_system, coordinate_checksum, sample_checksum, &
            spatial_dimension, sample_count, normalization, tolerances, &
            timing, runner_id, fortfem_commit, sister_repository_uri, &
            success, notes, status)
        type(oracle_manifest_t), intent(out) :: manifest
        character(len=*), intent(in) :: code_name, code_version
        character(len=*), intent(in) :: code_revision, code_license
        character(len=*), intent(in) :: case_name, case_revision
        character(len=*), intent(in) :: coordinate_system
        character(len=*), intent(in) :: coordinate_checksum, sample_checksum
        integer, intent(in) :: spatial_dimension, sample_count
        type(oracle_normalization_t), intent(in) :: normalization
        type(oracle_tolerance_t), intent(in) :: tolerances
        type(oracle_timing_t), intent(in) :: timing
        character(len=*), intent(in) :: runner_id, fortfem_commit
        character(len=*), intent(in), optional :: sister_repository_uri
        logical, intent(in), optional :: success
        character(len=*), intent(in), optional :: notes
        integer, intent(out) :: status

        manifest = oracle_manifest_t()
        manifest%code_name = code_name
        manifest%code_version = code_version
        manifest%code_revision = code_revision
        manifest%code_license = code_license
        manifest%case_name = case_name
        manifest%case_revision = case_revision
        manifest%coordinate_system = coordinate_system
        manifest%checksum_algorithm = "sha256"
        manifest%coordinate_checksum = coordinate_checksum
        manifest%sample_checksum = sample_checksum
        manifest%spatial_dimension = spatial_dimension
        manifest%sample_count = sample_count
        manifest%normalization = normalization
        manifest%tolerances = tolerances
        manifest%timing = timing
        manifest%runner_id = runner_id
        manifest%fortfem_commit = fortfem_commit
        if (present(sister_repository_uri)) then
            manifest%sister_repository_uri = sister_repository_uri
        end if
        if (present(success)) manifest%success = success
        if (present(notes)) manifest%notes = notes
        if (.not. validate_oracle_manifest(manifest, status)) then
            manifest = oracle_manifest_t()
            return
        end if
    end subroutine initialize_oracle_manifest

    logical function validate_oracle_manifest(manifest, status) result(valid)
        type(oracle_manifest_t), intent(in) :: manifest
        integer, intent(out) :: status
        real(dp) :: scales(5), tolerance_values(4), timing_values(6)

        valid = .false.
        status = 1
        if (manifest%schema_version /= oracle_manifest_schema_version) return
        if (len_trim(manifest%code_name) == 0 .or. &
            len_trim(manifest%code_version) == 0 .or. &
            len_trim(manifest%code_revision) == 0 .or. &
            len_trim(manifest%code_license) == 0 .or. &
            len_trim(manifest%case_name) == 0 .or. &
            len_trim(manifest%case_revision) == 0 .or. &
            len_trim(manifest%coordinate_system) == 0 .or. &
            len_trim(manifest%checksum_algorithm) == 0 .or. &
            len_trim(manifest%coordinate_checksum) == 0 .or. &
            len_trim(manifest%sample_checksum) == 0 .or. &
            len_trim(manifest%runner_id) == 0) return
        if (manifest%spatial_dimension < 1 .or. &
            manifest%spatial_dimension > 3 .or. manifest%sample_count < 1) return
        scales = [manifest%normalization%length_scale, &
            manifest%normalization%time_scale, &
            manifest%normalization%magnetic_field_scale, &
            manifest%normalization%pressure_scale, &
            manifest%normalization%current_scale]
        if (len_trim(manifest%normalization%normalization_name) == 0 .or. &
            len_trim(manifest%normalization%length_unit) == 0 .or. &
            len_trim(manifest%normalization%time_unit) == 0 .or. &
            len_trim(manifest%normalization%magnetic_field_unit) == 0 .or. &
            len_trim(manifest%normalization%pressure_unit) == 0 .or. &
            len_trim(manifest%normalization%current_unit) == 0 .or. &
            .not. all(ieee_is_finite(scales)) .or. any(scales <= 0.0_dp)) return
        tolerance_values = [manifest%tolerances%coordinate, &
            manifest%tolerances%absolute, manifest%tolerances%relative, &
            manifest%tolerances%residual]
        if (.not. all(ieee_is_finite(tolerance_values)) .or. &
            any(tolerance_values < 0.0_dp)) return
        timing_values = [manifest%timing%mesh_seconds, &
            manifest%timing%assembly_seconds, &
            manifest%timing%factorization_seconds, &
            manifest%timing%solve_seconds, manifest%timing%output_seconds, &
            manifest%timing%total_seconds]
        if (.not. all(ieee_is_finite(timing_values)) .or. &
            any(timing_values < 0.0_dp)) return
        if (manifest%timing%total_seconds < maxval(timing_values(:5))) return
        if (manifest%timing%peak_memory_bytes < 0_int64 .or. &
            manifest%timing%warmup_count < 0 .or. &
            manifest%timing%repetition_count < 1) return
        valid = .true.
        status = 0
    end function validate_oracle_manifest

    subroutine write_oracle_manifest(filename, manifest, status)
        character(len=*), intent(in) :: filename
        type(oracle_manifest_t), intent(in) :: manifest
        integer, intent(out) :: status
        integer :: unit, ios, validation_status

        status = 1
        if (.not. validate_oracle_manifest(manifest, validation_status)) return
        open(newunit=unit, file=filename, status="replace", action="write", &
            form="formatted", iostat=ios)
        if (ios /= 0) then
            status = 2
            return
        end if
        write(unit, '(A)', iostat=ios) oracle_manifest_schema_magic
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%schema_version)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%code_name)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%code_version)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%code_revision)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%code_license)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%case_name)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%case_revision)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%coordinate_system)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%checksum_algorithm)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%coordinate_checksum)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%sample_checksum)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%runner_id)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%fortfem_commit)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%sister_repository_uri)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(manifest%notes)
        if (ios == 0) write(unit, *, iostat=ios) manifest%spatial_dimension, &
            manifest%sample_count, merge(1, 0, manifest%success)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%normalization%normalization_name)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%normalization%length_unit)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%normalization%time_unit)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%normalization%magnetic_field_unit)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%normalization%pressure_unit)
        if (ios == 0) write(unit, '(A)', iostat=ios) &
            trim(manifest%normalization%current_unit)
        if (ios == 0) write(unit, *, iostat=ios) &
            manifest%normalization%length_scale, &
            manifest%normalization%time_scale, &
            manifest%normalization%magnetic_field_scale, &
            manifest%normalization%pressure_scale, &
            manifest%normalization%current_scale
        if (ios == 0) write(unit, *, iostat=ios) manifest%tolerances%coordinate, &
            manifest%tolerances%absolute, manifest%tolerances%relative, &
            manifest%tolerances%residual
        if (ios == 0) write(unit, *, iostat=ios) manifest%timing%mesh_seconds, &
            manifest%timing%assembly_seconds, manifest%timing%factorization_seconds, &
            manifest%timing%solve_seconds, manifest%timing%output_seconds, &
            manifest%timing%total_seconds, manifest%timing%peak_memory_bytes, &
            manifest%timing%warmup_count, manifest%timing%repetition_count
        close(unit, iostat=validation_status)
        if (ios /= 0 .or. validation_status /= 0) then
            status = 3
        else
            status = 0
        end if
    end subroutine write_oracle_manifest

    subroutine read_oracle_manifest(filename, manifest, status)
        character(len=*), intent(in) :: filename
        type(oracle_manifest_t), intent(out) :: manifest
        integer, intent(out) :: status
        character(len=512) :: line
        integer :: unit, ios, close_status, dimensions(3)

        manifest = oracle_manifest_t()
        status = 1
        open(newunit=unit, file=filename, status="old", action="read", &
            form="formatted", iostat=ios)
        if (ios /= 0) then
            status = 2
            return
        end if
        call read_line(unit, line, ios)
        if (ios /= 0 .or. trim(line) /= oracle_manifest_schema_magic) then
            call finish_read(unit, status, 3)
            return
        end if
        call read_line_into(unit, manifest%schema_version, ios)
        call read_line_into(unit, manifest%code_name, ios)
        call read_line_into(unit, manifest%code_version, ios)
        call read_line_into(unit, manifest%code_revision, ios)
        call read_line_into(unit, manifest%code_license, ios)
        call read_line_into(unit, manifest%case_name, ios)
        call read_line_into(unit, manifest%case_revision, ios)
        call read_line_into(unit, manifest%coordinate_system, ios)
        call read_line_into(unit, manifest%checksum_algorithm, ios)
        call read_line_into(unit, manifest%coordinate_checksum, ios)
        call read_line_into(unit, manifest%sample_checksum, ios)
        call read_line_into(unit, manifest%runner_id, ios)
        call read_line_into(unit, manifest%fortfem_commit, ios)
        call read_line_into(unit, manifest%sister_repository_uri, ios)
        call read_line_into(unit, manifest%notes, ios)
        if (ios == 0) read(unit, *, iostat=ios) dimensions
        if (ios == 0) then
            manifest%spatial_dimension = dimensions(1)
            manifest%sample_count = dimensions(2)
            if (dimensions(3) == 0) then
                manifest%success = .false.
            else if (dimensions(3) == 1) then
                manifest%success = .true.
            else
                ios = 1
            end if
        end if
        call read_line_into(unit, manifest%normalization%normalization_name, ios)
        call read_line_into(unit, manifest%normalization%length_unit, ios)
        call read_line_into(unit, manifest%normalization%time_unit, ios)
        call read_line_into(unit, manifest%normalization%magnetic_field_unit, ios)
        call read_line_into(unit, manifest%normalization%pressure_unit, ios)
        call read_line_into(unit, manifest%normalization%current_unit, ios)
        if (ios == 0) read(unit, *, iostat=ios) &
            manifest%normalization%length_scale, &
            manifest%normalization%time_scale, &
            manifest%normalization%magnetic_field_scale, &
            manifest%normalization%pressure_scale, &
            manifest%normalization%current_scale
        if (ios == 0) read(unit, *, iostat=ios) manifest%tolerances%coordinate, &
            manifest%tolerances%absolute, manifest%tolerances%relative, &
            manifest%tolerances%residual
        if (ios == 0) read(unit, *, iostat=ios) manifest%timing%mesh_seconds, &
            manifest%timing%assembly_seconds, manifest%timing%factorization_seconds, &
            manifest%timing%solve_seconds, manifest%timing%output_seconds, &
            manifest%timing%total_seconds, manifest%timing%peak_memory_bytes, &
            manifest%timing%warmup_count, manifest%timing%repetition_count
        close(unit, iostat=close_status)
        if (ios /= 0 .or. close_status /= 0) then
            status = 3
            return
        end if
        if (.not. validate_oracle_manifest(manifest, status)) then
            status = 4
            return
        end if
        status = 0
    end subroutine read_oracle_manifest

    subroutine read_line(unit, line, status)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: line
        integer, intent(out) :: status

        read(unit, '(A)', iostat=status) line
    end subroutine read_line

    subroutine read_line_into(unit, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: value
        integer, intent(inout) :: status
        character(len=512) :: line

        if (status /= 0) return
        call read_line(unit, line, status)
        if (status /= 0) return
        if (len_trim(line) > len(value)) then
            status = 1
            return
        end if
        value = trim(line)
    end subroutine read_line_into

    subroutine finish_read(unit, status, requested_status)
        integer, intent(in) :: unit, requested_status
        integer, intent(out) :: status
        integer :: close_status

        close(unit, iostat=close_status)
        status = merge(requested_status, 3, requested_status /= 0 .or. &
            close_status /= 0)
    end subroutine finish_read

end module fortfem_oracle_manifest
