module fortfem_linear_response_schema
    !! Versioned, neutral text interchange for linear response records.
    !!
    !! The schema contains metadata, modal labels, and complex blocks only.
    !! It deliberately does not know any application file format, units, or
    !! plasma normalization.  Reads are bounded before allocation so an
    !! untrusted metadata file cannot request an unbounded response matrix.
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortfem_linear_response_interchange, only: &
        initialize_linear_response_interchange, &
        linear_response_interchange_t, &
        validate_linear_response_interchange
    implicit none
    private

    character(len=*), parameter, public :: linear_response_schema_magic = &
        "FORTFEM_LINEAR_RESPONSE_TEXT 1"
    integer(int64), parameter :: max_complex_values = 4_int64*1000_int64*1000_int64

    public :: read_linear_response_interchange
    public :: write_linear_response_interchange

contains

    subroutine write_linear_response_interchange(filename, interchange, status)
        character(len=*), intent(in) :: filename
        type(linear_response_interchange_t), intent(in) :: interchange
        integer, intent(out) :: status

        integer :: unit, ios, validation_status, mode, response
        logical :: valid

        status = 1
        valid = validate_linear_response_interchange(interchange, validation_status)
        if (.not. valid) return
        open(newunit=unit, file=filename, status="replace", action="write", &
            form="formatted", iostat=ios)
        if (ios /= 0) then
            status = 2
            return
        end if
        write(unit, '(A)', iostat=ios) linear_response_schema_magic
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(interchange%schema_version)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(interchange%producer)
        if (ios == 0) write(unit, '(A)', iostat=ios) trim(interchange%provenance)
        if (ios == 0) write(unit, *, iostat=ios) real(interchange%frequency, dp), &
            aimag(interchange%frequency), interchange%normalization_scale
        if (ios == 0) write(unit, *, iostat=ios) interchange%state_count, &
            interchange%mode_count, interchange%response_count
        do mode = 1, interchange%mode_count
            if (ios /= 0) exit
            write(unit, *, iostat=ios) interchange%mode_numbers(:, mode)
            if (ios == 0) write(unit, '(A)', iostat=ios) trim(interchange%mode_labels(mode))
        end do
        do response = 1, interchange%response_count
            if (ios /= 0) exit
            write(unit, '(A)', iostat=ios) trim(interchange%response_labels(response))
        end do
        call write_complex_matrix(unit, interchange%equilibrium_block, ios)
        call write_complex_matrix(unit, interchange%inertia_block, ios)
        call write_complex_matrix(unit, interchange%resistive_block, ios)
        call write_complex_matrix(unit, interchange%vacuum_block, ios)
        call write_complex_matrix(unit, interchange%wall_block, ios)
        call write_complex_matrix(unit, interchange%response_matrix, ios)
        close(unit, iostat=validation_status)
        if (ios /= 0 .or. validation_status /= 0) then
            status = 3
        else
            status = 0
        end if
    end subroutine write_linear_response_interchange

    subroutine read_linear_response_interchange(filename, interchange, status)
        character(len=*), intent(in) :: filename
        type(linear_response_interchange_t), intent(out) :: interchange
        integer, intent(out) :: status

        character(len=256) :: line, schema_version
        character(len=64) :: producer
        character(len=128) :: provenance
        character(len=64), allocatable :: mode_labels(:), response_labels(:)
        complex(dp), allocatable :: equilibrium(:, :), inertia(:, :), resistive(:, :)
        complex(dp), allocatable :: vacuum(:, :), wall(:, :), response_matrix(:, :)
        complex(dp) :: frequency
        real(dp) :: frequency_real, frequency_imag, normalization_scale
        integer, allocatable :: mode_numbers(:, :)
        integer :: unit, ios, state_count, mode_count, response_count
        integer :: mode, response, validation_status
        logical :: valid

        status = 1
        open(newunit=unit, file=filename, status="old", action="read", &
            form="formatted", iostat=ios)
        if (ios /= 0) then
            status = 2
            return
        end if
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0 .or. trim(line) /= linear_response_schema_magic) then
            call close_after_read(unit, status, 3)
            return
        end if
        read(unit, '(A)', iostat=ios) schema_version
        read(unit, '(A)', iostat=ios) producer
        read(unit, '(A)', iostat=ios) provenance
        if (ios /= 0 .or. len_trim(schema_version) > len(schema_version) .or. &
            len_trim(producer) > len(producer) .or. len_trim(provenance) > len(provenance)) then
            call close_after_read(unit, status, 3)
            return
        end if
        read(unit, *, iostat=ios) frequency_real, frequency_imag, normalization_scale
        read(unit, *, iostat=ios) state_count, mode_count, response_count
        if (ios /= 0 .or. trim(schema_version) /= "fortfem-linear-response-1" .or. &
            .not. safe_dimensions(state_count, mode_count, response_count)) then
            call close_after_read(unit, status, 3)
            return
        end if
        allocate(mode_numbers(2, mode_count), mode_labels(mode_count), &
            response_labels(response_count))
        do mode = 1, mode_count
            read(unit, *, iostat=ios) mode_numbers(:, mode)
            if (ios == 0) read(unit, '(A)', iostat=ios) line
            if (ios /= 0 .or. len_trim(line) > len(mode_labels(mode))) then
                call close_after_read(unit, status, 3)
                return
            end if
            mode_labels(mode) = line
        end do
        do response = 1, response_count
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0 .or. len_trim(line) > len(response_labels(response))) then
                call close_after_read(unit, status, 3)
                return
            end if
            response_labels(response) = line
        end do
        allocate(equilibrium(state_count, state_count), inertia(state_count, state_count), &
            resistive(state_count, state_count), vacuum(state_count, state_count), &
            wall(state_count, state_count), response_matrix(response_count, response_count))
        call read_complex_matrix(unit, equilibrium, ios)
        call read_complex_matrix(unit, inertia, ios)
        call read_complex_matrix(unit, resistive, ios)
        call read_complex_matrix(unit, vacuum, ios)
        call read_complex_matrix(unit, wall, ios)
        call read_complex_matrix(unit, response_matrix, ios)
        call close_after_read(unit, status, 0)
        if (ios /= 0 .or. status /= 0) then
            status = 3
            return
        end if
        frequency = cmplx(frequency_real, frequency_imag, dp)
        call initialize_linear_response_interchange( &
            interchange, frequency, mode_numbers, mode_labels, equilibrium, inertia, &
            resistive, vacuum, wall, response_matrix, response_labels, &
            normalization_scale, validation_status)
        if (validation_status /= 0) then
            status = 4
            return
        end if
        interchange%schema_version = trim(schema_version)
        interchange%producer = trim(producer)
        interchange%provenance = trim(provenance)
        valid = validate_linear_response_interchange(interchange, validation_status)
        if (.not. valid) then
            status = 4
            return
        end if
        status = 0
    end subroutine read_linear_response_interchange

    logical function safe_dimensions(state_count, mode_count, response_count) result(valid)
        integer, intent(in) :: state_count, mode_count, response_count
        integer(int64) :: total

        valid = state_count > 0 .and. mode_count > 0 .and. response_count >= 0
        if (.not. valid) return
        total = 5_int64*int(state_count, int64)*int(state_count, int64) + &
            int(response_count, int64)*int(response_count, int64) + &
            2_int64*int(mode_count, int64)
        valid = total <= max_complex_values
    end function safe_dimensions

    subroutine write_complex_matrix(unit, matrix, status)
        integer, intent(in) :: unit
        complex(dp), intent(in) :: matrix(:, :)
        integer, intent(inout) :: status
        integer :: row, column

        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                if (status /= 0) return
                write(unit, *, iostat=status) real(matrix(row, column), dp), &
                    aimag(matrix(row, column))
            end do
        end do
    end subroutine write_complex_matrix

    subroutine read_complex_matrix(unit, matrix, status)
        integer, intent(in) :: unit
        complex(dp), intent(out) :: matrix(:, :)
        integer, intent(inout) :: status
        integer :: row, column
        real(dp) :: real_part, imaginary_part

        do column = 1, size(matrix, 2)
            do row = 1, size(matrix, 1)
                if (status /= 0) return
                read(unit, *, iostat=status) real_part, imaginary_part
                if (status == 0) matrix(row, column) = cmplx(real_part, imaginary_part, dp)
            end do
        end do
    end subroutine read_complex_matrix

    subroutine close_after_read(unit, status, requested_status)
        integer, intent(in) :: unit, requested_status
        integer, intent(out) :: status
        integer :: close_status

        close(unit, iostat=close_status)
        if (requested_status /= 0 .or. close_status /= 0) then
            status = merge(requested_status, 3, requested_status /= 0)
        else
            status = 0
        end if
    end subroutine close_after_read

end module fortfem_linear_response_schema
