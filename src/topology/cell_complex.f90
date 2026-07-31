module fortfem_cell_complex
    !! Oriented finite cell-complex metadata and homology diagnostics.
    !!
    !! The boundary matrices use the chain convention
    !!
    !!     C_3 -- boundary_3 --> C_2 -- boundary_2 --> C_1
    !!       -- boundary_1 --> C_0.
    !!
    !! Columns carry oriented cells and rows carry their oriented boundary
    !! incidences.  This module deliberately owns no mesh or application data:
    !! it is the small topological contract consumed by compatible spaces,
    !! gauges, cuts, and interface graphs.
    use iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortfem_cell_identification, only: cell_identification_t, &
        identify_boundary_matrix
    implicit none
    private

    type, public :: cell_complex_t
        private
        integer :: vertex_count = 0
        integer, allocatable :: boundary_1(:, :)
        integer, allocatable :: boundary_2(:, :)
        integer, allocatable :: boundary_3(:, :)
    end type cell_complex_t

    public :: initialize_cell_complex
    public :: validate_cell_complex
    public :: cell_complex_euler_characteristic
    public :: cell_complex_betti_numbers
    public :: cell_complex_cycle_basis
    public :: quotient_cell_complex

contains

    subroutine initialize_cell_complex( &
            complex, vertex_count, boundary_1, boundary_2, boundary_3, status)
        type(cell_complex_t), intent(inout) :: complex
        integer, intent(in) :: vertex_count
        integer, intent(in) :: boundary_1(:, :)
        integer, intent(in), optional :: boundary_2(:, :), boundary_3(:, :)
        integer, intent(out) :: status
        integer :: edge_count, face_count

        call clear_cell_complex(complex)
        status = 1
        if (vertex_count < 0) return
        if (size(boundary_1, 1) /= vertex_count) then
            status = 2
            return
        end if
        edge_count = size(boundary_1, 2)
        face_count = 0
        if (present(boundary_2)) then
            if (size(boundary_2, 1) /= edge_count) then
                status = 3
                return
            end if
            face_count = size(boundary_2, 2)
        end if
        if (present(boundary_3)) then
            if (.not. present(boundary_2)) then
                status = 4
                return
            end if
            if (size(boundary_3, 1) /= face_count) then
                status = 5
                return
            end if
        end if

        complex%vertex_count = vertex_count
        allocate(complex%boundary_1(vertex_count, edge_count))
        complex%boundary_1 = boundary_1
        if (present(boundary_2)) then
            allocate(complex%boundary_2(edge_count, face_count))
            complex%boundary_2 = boundary_2
        else
            allocate(complex%boundary_2(edge_count, 0))
        end if
        if (present(boundary_3)) then
            allocate(complex%boundary_3(face_count, size(boundary_3, 2)))
            complex%boundary_3 = boundary_3
        else
            allocate(complex%boundary_3(face_count, 0))
        end if
        status = 0
    end subroutine initialize_cell_complex

    subroutine validate_cell_complex(complex, status)
        type(cell_complex_t), intent(in) :: complex
        integer, intent(out) :: status
        integer(int64), allocatable :: composition(:, :)

        status = 1
        if (complex%vertex_count < 0) return
        if (.not. allocated(complex%boundary_1) .or. &
            .not. allocated(complex%boundary_2) .or. &
            .not. allocated(complex%boundary_3)) return
        if (size(complex%boundary_1, 1) /= complex%vertex_count) then
            status = 2
            return
        end if
        if (size(complex%boundary_2, 1) /= size(complex%boundary_1, 2)) then
            status = 3
            return
        end if
        if (size(complex%boundary_3, 1) /= size(complex%boundary_2, 2)) then
            status = 4
            return
        end if

        if (size(complex%boundary_2, 2) > 0 .and. &
            size(complex%boundary_1, 2) > 0) then
            composition = matmul( &
                int(complex%boundary_1, int64), &
                int(complex%boundary_2, int64))
            if (any(composition /= 0_int64)) then
                status = 5
                return
            end if
        end if
        if (size(complex%boundary_3, 2) > 0 .and. &
            size(complex%boundary_2, 2) > 0) then
            composition = matmul( &
                int(complex%boundary_2, int64), &
                int(complex%boundary_3, int64))
            if (any(composition /= 0_int64)) then
                status = 6
                return
            end if
        end if
        status = 0
    end subroutine validate_cell_complex

    subroutine cell_complex_euler_characteristic(complex, euler, status)
        type(cell_complex_t), intent(in) :: complex
        integer, intent(out) :: euler
        integer, intent(out) :: status

        euler = 0
        call validate_cell_complex(complex, status)
        if (status /= 0) return
        euler = complex%vertex_count - size(complex%boundary_1, 2) + &
            size(complex%boundary_2, 2) - size(complex%boundary_3, 2)
    end subroutine cell_complex_euler_characteristic

    subroutine cell_complex_betti_numbers(complex, betti, status)
        type(cell_complex_t), intent(in) :: complex
        integer, intent(out) :: betti(4)
        integer, intent(out) :: status
        integer :: rank_1, rank_2, rank_3

        betti = 0
        call validate_cell_complex(complex, status)
        if (status /= 0) return
        call matrix_rank(complex%boundary_1, rank_1)
        call matrix_rank(complex%boundary_2, rank_2)
        call matrix_rank(complex%boundary_3, rank_3)
        betti = [ &
            complex%vertex_count - rank_1, &
            size(complex%boundary_1, 2) - rank_1 - rank_2, &
            size(complex%boundary_2, 2) - rank_2 - rank_3, &
            size(complex%boundary_3, 2) - rank_3]
        if (any(betti < 0)) status = 7
    end subroutine cell_complex_betti_numbers

    subroutine cell_complex_cycle_basis( &
            complex, cycles, cycle_count, status)
        !! Return a real basis of the one-cycle space ker(boundary_1).
        !!
        !! The columns are coefficient vectors on oriented edges and satisfy
        !! boundary_1*cycles = 0 up to the scale-aware elimination tolerance.
        !! These are cycles, not yet homology or metric-harmonic
        !! representatives: boundaries of faces and Hodge factors remain
        !! explicit in higher layers.
        type(cell_complex_t), intent(in) :: complex
        real(dp), allocatable, intent(out) :: cycles(:, :)
        integer, intent(out) :: cycle_count
        integer, intent(out) :: status

        if (allocated(cycles)) deallocate(cycles)
        cycle_count = 0
        call validate_cell_complex(complex, status)
        if (status /= 0) return
        call matrix_nullspace(complex%boundary_1, cycles, cycle_count)
        status = 0
    end subroutine cell_complex_cycle_basis

    subroutine quotient_cell_complex( &
            complex, vertex_identification, edge_identification, quotient, &
            status, face_identification, volume_identification)
        !! Build a quotient complex from signed representative maps.
        !!
        !! Each supplied identification is a chain-level map: the lower
        !! boundary of every identified upper cell must agree with the
        !! declared orientation of its representative.  This routine only
        !! composes those maps; geometry, metric Hodge factors, and mesh
        !! ownership remain outside the topology module.
        type(cell_complex_t), intent(in) :: complex
        type(cell_identification_t), intent(in) :: vertex_identification
        type(cell_identification_t), intent(in) :: edge_identification
        type(cell_complex_t), intent(inout) :: quotient
        integer, intent(out) :: status
        type(cell_identification_t), intent(in), optional :: face_identification
        type(cell_identification_t), intent(in), optional :: volume_identification
        integer, allocatable :: quotient_boundary_1(:, :)
        integer, allocatable :: quotient_boundary_2(:, :)
        integer, allocatable :: quotient_boundary_3(:, :)
        integer :: face_count, volume_count, validation_status

        call clear_cell_complex(quotient)
        status = 1
        call validate_cell_complex(complex, validation_status)
        if (validation_status /= 0) then
            status = 2
            return
        end if

        call identify_boundary_matrix( &
            vertex_identification, edge_identification, &
            complex%boundary_1, quotient_boundary_1, status)
        if (status /= 0) then
            status = 3
            return
        end if

        face_count = size(complex%boundary_2, 2)
        if (face_count > 0) then
            if (.not. present(face_identification)) then
                status = 4
                return
            end if
            call identify_boundary_matrix( &
                edge_identification, face_identification, &
                complex%boundary_2, quotient_boundary_2, status)
            if (status /= 0) then
                status = 5
                return
            end if
        else
            allocate(quotient_boundary_2(size(quotient_boundary_1, 1), 0))
        end if

        volume_count = size(complex%boundary_3, 2)
        if (volume_count > 0) then
            if (.not. present(face_identification) .or. &
                .not. present(volume_identification)) then
                status = 6
                return
            end if
            call identify_boundary_matrix( &
                face_identification, volume_identification, &
                complex%boundary_3, quotient_boundary_3, status)
            if (status /= 0) then
                status = 7
                return
            end if
        else
            allocate(quotient_boundary_3(size(quotient_boundary_2, 2), 0))
        end if

        call initialize_cell_complex(quotient, size(quotient_boundary_1, 1), &
            quotient_boundary_1, quotient_boundary_2, quotient_boundary_3, &
            status)
        if (status /= 0) then
            status = 8
            return
        end if
        call validate_cell_complex(quotient, validation_status)
        if (validation_status /= 0) then
            call clear_cell_complex(quotient)
            status = 9
            return
        end if
        status = 0
    end subroutine quotient_cell_complex

    subroutine clear_cell_complex(complex)
        type(cell_complex_t), intent(inout) :: complex

        if (allocated(complex%boundary_1)) deallocate(complex%boundary_1)
        if (allocated(complex%boundary_2)) deallocate(complex%boundary_2)
        if (allocated(complex%boundary_3)) deallocate(complex%boundary_3)
        complex%vertex_count = 0
    end subroutine clear_cell_complex

    subroutine matrix_rank(matrix, rank)
        integer, intent(in) :: matrix(:, :)
        integer, intent(out) :: rank
        real(dp), allocatable :: work(:, :)
        real(dp) :: scale, tolerance, pivot_value, factor, value
        integer :: row, column, pivot_row, pivot, rows, columns

        rows = size(matrix, 1)
        columns = size(matrix, 2)
        rank = 0
        if (rows == 0 .or. columns == 0) return
        allocate(work(rows, columns))
        work = real(matrix, dp)
        scale = max(1.0_dp, maxval(abs(work)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale* &
            real(max(rows, columns), dp)
        pivot_row = 1
        do column = 1, columns
            if (pivot_row > rows) exit
            pivot = pivot_row
            pivot_value = abs(work(pivot_row, column))
            do row = pivot_row + 1, rows
                value = abs(work(row, column))
                if (value > pivot_value) then
                    pivot = row
                    pivot_value = value
                end if
            end do
            if (pivot_value <= tolerance) cycle
            if (pivot /= pivot_row) then
                work([pivot, pivot_row], :) = work([pivot_row, pivot], :)
            end if
            pivot_value = work(pivot_row, column)
            work(pivot_row, :) = work(pivot_row, :)/pivot_value
            do row = 1, rows
                if (row == pivot_row) cycle
                factor = work(row, column)
                if (abs(factor) > tolerance) then
                    work(row, :) = work(row, :) - &
                        factor*work(pivot_row, :)
                end if
            end do
            rank = rank + 1
            pivot_row = pivot_row + 1
        end do
    end subroutine matrix_rank

    subroutine matrix_nullspace(matrix, basis, nullity)
        integer, intent(in) :: matrix(:, :)
        real(dp), allocatable, intent(out) :: basis(:, :)
        integer, intent(out) :: nullity
        real(dp), allocatable :: work(:, :)
        real(dp) :: scale, tolerance, pivot_value, factor, value
        integer, allocatable :: pivot_columns(:)
        integer :: row, column, pivot_row, pivot, rows, columns
        integer :: free_column, free_index

        if (allocated(basis)) deallocate(basis)
        rows = size(matrix, 1)
        columns = size(matrix, 2)
        nullity = columns
        allocate(basis(columns, columns))
        basis = 0.0_dp
        if (columns == 0) then
            deallocate(basis)
            allocate(basis(0, 0))
            nullity = 0
            return
        end if

        allocate(work(rows, columns), pivot_columns(rows))
        work = real(matrix, dp)
        pivot_columns = 0
        pivot_row = 0
        scale = 1.0_dp
        if (rows > 0) scale = max(scale, maxval(abs(work)))
        tolerance = 128.0_dp*epsilon(1.0_dp)*scale* &
            real(max(rows, columns), dp)
        do column = 1, columns
            if (pivot_row >= rows) exit
            pivot = pivot_row + 1
            pivot_value = abs(work(pivot, column))
            do row = pivot + 1, rows
                value = abs(work(row, column))
                if (value > pivot_value) then
                    pivot = row
                    pivot_value = value
                end if
            end do
            if (pivot_value <= tolerance) cycle
            pivot_row = pivot_row + 1
            if (pivot /= pivot_row) then
                work([pivot, pivot_row], :) = work([pivot_row, pivot], :)
            end if
            pivot_value = work(pivot_row, column)
            work(pivot_row, :) = work(pivot_row, :)/pivot_value
            do row = 1, rows
                if (row == pivot_row) cycle
                factor = work(row, column)
                if (abs(factor) > tolerance) then
                    work(row, :) = work(row, :) - &
                        factor*work(pivot_row, :)
                end if
            end do
            pivot_columns(pivot_row) = column
        end do

        nullity = columns - pivot_row
        deallocate(basis)
        allocate(basis(columns, nullity))
        basis = 0.0_dp
        free_index = 0
        do free_column = 1, columns
            if (pivot_row > 0) then
                if (any(pivot_columns(:pivot_row) == free_column)) cycle
            end if
            free_index = free_index + 1
            basis(free_column, free_index) = 1.0_dp
            do row = 1, pivot_row
                basis(pivot_columns(row), free_index) = -work(row, free_column)
            end do
        end do
    end subroutine matrix_nullspace

end module fortfem_cell_complex
