module fortfem_assembly_bspline_polar_2d
    !! Sparse Galerkin restriction from periodic tensor to polar bases.
    use fortfem_kinds, only: dp
    use fortsparse, only: &
        csc_from_triplet, csc_matmul, csc_t, csc_transpose, &
        FORTSPARSE_INVALID_MATRIX, fortsparse_status_t, status_set
    implicit none
    private

    public :: restrict_bspline_polar_operator_csc

contains

    subroutine restrict_bspline_polar_operator_csc( &
            extraction, tensor_operator, polar_operator, status)
        !! Compute E A E^T for a tensor operator A and polar basis extraction E.
        real(dp), intent(in) :: extraction(:, :)
        type(csc_t), intent(in) :: tensor_operator
        type(csc_t), intent(out) :: polar_operator
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: values(:)
        type(csc_t) :: extraction_csc, extraction_transpose, work
        integer :: column, entry, row

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Polar Galerkin restriction requires compatible dimensions")
        if (size(extraction, 2) /= tensor_operator%nrow .or. &
            tensor_operator%nrow /= tensor_operator%ncol) return
        allocate( &
            rows(count(extraction /= 0.0_dp)), &
            columns(count(extraction /= 0.0_dp)), &
            values(count(extraction /= 0.0_dp)))
        entry = 0
        do column = 1, size(extraction, 2)
            do row = 1, size(extraction, 1)
                if (extraction(row, column) == 0.0_dp) cycle
                entry = entry + 1
                rows(entry) = row
                columns(entry) = column
                values(entry) = extraction(row, column)
            end do
        end do
        call csc_from_triplet( &
            size(extraction, 1), size(extraction, 2), rows, columns, values, &
            extraction_csc, status)
        if (status%code /= 0) return
        call csc_matmul(extraction_csc, tensor_operator, work, status)
        if (status%code /= 0) return
        call csc_transpose(extraction_csc, extraction_transpose, status)
        if (status%code /= 0) return
        call csc_matmul(work, extraction_transpose, polar_operator, status)
    end subroutine restrict_bspline_polar_operator_csc

end module fortfem_assembly_bspline_polar_2d
