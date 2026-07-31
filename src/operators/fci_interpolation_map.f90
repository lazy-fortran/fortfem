module fortfem_fci_interpolation_map
    !! Small geometry-side interpolation contracts for FCI support operators.
    !!
    !! The routine here deliberately handles only one-dimensional, piecewise
    !! linear stencils.  A field-line/map service can use it for a coordinate
    !! slice or as an oracle for a higher-dimensional interpolation builder;
    !! it does not make assumptions about magnetic-field storage.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: build_fci_linear_interpolation_map_1d

contains

    subroutine build_fci_linear_interpolation_map_1d( &
            source_coordinates, target_coordinates, interpolation_map, status)
        !! Build row-wise linear interpolation weights.
        !!
        !! `interpolation_map(row, :)` maps values at strictly increasing
        !! `source_coordinates` to `target_coordinates(row)`.  Targets must
        !! lie in the closed source interval.  Endpoint rows are represented
        !! exactly, and every valid row has at most two nonzero weights whose
        !! sum is one.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        real(dp), intent(out) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row, left
        real(dp) :: alpha, denominator

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI interpolation received incompatible coordinates")
        interpolation_map = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 2 .or. target_count < 1) return
        if (size(interpolation_map, 1) /= target_count .or. &
            size(interpolation_map, 2) /= source_count) return
        if (any(.not. ieee_is_finite(source_coordinates)) .or. &
            any(.not. ieee_is_finite(target_coordinates))) return
        do left = 1, source_count - 1
            if (source_coordinates(left + 1) <= source_coordinates(left)) return
        end do

        do row = 1, target_count
            if (target_coordinates(row) < source_coordinates(1) .or. &
                target_coordinates(row) > source_coordinates(source_count)) then
                return
            end if
            if (target_coordinates(row) == source_coordinates(1)) then
                interpolation_map(row, 1) = 1.0_dp
            else if (target_coordinates(row) == &
                    source_coordinates(source_count)) then
                interpolation_map(row, source_count) = 1.0_dp
            else
                left = 0
                do while (left < source_count - 1)
                    left = left + 1
                    if (target_coordinates(row) <= source_coordinates(left + 1)) &
                        exit
                end do
                if (left < 1 .or. left >= source_count) return
                denominator = source_coordinates(left + 1) - &
                    source_coordinates(left)
                if (denominator <= 0.0_dp) return
                alpha = (target_coordinates(row) - source_coordinates(left))/ &
                    denominator
                interpolation_map(row, left) = 1.0_dp - alpha
                interpolation_map(row, left + 1) = alpha
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_linear_interpolation_map_1d

end module fortfem_fci_interpolation_map
