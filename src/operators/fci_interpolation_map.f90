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
    public :: build_fci_linear_interpolation_map_1d_jvp
    public :: build_fci_linear_interpolation_map_1d_vjp

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

    subroutine build_fci_linear_interpolation_map_1d_jvp( &
            source_coordinates, target_coordinates, source_coordinates_dot, &
            target_coordinates_dot, interpolation_map_dot, status)
        !! Differentiate the map on a fixed interpolation-cell topology.
        !!
        !! A target exactly on a source node is a topology event and is
        !! rejected rather than assigned a misleading one-sided derivative.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        real(dp), intent(in) :: source_coordinates_dot(:)
        real(dp), intent(in) :: target_coordinates_dot(:)
        real(dp), intent(out) :: interpolation_map_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row, left
        real(dp) :: numerator, denominator, numerator_dot, denominator_dot
        real(dp) :: alpha_dot
        real(dp), allocatable :: interpolation_map(:, :)
        logical :: valid_interval

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI interpolation JVP received incompatible arrays")
        interpolation_map_dot = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 2 .or. target_count < 1) return
        allocate(interpolation_map(target_count, source_count))
        call build_fci_linear_interpolation_map_1d( &
            source_coordinates, target_coordinates, interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(source_coordinates_dot) /= source_count .or. &
            size(target_coordinates_dot) /= target_count .or. &
            size(interpolation_map_dot, 1) /= target_count .or. &
            size(interpolation_map_dot, 2) /= source_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI interpolation JVP received incompatible tangents")
            return
        end if

        do row = 1, target_count
            call locate_fixed_interval( &
                source_coordinates, target_coordinates(row), left, &
                valid_interval)
            if (.not. valid_interval) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI interpolation JVP crosses a stencil topology event")
                return
            end if
            numerator = target_coordinates(row) - source_coordinates(left)
            denominator = source_coordinates(left + 1) - &
                source_coordinates(left)
            numerator_dot = target_coordinates_dot(row) - &
                source_coordinates_dot(left)
            denominator_dot = source_coordinates_dot(left + 1) - &
                source_coordinates_dot(left)
            alpha_dot = (numerator_dot*denominator - &
                numerator*denominator_dot)/(denominator*denominator)
            interpolation_map_dot(row, left) = -alpha_dot
            interpolation_map_dot(row, left + 1) = alpha_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_linear_interpolation_map_1d_jvp

    subroutine build_fci_linear_interpolation_map_1d_vjp( &
            source_coordinates, target_coordinates, interpolation_map_bar, &
            source_coordinates_bar, target_coordinates_bar, status)
        !! Apply the VJP of the fixed-topology interpolation map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        real(dp), intent(in) :: interpolation_map_bar(:, :)
        real(dp), intent(out) :: source_coordinates_bar(:)
        real(dp), intent(out) :: target_coordinates_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row, left
        real(dp) :: numerator, denominator, alpha_bar
        real(dp), allocatable :: interpolation_map(:, :)
        logical :: valid_interval

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI interpolation VJP received incompatible arrays")
        source_coordinates_bar = 0.0_dp
        target_coordinates_bar = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 2 .or. target_count < 1) return
        allocate(interpolation_map(target_count, source_count))
        call build_fci_linear_interpolation_map_1d( &
            source_coordinates, target_coordinates, interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_bar, 1) /= target_count .or. &
            size(interpolation_map_bar, 2) /= source_count .or. &
            size(source_coordinates_bar) /= source_count .or. &
            size(target_coordinates_bar) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI interpolation VJP received incompatible cotangents")
            return
        end if

        do row = 1, target_count
            call locate_fixed_interval( &
                source_coordinates, target_coordinates(row), left, &
                valid_interval)
            if (.not. valid_interval) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI interpolation VJP crosses a stencil topology event")
                return
            end if
            numerator = target_coordinates(row) - source_coordinates(left)
            denominator = source_coordinates(left + 1) - &
                source_coordinates(left)
            alpha_bar = interpolation_map_bar(row, left + 1) - &
                interpolation_map_bar(row, left)
            target_coordinates_bar(row) = target_coordinates_bar(row) + &
                alpha_bar/denominator
            source_coordinates_bar(left) = source_coordinates_bar(left) + &
                alpha_bar*(numerator - denominator)/(denominator*denominator)
            source_coordinates_bar(left + 1) = &
                source_coordinates_bar(left + 1) - &
                alpha_bar*numerator/(denominator*denominator)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_linear_interpolation_map_1d_vjp

    subroutine locate_fixed_interval( &
            source_coordinates, target_coordinate, left, valid_interval)
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinate
        integer, intent(out) :: left
        logical, intent(out) :: valid_interval

        left = 0
        valid_interval = .false.
        do left = 1, size(source_coordinates) - 1
            if (target_coordinate > source_coordinates(left) .and. &
                target_coordinate < source_coordinates(left + 1)) then
                valid_interval = .true.
                return
            end if
        end do
        left = 0
    end subroutine locate_fixed_interval

end module fortfem_fci_interpolation_map
