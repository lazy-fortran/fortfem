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
    public :: build_fci_bilinear_interpolation_map_2d
    public :: build_fci_bilinear_interpolation_map_2d_jvp
    public :: build_fci_bilinear_interpolation_map_2d_vjp

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

    subroutine build_fci_bilinear_interpolation_map_2d( &
            source_x, source_y, target_x, target_y, interpolation_map, status)
        !! Build tensor-product bilinear map rows on a Cartesian source grid.
        !!
        !! Source columns use x as the fastest index:
        !! `column = (y_index - 1)*size(source_x) + x_index`.
        !! Targets on source grid lines are represented exactly; all valid
        !! rows have at most four nonzero weights and sum to one.
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: target_x(:)
        real(dp), intent(in) :: target_y(:)
        real(dp), intent(out) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nx, ny, target_count, row, ix, iy, column
        real(dp) :: alpha_x, alpha_y
        logical :: valid_x, valid_y

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI bilinear interpolation received incompatible arrays")
        interpolation_map = 0.0_dp
        nx = size(source_x)
        ny = size(source_y)
        target_count = size(target_x)
        if (nx < 2 .or. ny < 2 .or. target_count < 1) return
        if (size(target_y) /= target_count) return
        if (size(interpolation_map, 1) /= target_count .or. &
            size(interpolation_map, 2) /= nx*ny) return
        if (any(.not. ieee_is_finite(source_x)) .or. &
            any(.not. ieee_is_finite(source_y)) .or. &
            any(.not. ieee_is_finite(target_x)) .or. &
            any(.not. ieee_is_finite(target_y))) return
        if (.not. strictly_increasing(source_x) .or. &
            .not. strictly_increasing(source_y)) return

        do row = 1, target_count
            call find_linear_bracket( &
                source_x, target_x(row), ix, alpha_x, valid_x)
            call find_linear_bracket( &
                source_y, target_y(row), iy, alpha_y, valid_y)
            if (.not. valid_x .or. .not. valid_y) return
            column = (iy - 1)*nx + ix
            interpolation_map(row, column) = &
                interpolation_map(row, column) + (1.0_dp - alpha_x)* &
                (1.0_dp - alpha_y)
            interpolation_map(row, column + 1) = &
                interpolation_map(row, column + 1) + alpha_x* &
                (1.0_dp - alpha_y)
            interpolation_map(row, column + nx) = &
                interpolation_map(row, column + nx) + (1.0_dp - alpha_x)* &
                alpha_y
            interpolation_map(row, column + nx + 1) = &
                interpolation_map(row, column + nx + 1) + alpha_x*alpha_y
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_bilinear_interpolation_map_2d

    subroutine build_fci_bilinear_interpolation_map_2d_jvp( &
            source_x, source_y, target_x, target_y, source_x_dot, source_y_dot, &
            target_x_dot, target_y_dot, interpolation_map_dot, status)
        !! Differentiate bilinear weights on a fixed Cartesian cell topology.
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: target_x(:)
        real(dp), intent(in) :: target_y(:)
        real(dp), intent(in) :: source_x_dot(:)
        real(dp), intent(in) :: source_y_dot(:)
        real(dp), intent(in) :: target_x_dot(:)
        real(dp), intent(in) :: target_y_dot(:)
        real(dp), intent(out) :: interpolation_map_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nx, ny, target_count, row, ix, iy, column
        real(dp) :: alpha_x, alpha_y, alpha_x_dot, alpha_y_dot
        real(dp) :: numerator_x, numerator_y, denominator_x, denominator_y
        real(dp) :: numerator_x_dot, numerator_y_dot
        real(dp) :: denominator_x_dot, denominator_y_dot
        real(dp), allocatable :: interpolation_map(:, :)
        logical :: valid_x, valid_y

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI bilinear JVP received incompatible arrays")
        interpolation_map_dot = 0.0_dp
        nx = size(source_x)
        ny = size(source_y)
        target_count = size(target_x)
        if (nx < 2 .or. ny < 2 .or. target_count < 1) return
        allocate(interpolation_map(target_count, nx*ny))
        call build_fci_bilinear_interpolation_map_2d( &
            source_x, source_y, target_x, target_y, interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(source_x_dot) /= nx .or. size(source_y_dot) /= ny .or. &
            size(target_x_dot) /= target_count .or. &
            size(target_y_dot) /= target_count .or. &
            size(interpolation_map_dot, 1) /= target_count .or. &
            size(interpolation_map_dot, 2) /= nx*ny) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI bilinear JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(source_x_dot)) .or. &
            any(.not. ieee_is_finite(source_y_dot)) .or. &
            any(.not. ieee_is_finite(target_x_dot)) .or. &
            any(.not. ieee_is_finite(target_y_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI bilinear JVP received non-finite tangents")
            return
        end if

        do row = 1, target_count
            call locate_fixed_interval( &
                source_x, target_x(row), ix, valid_x)
            call locate_fixed_interval( &
                source_y, target_y(row), iy, valid_y)
            if (.not. valid_x .or. .not. valid_y) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI bilinear JVP crosses a grid-line topology event")
                return
            end if
            numerator_x = target_x(row) - source_x(ix)
            denominator_x = source_x(ix + 1) - source_x(ix)
            numerator_x_dot = target_x_dot(row) - source_x_dot(ix)
            denominator_x_dot = source_x_dot(ix + 1) - source_x_dot(ix)
            alpha_x = numerator_x/denominator_x
            alpha_x_dot = (numerator_x_dot*denominator_x - &
                numerator_x*denominator_x_dot)/(denominator_x*denominator_x)
            numerator_y = target_y(row) - source_y(iy)
            denominator_y = source_y(iy + 1) - source_y(iy)
            numerator_y_dot = target_y_dot(row) - source_y_dot(iy)
            denominator_y_dot = source_y_dot(iy + 1) - source_y_dot(iy)
            alpha_y = numerator_y/denominator_y
            alpha_y_dot = (numerator_y_dot*denominator_y - &
                numerator_y*denominator_y_dot)/(denominator_y*denominator_y)
            column = (iy - 1)*nx + ix
            interpolation_map_dot(row, column) = &
                -alpha_x_dot*(1.0_dp - alpha_y) - &
                alpha_y_dot*(1.0_dp - alpha_x)
            interpolation_map_dot(row, column + 1) = &
                alpha_x_dot*(1.0_dp - alpha_y) - alpha_x*alpha_y_dot
            interpolation_map_dot(row, column + nx) = &
                -alpha_x_dot*alpha_y + (1.0_dp - alpha_x)*alpha_y_dot
            interpolation_map_dot(row, column + nx + 1) = &
                alpha_x_dot*alpha_y + alpha_x*alpha_y_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_bilinear_interpolation_map_2d_jvp

    subroutine build_fci_bilinear_interpolation_map_2d_vjp( &
            source_x, source_y, target_x, target_y, interpolation_map_bar, &
            source_x_bar, source_y_bar, target_x_bar, target_y_bar, status)
        !! Apply the VJP of fixed-topology bilinear interpolation weights.
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: target_x(:)
        real(dp), intent(in) :: target_y(:)
        real(dp), intent(in) :: interpolation_map_bar(:, :)
        real(dp), intent(out) :: source_x_bar(:)
        real(dp), intent(out) :: source_y_bar(:)
        real(dp), intent(out) :: target_x_bar(:)
        real(dp), intent(out) :: target_y_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nx, ny, target_count, row, ix, iy, column
        real(dp) :: alpha_x, alpha_y, alpha_x_bar, alpha_y_bar
        real(dp) :: numerator_x, numerator_y, denominator_x, denominator_y
        real(dp), allocatable :: interpolation_map(:, :)
        logical :: valid_x, valid_y

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI bilinear VJP received incompatible arrays")
        source_x_bar = 0.0_dp
        source_y_bar = 0.0_dp
        target_x_bar = 0.0_dp
        target_y_bar = 0.0_dp
        nx = size(source_x)
        ny = size(source_y)
        target_count = size(target_x)
        if (nx < 2 .or. ny < 2 .or. target_count < 1) return
        allocate(interpolation_map(target_count, nx*ny))
        call build_fci_bilinear_interpolation_map_2d( &
            source_x, source_y, target_x, target_y, interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_bar, 1) /= target_count .or. &
            size(interpolation_map_bar, 2) /= nx*ny .or. &
            size(source_x_bar) /= nx .or. size(source_y_bar) /= ny .or. &
            size(target_x_bar) /= target_count .or. &
            size(target_y_bar) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI bilinear VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(interpolation_map_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI bilinear VJP received non-finite cotangents")
            return
        end if

        do row = 1, target_count
            call locate_fixed_interval( &
                source_x, target_x(row), ix, valid_x)
            call locate_fixed_interval( &
                source_y, target_y(row), iy, valid_y)
            if (.not. valid_x .or. .not. valid_y) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI bilinear VJP crosses a grid-line topology event")
                return
            end if
            numerator_x = target_x(row) - source_x(ix)
            denominator_x = source_x(ix + 1) - source_x(ix)
            alpha_x = numerator_x/denominator_x
            numerator_y = target_y(row) - source_y(iy)
            denominator_y = source_y(iy + 1) - source_y(iy)
            alpha_y = numerator_y/denominator_y
            column = (iy - 1)*nx + ix
            alpha_x_bar = (interpolation_map_bar(row, column + 1) - &
                interpolation_map_bar(row, column))*(1.0_dp - alpha_y) + &
                (interpolation_map_bar(row, column + nx + 1) - &
                interpolation_map_bar(row, column + nx))*alpha_y
            alpha_y_bar = (interpolation_map_bar(row, column + nx) - &
                interpolation_map_bar(row, column))*(1.0_dp - alpha_x) + &
                (interpolation_map_bar(row, column + nx + 1) - &
                interpolation_map_bar(row, column + 1))*alpha_x
            target_x_bar(row) = target_x_bar(row) + &
                alpha_x_bar/denominator_x
            source_x_bar(ix) = source_x_bar(ix) + alpha_x_bar* &
                (numerator_x - denominator_x)/(denominator_x*denominator_x)
            source_x_bar(ix + 1) = source_x_bar(ix + 1) - alpha_x_bar* &
                numerator_x/(denominator_x*denominator_x)
            target_y_bar(row) = target_y_bar(row) + &
                alpha_y_bar/denominator_y
            source_y_bar(iy) = source_y_bar(iy) + alpha_y_bar* &
                (numerator_y - denominator_y)/(denominator_y*denominator_y)
            source_y_bar(iy + 1) = source_y_bar(iy + 1) - alpha_y_bar* &
                numerator_y/(denominator_y*denominator_y)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_bilinear_interpolation_map_2d_vjp

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

    logical function strictly_increasing(coordinates)
        real(dp), intent(in) :: coordinates(:)
        integer :: index

        strictly_increasing = size(coordinates) >= 2
        if (.not. strictly_increasing) return
        do index = 1, size(coordinates) - 1
            if (coordinates(index + 1) <= coordinates(index)) then
                strictly_increasing = .false.
                return
            end if
        end do
    end function strictly_increasing

    subroutine find_linear_bracket( &
            coordinates, coordinate, left, alpha, valid)
        real(dp), intent(in) :: coordinates(:)
        real(dp), intent(in) :: coordinate
        integer, intent(out) :: left
        real(dp), intent(out) :: alpha
        logical, intent(out) :: valid

        left = 0
        alpha = 0.0_dp
        valid = .false.
        if (coordinate < coordinates(1) .or. &
            coordinate > coordinates(size(coordinates))) return
        if (coordinate == coordinates(size(coordinates))) then
            left = size(coordinates) - 1
            alpha = 1.0_dp
        else
            do left = 1, size(coordinates) - 1
                if (coordinate <= coordinates(left + 1)) exit
            end do
            if (left < 1 .or. left >= size(coordinates)) then
                left = 0
                return
            end if
            alpha = (coordinate - coordinates(left))/ &
                (coordinates(left + 1) - coordinates(left))
        end if
        valid = .true.
    end subroutine find_linear_bracket

end module fortfem_fci_interpolation_map
