module fortfem_fci_interpolation_map
    !! Small geometry-side interpolation contracts for FCI support operators.
    !!
    !! The routines here handle one-dimensional piecewise-linear and explicit
    !! quadratic/cubic stencils. A field-line/map service can use them for a
    !! coordinate slice or as an oracle for a higher-dimensional interpolation
    !! builder; they do not make assumptions about magnetic-field storage.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_generated_fci_quadratic_lagrange, only: &
        generated_fci_quadratic_lagrange_weights
    use fortfem_generated_fci_quadratic_lagrange_jvp, only: &
        generated_fci_quadratic_lagrange_weights_jvp
    use fortfem_generated_fci_quadratic_lagrange_vjp, only: &
        generated_fci_quadratic_lagrange_weights_vjp
    use fortfem_generated_fci_cubic_lagrange, only: &
        generated_fci_cubic_lagrange_weights
    use fortfem_generated_fci_cubic_lagrange_jvp, only: &
        generated_fci_cubic_lagrange_weights_jvp
    use fortfem_generated_fci_cubic_lagrange_vjp, only: &
        generated_fci_cubic_lagrange_weights_vjp
    use fortfem_generated_fci_quartic_lagrange, only: &
        generated_fci_quartic_lagrange_weights
    use fortfem_generated_fci_quartic_lagrange_jvp, only: &
        generated_fci_quartic_lagrange_weights_jvp
    use fortfem_generated_fci_quartic_lagrange_vjp, only: &
        generated_fci_quartic_lagrange_weights_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: build_fci_linear_interpolation_map_1d
    public :: build_fci_linear_interpolation_map_1d_jvp
    public :: build_fci_linear_interpolation_map_1d_vjp
    public :: build_fci_quadratic_interpolation_map_1d
    public :: build_fci_quadratic_interpolation_map_1d_jvp
    public :: build_fci_quadratic_interpolation_map_1d_vjp
    public :: build_fci_cubic_interpolation_map_1d
    public :: build_fci_cubic_interpolation_map_1d_jvp
    public :: build_fci_cubic_interpolation_map_1d_vjp
    public :: build_fci_quartic_interpolation_map_1d
    public :: build_fci_quartic_interpolation_map_1d_jvp
    public :: build_fci_quartic_interpolation_map_1d_vjp
    public :: build_fci_quadratic_interpolation_maps_1d
    public :: build_fci_quadratic_interpolation_maps_1d_jvp
    public :: build_fci_quadratic_interpolation_maps_1d_vjp
    public :: build_fci_bilinear_interpolation_map_2d
    public :: build_fci_bilinear_interpolation_map_2d_jvp
    public :: build_fci_bilinear_interpolation_map_2d_vjp
    public :: build_fci_bilinear_interpolation_maps_2d
    public :: build_fci_bilinear_interpolation_maps_2d_jvp
    public :: build_fci_bilinear_interpolation_maps_2d_vjp
    public :: build_fci_triangle_interpolation_map_2d
    public :: build_fci_triangle_interpolation_map_2d_jvp
    public :: build_fci_triangle_interpolation_map_2d_vjp
    public :: build_fci_triangle_interpolation_maps_2d
    public :: build_fci_triangle_interpolation_maps_2d_jvp
    public :: build_fci_triangle_interpolation_maps_2d_vjp

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

    subroutine build_fci_quadratic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        !! Build row-wise quadratic Lagrange interpolation weights.
        !!
        !! `stencil_indices(:, row)` supplies three distinct source nodes for
        !! each target.  The target must lie in the closed interval spanned by
        !! its local stencil; this keeps the builder an interpolation contract
        !! while allowing nonuniform and logically unstructured source slices.
        !! The local three-node weight kernel is generated by FortSym.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(out) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: first_node, second_node, third_node
        real(dp) :: weight_1, weight_2, weight_3
        real(dp) :: nodes(3)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI quadratic interpolation received incompatible coordinates")
        interpolation_map = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 3 .or. target_count < 1) return
        if (size(stencil_indices, 1) /= 3 .or. &
            size(stencil_indices, 2) /= target_count) return
        if (size(interpolation_map, 1) /= target_count .or. &
            size(interpolation_map, 2) /= source_count) return
        if (any(.not. ieee_is_finite(source_coordinates)) .or. &
            any(.not. ieee_is_finite(target_coordinates))) return

        do row = 1, target_count
            first_node = stencil_indices(1, row)
            second_node = stencil_indices(2, row)
            third_node = stencil_indices(3, row)
            if (first_node < 1 .or. first_node > source_count .or. &
                second_node < 1 .or. second_node > source_count .or. &
                third_node < 1 .or. third_node > source_count) return
            if (first_node == second_node .or. first_node == third_node .or. &
                second_node == third_node) return
            nodes = source_coordinates([first_node, second_node, third_node])
            if (nodes(1) == nodes(2) .or. nodes(1) == nodes(3) .or. &
                nodes(2) == nodes(3)) return
            if (target_coordinates(row) < minval(nodes) .or. &
                target_coordinates(row) > maxval(nodes)) return
            call generated_fci_quadratic_lagrange_weights( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), &
                weight_1, weight_2, weight_3)
            interpolation_map(row, first_node) = weight_1
            interpolation_map(row, second_node) = weight_2
            interpolation_map(row, third_node) = weight_3
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quadratic_interpolation_map_1d

    subroutine build_fci_quadratic_interpolation_map_1d_jvp( &
            source_coordinates, target_coordinates, stencil_indices, &
            source_coordinates_dot, target_coordinates_dot, &
            interpolation_map_dot, status)
        !! Apply the fixed-stencil JVP of the quadratic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(in) :: source_coordinates_dot(:)
        real(dp), intent(in) :: target_coordinates_dot(:)
        real(dp), intent(out) :: interpolation_map_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: first_node, second_node, third_node
        real(dp) :: nodes(3), nodes_dot(3)
        real(dp) :: weight_1_dot, weight_2_dot, weight_3_dot
        real(dp), allocatable :: interpolation_map(:, :)

        interpolation_map_dot = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 3 .or. target_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadratic JVP received incompatible coordinates")
            return
        end if
        allocate(interpolation_map(target_count, source_count))
        call build_fci_quadratic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(source_coordinates_dot) /= source_count .or. &
            size(target_coordinates_dot) /= target_count .or. &
            size(interpolation_map_dot, 1) /= target_count .or. &
            size(interpolation_map_dot, 2) /= source_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadratic JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(source_coordinates_dot)) .or. &
            any(.not. ieee_is_finite(target_coordinates_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadratic JVP received non-finite tangents")
            return
        end if

        do row = 1, target_count
            first_node = stencil_indices(1, row)
            second_node = stencil_indices(2, row)
            third_node = stencil_indices(3, row)
            nodes = source_coordinates([first_node, second_node, third_node])
            nodes_dot = source_coordinates_dot( &
                [first_node, second_node, third_node])
            call generated_fci_quadratic_lagrange_weights_jvp( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), &
                target_coordinates_dot(row), nodes_dot(1), nodes_dot(2), &
                nodes_dot(3), weight_1_dot, weight_2_dot, weight_3_dot)
            interpolation_map_dot(row, first_node) = weight_1_dot
            interpolation_map_dot(row, second_node) = weight_2_dot
            interpolation_map_dot(row, third_node) = weight_3_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quadratic_interpolation_map_1d_jvp

    subroutine build_fci_quadratic_interpolation_map_1d_vjp( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map_bar, source_coordinates_bar, &
            target_coordinates_bar, status)
        !! Apply the real VJP of the fixed-stencil quadratic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(in) :: interpolation_map_bar(:, :)
        real(dp), intent(out) :: source_coordinates_bar(:)
        real(dp), intent(out) :: target_coordinates_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: first_node, second_node, third_node
        real(dp) :: nodes(3), weight_bar(3)
        real(dp) :: target_bar, node_1_bar, node_2_bar, node_3_bar
        real(dp), allocatable :: interpolation_map(:, :)

        source_coordinates_bar = 0.0_dp
        target_coordinates_bar = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 3 .or. target_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadratic VJP received incompatible coordinates")
            return
        end if
        allocate(interpolation_map(target_count, source_count))
        call build_fci_quadratic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_bar, 1) /= target_count .or. &
            size(interpolation_map_bar, 2) /= source_count .or. &
            size(source_coordinates_bar) /= source_count .or. &
            size(target_coordinates_bar) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadratic VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(interpolation_map_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quadratic VJP received non-finite cotangents")
            return
        end if

        do row = 1, target_count
            first_node = stencil_indices(1, row)
            second_node = stencil_indices(2, row)
            third_node = stencil_indices(3, row)
            nodes = source_coordinates([first_node, second_node, third_node])
            weight_bar = [ &
                interpolation_map_bar(row, first_node), &
                interpolation_map_bar(row, second_node), &
                interpolation_map_bar(row, third_node)]
            call generated_fci_quadratic_lagrange_weights_vjp( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), &
                weight_bar(1), weight_bar(2), weight_bar(3), target_bar, &
                node_1_bar, node_2_bar, node_3_bar)
            target_coordinates_bar(row) = target_coordinates_bar(row) + &
                target_bar
            source_coordinates_bar(first_node) = &
                source_coordinates_bar(first_node) + node_1_bar
            source_coordinates_bar(second_node) = &
                source_coordinates_bar(second_node) + node_2_bar
            source_coordinates_bar(third_node) = &
                source_coordinates_bar(third_node) + node_3_bar
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quadratic_interpolation_map_1d_vjp

    subroutine build_fci_cubic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        !! Build row-wise cubic Lagrange interpolation weights.
        !!
        !! `stencil_indices(:, row)` supplies four distinct source nodes for
        !! each target.  The local weight kernel is generated by FortSym.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(out) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: first_node, second_node, third_node, fourth_node
        real(dp) :: nodes(4), weights(4)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI cubic interpolation received incompatible coordinates")
        interpolation_map = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 4 .or. target_count < 1) return
        if (size(stencil_indices, 1) /= 4 .or. &
            size(stencil_indices, 2) /= target_count) return
        if (size(interpolation_map, 1) /= target_count .or. &
            size(interpolation_map, 2) /= source_count) return
        if (any(.not. ieee_is_finite(source_coordinates)) .or. &
            any(.not. ieee_is_finite(target_coordinates))) return

        do row = 1, target_count
            first_node = stencil_indices(1, row)
            second_node = stencil_indices(2, row)
            third_node = stencil_indices(3, row)
            fourth_node = stencil_indices(4, row)
            if (first_node < 1 .or. first_node > source_count .or. &
                second_node < 1 .or. second_node > source_count .or. &
                third_node < 1 .or. third_node > source_count .or. &
                fourth_node < 1 .or. fourth_node > source_count) return
            if (first_node == second_node .or. first_node == third_node .or. &
                first_node == fourth_node .or. second_node == third_node .or. &
                second_node == fourth_node .or. third_node == fourth_node) return
            nodes = source_coordinates([first_node, second_node, third_node, &
                fourth_node])
            if (nodes(1) == nodes(2) .or. nodes(1) == nodes(3) .or. &
                nodes(1) == nodes(4) .or. nodes(2) == nodes(3) .or. &
                nodes(2) == nodes(4) .or. nodes(3) == nodes(4)) return
            if (target_coordinates(row) < minval(nodes) .or. &
                target_coordinates(row) > maxval(nodes)) return
            call generated_fci_cubic_lagrange_weights( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), nodes(4), &
                weights(1), weights(2), weights(3), weights(4))
            interpolation_map(row, first_node) = weights(1)
            interpolation_map(row, second_node) = weights(2)
            interpolation_map(row, third_node) = weights(3)
            interpolation_map(row, fourth_node) = weights(4)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_cubic_interpolation_map_1d

    subroutine build_fci_cubic_interpolation_map_1d_jvp( &
            source_coordinates, target_coordinates, stencil_indices, &
            source_coordinates_dot, target_coordinates_dot, &
            interpolation_map_dot, status)
        !! Apply the fixed-stencil JVP of the cubic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(in) :: source_coordinates_dot(:)
        real(dp), intent(in) :: target_coordinates_dot(:)
        real(dp), intent(out) :: interpolation_map_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: first_node, second_node, third_node, fourth_node
        real(dp) :: nodes(4), nodes_dot(4), weights_dot(4)
        real(dp), allocatable :: interpolation_map(:, :)

        interpolation_map_dot = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 4 .or. target_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI cubic JVP received incompatible coordinates")
            return
        end if
        allocate(interpolation_map(target_count, source_count))
        call build_fci_cubic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(source_coordinates_dot) /= source_count .or. &
            size(target_coordinates_dot) /= target_count .or. &
            size(interpolation_map_dot, 1) /= target_count .or. &
            size(interpolation_map_dot, 2) /= source_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI cubic JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(source_coordinates_dot)) .or. &
            any(.not. ieee_is_finite(target_coordinates_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI cubic JVP received non-finite tangents")
            return
        end if

        do row = 1, target_count
            first_node = stencil_indices(1, row)
            second_node = stencil_indices(2, row)
            third_node = stencil_indices(3, row)
            fourth_node = stencil_indices(4, row)
            nodes = source_coordinates([first_node, second_node, third_node, &
                fourth_node])
            nodes_dot = source_coordinates_dot( &
                [first_node, second_node, third_node, fourth_node])
            call generated_fci_cubic_lagrange_weights_jvp( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), nodes(4), &
                target_coordinates_dot(row), nodes_dot(1), nodes_dot(2), &
                nodes_dot(3), nodes_dot(4), weights_dot(1), weights_dot(2), &
                weights_dot(3), weights_dot(4))
            interpolation_map_dot(row, first_node) = weights_dot(1)
            interpolation_map_dot(row, second_node) = weights_dot(2)
            interpolation_map_dot(row, third_node) = weights_dot(3)
            interpolation_map_dot(row, fourth_node) = weights_dot(4)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_cubic_interpolation_map_1d_jvp

    subroutine build_fci_cubic_interpolation_map_1d_vjp( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map_bar, source_coordinates_bar, &
            target_coordinates_bar, status)
        !! Apply the real VJP of the fixed-stencil cubic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(in) :: interpolation_map_bar(:, :)
        real(dp), intent(out) :: source_coordinates_bar(:)
        real(dp), intent(out) :: target_coordinates_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: first_node, second_node, third_node, fourth_node
        real(dp) :: nodes(4), weight_bar(4), node_bar(4), target_bar
        real(dp), allocatable :: interpolation_map(:, :)

        source_coordinates_bar = 0.0_dp
        target_coordinates_bar = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 4 .or. target_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI cubic VJP received incompatible coordinates")
            return
        end if
        allocate(interpolation_map(target_count, source_count))
        call build_fci_cubic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_bar, 1) /= target_count .or. &
            size(interpolation_map_bar, 2) /= source_count .or. &
            size(source_coordinates_bar) /= source_count .or. &
            size(target_coordinates_bar) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI cubic VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(interpolation_map_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI cubic VJP received non-finite cotangents")
            return
        end if

        do row = 1, target_count
            first_node = stencil_indices(1, row)
            second_node = stencil_indices(2, row)
            third_node = stencil_indices(3, row)
            fourth_node = stencil_indices(4, row)
            nodes = source_coordinates([first_node, second_node, third_node, &
                fourth_node])
            weight_bar = [ &
                interpolation_map_bar(row, first_node), &
                interpolation_map_bar(row, second_node), &
                interpolation_map_bar(row, third_node), &
                interpolation_map_bar(row, fourth_node)]
            call generated_fci_cubic_lagrange_weights_vjp( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), nodes(4), &
                weight_bar(1), weight_bar(2), weight_bar(3), weight_bar(4), &
                target_bar, node_bar(1), node_bar(2), node_bar(3), node_bar(4))
            target_coordinates_bar(row) = target_coordinates_bar(row) + &
                target_bar
            source_coordinates_bar(first_node) = &
                source_coordinates_bar(first_node) + node_bar(1)
            source_coordinates_bar(second_node) = &
                source_coordinates_bar(second_node) + node_bar(2)
            source_coordinates_bar(third_node) = &
                source_coordinates_bar(third_node) + node_bar(3)
            source_coordinates_bar(fourth_node) = &
                source_coordinates_bar(fourth_node) + node_bar(4)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_cubic_interpolation_map_1d_vjp

    subroutine build_fci_quartic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        !! Build row-wise quartic Lagrange interpolation weights.
        !!
        !! `stencil_indices(:, row)` supplies five distinct source nodes for
        !! each target.  The local weight kernel is generated by FortSym.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(out) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row
        integer :: node_index(5), node, other
        real(dp) :: nodes(5), weights(5)

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI quartic interpolation received incompatible coordinates")
        interpolation_map = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 5 .or. target_count < 1) return
        if (size(stencil_indices, 1) /= 5 .or. &
            size(stencil_indices, 2) /= target_count) return
        if (size(interpolation_map, 1) /= target_count .or. &
            size(interpolation_map, 2) /= source_count) return
        if (any(.not. ieee_is_finite(source_coordinates)) .or. &
            any(.not. ieee_is_finite(target_coordinates))) return

        do row = 1, target_count
            node_index = stencil_indices(:, row)
            if (any(node_index < 1) .or. any(node_index > source_count)) return
            do node = 1, 5
                do other = node + 1, 5
                    if (node_index(node) == node_index(other)) return
                end do
            end do
            nodes = source_coordinates(node_index)
            do node = 1, 5
                do other = node + 1, 5
                    if (nodes(node) == nodes(other)) return
                end do
            end do
            if (target_coordinates(row) < minval(nodes) .or. &
                target_coordinates(row) > maxval(nodes)) return
            call generated_fci_quartic_lagrange_weights( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), &
                nodes(4), nodes(5), weights(1), weights(2), weights(3), &
                weights(4), weights(5))
            do node = 1, 5
                interpolation_map(row, node_index(node)) = weights(node)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quartic_interpolation_map_1d

    subroutine build_fci_quartic_interpolation_map_1d_jvp( &
            source_coordinates, target_coordinates, stencil_indices, &
            source_coordinates_dot, target_coordinates_dot, &
            interpolation_map_dot, status)
        !! Apply the fixed-stencil JVP of the quartic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(in) :: source_coordinates_dot(:)
        real(dp), intent(in) :: target_coordinates_dot(:)
        real(dp), intent(out) :: interpolation_map_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row, node
        integer :: node_index(5)
        real(dp) :: nodes(5), nodes_dot(5), weights_dot(5)
        real(dp), allocatable :: interpolation_map(:, :)

        interpolation_map_dot = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 5 .or. target_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quartic JVP received incompatible coordinates")
            return
        end if
        allocate(interpolation_map(target_count, source_count))
        call build_fci_quartic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(source_coordinates_dot) /= source_count .or. &
            size(target_coordinates_dot) /= target_count .or. &
            size(interpolation_map_dot, 1) /= target_count .or. &
            size(interpolation_map_dot, 2) /= source_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quartic JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(source_coordinates_dot)) .or. &
            any(.not. ieee_is_finite(target_coordinates_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quartic JVP received non-finite tangents")
            return
        end if

        do row = 1, target_count
            node_index = stencil_indices(:, row)
            nodes = source_coordinates(node_index)
            nodes_dot = source_coordinates_dot(node_index)
            call generated_fci_quartic_lagrange_weights_jvp( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), &
                nodes(4), nodes(5), target_coordinates_dot(row), &
                nodes_dot(1), nodes_dot(2), nodes_dot(3), nodes_dot(4), &
                nodes_dot(5), weights_dot(1), weights_dot(2), weights_dot(3), &
                weights_dot(4), weights_dot(5))
            do node = 1, 5
                interpolation_map_dot(row, node_index(node)) = weights_dot(node)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quartic_interpolation_map_1d_jvp

    subroutine build_fci_quartic_interpolation_map_1d_vjp( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map_bar, source_coordinates_bar, &
            target_coordinates_bar, status)
        !! Apply the real VJP of the fixed-stencil quartic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:)
        integer, intent(in) :: stencil_indices(:, :)
        real(dp), intent(in) :: interpolation_map_bar(:, :)
        real(dp), intent(out) :: source_coordinates_bar(:)
        real(dp), intent(out) :: target_coordinates_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, row, node
        integer :: node_index(5)
        real(dp) :: nodes(5), weight_bar(5), node_bar(5), target_bar
        real(dp), allocatable :: interpolation_map(:, :)

        source_coordinates_bar = 0.0_dp
        target_coordinates_bar = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates)
        if (source_count < 5 .or. target_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quartic VJP received incompatible coordinates")
            return
        end if
        allocate(interpolation_map(target_count, source_count))
        call build_fci_quartic_interpolation_map_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_map, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_bar, 1) /= target_count .or. &
            size(interpolation_map_bar, 2) /= source_count .or. &
            size(source_coordinates_bar) /= source_count .or. &
            size(target_coordinates_bar) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quartic VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(interpolation_map_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI quartic VJP received non-finite cotangents")
            return
        end if

        do row = 1, target_count
            node_index = stencil_indices(:, row)
            nodes = source_coordinates(node_index)
            do node = 1, 5
                weight_bar(node) = interpolation_map_bar(row, node_index(node))
            end do
            call generated_fci_quartic_lagrange_weights_vjp( &
                target_coordinates(row), nodes(1), nodes(2), nodes(3), &
                nodes(4), nodes(5), weight_bar(1), weight_bar(2), &
                weight_bar(3), weight_bar(4), weight_bar(5), target_bar, &
                node_bar(1), node_bar(2), node_bar(3), node_bar(4), node_bar(5))
            target_coordinates_bar(row) = target_coordinates_bar(row) + &
                target_bar
            do node = 1, 5
                source_coordinates_bar(node_index(node)) = &
                    source_coordinates_bar(node_index(node)) + node_bar(node)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quartic_interpolation_map_1d_vjp

    subroutine build_fci_quadratic_interpolation_maps_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_maps, status)
        !! Build quadratic maps for several fixed-topology FCI segments.
        !!
        !! `target_coordinates(:, segment)` and
        !! `stencil_indices(:, :, segment)` describe one mapped slice.  The
        !! output has shape `(n_target, n_source, n_segment)` and uses the
        !! same source-column ordering as the single-slice builder.  Each
        !! segment is delegated to the generated local Lagrange kernel through
        !! `build_fci_quadratic_interpolation_map_1d`.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:, :)
        integer, intent(in) :: stencil_indices(:, :, :)
        real(dp), intent(out) :: interpolation_maps(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, segment_count, segment

        interpolation_maps = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI batched quadratic interpolation received incompatible arrays")
        source_count = size(source_coordinates)
        target_count = size(target_coordinates, 1)
        segment_count = size(target_coordinates, 2)
        if (source_count < 3 .or. target_count < 1 .or. segment_count < 1) return
        if (size(stencil_indices, 1) /= 3 .or. &
            size(stencil_indices, 2) /= target_count .or. &
            size(stencil_indices, 3) /= segment_count) return
        if (size(interpolation_maps, 1) /= target_count .or. &
            size(interpolation_maps, 2) /= source_count .or. &
            size(interpolation_maps, 3) /= segment_count) return
        if (any(.not. ieee_is_finite(source_coordinates)) .or. &
            any(.not. ieee_is_finite(target_coordinates))) return

        do segment = 1, segment_count
            call build_fci_quadratic_interpolation_map_1d( &
                source_coordinates, target_coordinates(:, segment), &
                stencil_indices(:, :, segment), interpolation_maps(:, :, segment), &
                status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quadratic_interpolation_maps_1d

    subroutine build_fci_quadratic_interpolation_maps_1d_jvp( &
            source_coordinates, target_coordinates, stencil_indices, &
            source_coordinates_dot, target_coordinates_dot, &
            interpolation_maps_dot, status)
        !! Apply the fixed-stencil JVP of the batched quadratic map builder.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:, :)
        integer, intent(in) :: stencil_indices(:, :, :)
        real(dp), intent(in) :: source_coordinates_dot(:)
        real(dp), intent(in) :: target_coordinates_dot(:, :)
        real(dp), intent(out) :: interpolation_maps_dot(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, segment_count, segment

        interpolation_maps_dot = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates, 1)
        segment_count = size(target_coordinates, 2)
        if (source_count < 3 .or. target_count < 1 .or. segment_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched quadratic JVP received incompatible coordinates")
            return
        end if
        if (size(stencil_indices, 1) /= 3 .or. &
            size(stencil_indices, 2) /= target_count .or. &
            size(stencil_indices, 3) /= segment_count .or. &
            size(source_coordinates_dot) /= source_count .or. &
            any(shape(target_coordinates_dot) /= shape(target_coordinates)) .or. &
            any(shape(interpolation_maps_dot) /= [target_count, source_count, &
            segment_count])) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched quadratic JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(source_coordinates_dot)) .or. &
            any(.not. ieee_is_finite(target_coordinates_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched quadratic JVP received non-finite tangents")
            return
        end if
        do segment = 1, segment_count
            call build_fci_quadratic_interpolation_map_1d_jvp( &
                source_coordinates, target_coordinates(:, segment), &
                stencil_indices(:, :, segment), source_coordinates_dot, &
                target_coordinates_dot(:, segment), &
                interpolation_maps_dot(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quadratic_interpolation_maps_1d_jvp

    subroutine build_fci_quadratic_interpolation_maps_1d_vjp( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_maps_bar, source_coordinates_bar, &
            target_coordinates_bar, status)
        !! Apply the real VJP of the fixed-stencil batched quadratic map.
        real(dp), intent(in) :: source_coordinates(:)
        real(dp), intent(in) :: target_coordinates(:, :)
        integer, intent(in) :: stencil_indices(:, :, :)
        real(dp), intent(in) :: interpolation_maps_bar(:, :, :)
        real(dp), intent(out) :: source_coordinates_bar(:)
        real(dp), intent(out) :: target_coordinates_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: source_count, target_count, segment_count, segment
        real(dp), allocatable :: interpolation_maps(:, :, :)
        real(dp), allocatable :: source_coordinates_bar_local(:)
        real(dp), allocatable :: target_coordinates_bar_local(:)

        source_coordinates_bar = 0.0_dp
        target_coordinates_bar = 0.0_dp
        source_count = size(source_coordinates)
        target_count = size(target_coordinates, 1)
        segment_count = size(target_coordinates, 2)
        if (source_count < 3 .or. target_count < 1 .or. segment_count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched quadratic VJP received incompatible coordinates")
            return
        end if
        allocate(interpolation_maps(target_count, source_count, segment_count))
        call build_fci_quadratic_interpolation_maps_1d( &
            source_coordinates, target_coordinates, stencil_indices, &
            interpolation_maps, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(interpolation_maps_bar) /= &
            [target_count, source_count, segment_count]) .or. &
            size(source_coordinates_bar) /= source_count .or. &
            any(shape(target_coordinates_bar) /= [target_count, segment_count])) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched quadratic VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(interpolation_maps_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched quadratic VJP received non-finite cotangents")
            return
        end if
        allocate(source_coordinates_bar_local(source_count))
        allocate(target_coordinates_bar_local(target_count))
        do segment = 1, segment_count
            call build_fci_quadratic_interpolation_map_1d_vjp( &
                source_coordinates, target_coordinates(:, segment), &
                stencil_indices(:, :, segment), interpolation_maps_bar(:, :, segment), &
                source_coordinates_bar_local, target_coordinates_bar_local, status)
            if (status%code /= FORTSPARSE_OK) return
            source_coordinates_bar = source_coordinates_bar + &
                source_coordinates_bar_local
            target_coordinates_bar(:, segment) = target_coordinates_bar_local
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_quadratic_interpolation_maps_1d_vjp

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

    subroutine build_fci_triangle_interpolation_map_2d( &
            vertices, triangles, target_points, target_cells, interpolation_map, &
            status)
        !! Build barycentric interpolation weights on a fixed triangle mesh.
        !!
        !! `vertices` has shape `(2, n_vertices)`, `triangles` has shape
        !! `(3, n_triangles)`, and `target_cells(row)` identifies the triangle
        !! containing `target_points(:, row)`.  Connectivity is deliberately
        !! discrete; the derivative routines keep it fixed.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: target_points(:, :)
        integer, intent(in) :: target_cells(:)
        real(dp), intent(out) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_vertices, n_triangles, target_count, row
        integer :: vertex_1, vertex_2, vertex_3
        real(dp) :: determinant, numerator_1, numerator_2
        real(dp) :: weights(3)

        interpolation_map = 0.0_dp
        call validate_triangle_map_inputs( &
            vertices, triangles, target_points, target_cells, interpolation_map, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        n_vertices = size(vertices, 2)
        n_triangles = size(triangles, 2)
        target_count = size(target_points, 2)
        do row = 1, target_count
            call triangle_weights( &
                vertices, triangles, target_points, target_cells, row, vertex_1, &
                vertex_2, vertex_3, determinant, numerator_1, numerator_2, weights)
            if (abs(determinant) <= tiny(determinant)) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI triangle interpolation received a degenerate cell")
                return
            end if
            if (minval(weights) < -1.0e-12_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI triangle interpolation received an outside target")
                return
            end if
            interpolation_map(row, vertex_1) = weights(1)
            interpolation_map(row, vertex_2) = weights(2)
            interpolation_map(row, vertex_3) = weights(3)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_triangle_interpolation_map_2d

    subroutine build_fci_triangle_interpolation_map_2d_jvp( &
            vertices, triangles, target_points, target_cells, vertices_dot, &
            target_points_dot, interpolation_map_dot, status)
        !! Apply the fixed-cell JVP of barycentric interpolation weights.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: target_points(:, :)
        integer, intent(in) :: target_cells(:)
        real(dp), intent(in) :: vertices_dot(:, :)
        real(dp), intent(in) :: target_points_dot(:, :)
        real(dp), intent(out) :: interpolation_map_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_vertices, target_count, row
        integer :: vertex_1, vertex_2, vertex_3
        real(dp) :: determinant, numerator_1, numerator_2, determinant_dot
        real(dp) :: numerator_1_dot, numerator_2_dot
        real(dp) :: weights(3), weights_dot(3)
        real(dp) :: edge_1(2), edge_2(2), edge_1_dot(2), edge_2_dot(2)
        real(dp) :: target_edge_1(2), target_edge_2(2)
        real(dp) :: target_edge_1_dot(2), target_edge_2_dot(2)
        real(dp), allocatable :: interpolation_map(:, :)

        interpolation_map_dot = 0.0_dp
        n_vertices = size(vertices, 2)
        target_count = size(target_points, 2)
        if (size(vertices_dot, 1) /= 2 .or. size(vertices_dot, 2) /= n_vertices .or. &
            size(target_points_dot, 1) /= 2 .or. &
            size(target_points_dot, 2) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle JVP received incompatible tangents")
            return
        end if
        allocate(interpolation_map(target_count, n_vertices))
        call build_fci_triangle_interpolation_map_2d( &
            vertices, triangles, target_points, target_cells, interpolation_map, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_dot, 1) /= target_count .or. &
            size(interpolation_map_dot, 2) /= n_vertices) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle JVP received an incompatible output")
            return
        end if
        if (any(.not. ieee_is_finite(vertices_dot)) .or. &
            any(.not. ieee_is_finite(target_points_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle JVP received non-finite tangents")
            return
        end if

        do row = 1, target_count
            call triangle_weights( &
                vertices, triangles, target_points, target_cells, row, vertex_1, &
                vertex_2, vertex_3, determinant, numerator_1, numerator_2, weights)
            if (minval(weights) <= 1.0e-12_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI triangle JVP crosses a cell-boundary topology event")
                return
            end if
            target_edge_1 = vertices(:, vertex_2) - target_points(:, row)
            target_edge_2 = vertices(:, vertex_3) - target_points(:, row)
            target_edge_1_dot = vertices_dot(:, vertex_2) - &
                target_points_dot(:, row)
            target_edge_2_dot = vertices_dot(:, vertex_3) - &
                target_points_dot(:, row)
            numerator_1_dot = cross_2d_dot( &
                target_edge_1, target_edge_2, target_edge_1_dot, target_edge_2_dot)
            target_edge_1 = vertices(:, vertex_3) - target_points(:, row)
            target_edge_2 = vertices(:, vertex_1) - target_points(:, row)
            target_edge_1_dot = vertices_dot(:, vertex_3) - &
                target_points_dot(:, row)
            target_edge_2_dot = vertices_dot(:, vertex_1) - &
                target_points_dot(:, row)
            numerator_2_dot = cross_2d_dot( &
                target_edge_1, target_edge_2, target_edge_1_dot, target_edge_2_dot)
            edge_1 = vertices(:, vertex_2) - vertices(:, vertex_1)
            edge_2 = vertices(:, vertex_3) - vertices(:, vertex_1)
            edge_1_dot = vertices_dot(:, vertex_2) - vertices_dot(:, vertex_1)
            edge_2_dot = vertices_dot(:, vertex_3) - vertices_dot(:, vertex_1)
            determinant_dot = cross_2d_dot(edge_1, edge_2, edge_1_dot, edge_2_dot)
            weights_dot(1) = (numerator_1_dot*determinant - &
                numerator_1*determinant_dot)/(determinant*determinant)
            weights_dot(2) = (numerator_2_dot*determinant - &
                numerator_2*determinant_dot)/(determinant*determinant)
            weights_dot(3) = -weights_dot(1) - weights_dot(2)
            interpolation_map_dot(row, vertex_1) = weights_dot(1)
            interpolation_map_dot(row, vertex_2) = weights_dot(2)
            interpolation_map_dot(row, vertex_3) = weights_dot(3)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_triangle_interpolation_map_2d_jvp

    subroutine build_fci_triangle_interpolation_map_2d_vjp( &
            vertices, triangles, target_points, target_cells, interpolation_map_bar, &
            vertices_bar, target_points_bar, status)
        !! Apply the fixed-cell VJP of barycentric interpolation weights.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: target_points(:, :)
        integer, intent(in) :: target_cells(:)
        real(dp), intent(in) :: interpolation_map_bar(:, :)
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: target_points_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_vertices, target_count, row
        integer :: vertex_1, vertex_2, vertex_3
        real(dp) :: determinant, numerator_1, numerator_2, determinant_bar
        real(dp) :: numerator_1_bar, numerator_2_bar
        real(dp) :: weights(3), weights_bar(3)
        real(dp) :: edge_1(2), edge_2(2), gradient_1(2), gradient_2(2)
        real(dp), allocatable :: interpolation_map(:, :)

        vertices_bar = 0.0_dp
        target_points_bar = 0.0_dp
        n_vertices = size(vertices, 2)
        target_count = size(target_points, 2)
        allocate(interpolation_map(target_count, n_vertices))
        call build_fci_triangle_interpolation_map_2d( &
            vertices, triangles, target_points, target_cells, interpolation_map, &
            status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(interpolation_map_bar, 1) /= target_count .or. &
            size(interpolation_map_bar, 2) /= n_vertices .or. &
            size(vertices_bar, 1) /= 2 .or. size(vertices_bar, 2) /= n_vertices .or. &
            size(target_points_bar, 1) /= 2 .or. &
            size(target_points_bar, 2) /= target_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle VJP received incompatible cotangents")
            return
        end if
        if (any(.not. ieee_is_finite(interpolation_map_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle VJP received non-finite cotangents")
            return
        end if

        do row = 1, target_count
            call triangle_weights( &
                vertices, triangles, target_points, target_cells, row, vertex_1, &
                vertex_2, vertex_3, determinant, numerator_1, numerator_2, weights)
            if (minval(weights) <= 1.0e-12_dp) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "FCI triangle VJP crosses a cell-boundary topology event")
                return
            end if
            weights_bar = [ &
                interpolation_map_bar(row, vertex_1), &
                interpolation_map_bar(row, vertex_2), &
                interpolation_map_bar(row, vertex_3)]
            weights_bar(1) = weights_bar(1) - weights_bar(3)
            weights_bar(2) = weights_bar(2) - weights_bar(3)
            weights_bar(3) = 0.0_dp
            numerator_1_bar = weights_bar(1)/determinant
            numerator_2_bar = weights_bar(2)/determinant
            determinant_bar = -(weights_bar(1)*numerator_1 + &
                weights_bar(2)*numerator_2)/(determinant*determinant)

            edge_1 = vertices(:, vertex_2) - target_points(:, row)
            edge_2 = vertices(:, vertex_3) - target_points(:, row)
            gradient_1 = [edge_2(2), -edge_2(1)]
            gradient_2 = [-edge_1(2), edge_1(1)]
            vertices_bar(:, vertex_2) = vertices_bar(:, vertex_2) + &
                numerator_1_bar*gradient_1
            vertices_bar(:, vertex_3) = vertices_bar(:, vertex_3) + &
                numerator_1_bar*gradient_2
            target_points_bar(:, row) = target_points_bar(:, row) - &
                numerator_1_bar*(gradient_1 + gradient_2)

            edge_1 = vertices(:, vertex_3) - target_points(:, row)
            edge_2 = vertices(:, vertex_1) - target_points(:, row)
            gradient_1 = [edge_2(2), -edge_2(1)]
            gradient_2 = [-edge_1(2), edge_1(1)]
            vertices_bar(:, vertex_3) = vertices_bar(:, vertex_3) + &
                numerator_2_bar*gradient_1
            vertices_bar(:, vertex_1) = vertices_bar(:, vertex_1) + &
                numerator_2_bar*gradient_2
            target_points_bar(:, row) = target_points_bar(:, row) - &
                numerator_2_bar*(gradient_1 + gradient_2)

            edge_1 = vertices(:, vertex_2) - vertices(:, vertex_1)
            edge_2 = vertices(:, vertex_3) - vertices(:, vertex_1)
            gradient_1 = [edge_2(2), -edge_2(1)]
            gradient_2 = [-edge_1(2), edge_1(1)]
            vertices_bar(:, vertex_2) = vertices_bar(:, vertex_2) + &
                determinant_bar*gradient_1
            vertices_bar(:, vertex_3) = vertices_bar(:, vertex_3) + &
                determinant_bar*gradient_2
            vertices_bar(:, vertex_1) = vertices_bar(:, vertex_1) - &
                determinant_bar*(gradient_1 + gradient_2)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_triangle_interpolation_map_2d_vjp

    subroutine validate_triangle_map_inputs( &
            vertices, triangles, target_points, target_cells, interpolation_map, &
            status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: target_points(:, :)
        integer, intent(in) :: target_cells(:)
        real(dp), intent(in) :: interpolation_map(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_vertices, n_triangles, target_count, row, cell

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI triangle interpolation received incompatible arrays")
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3 .or. &
            size(target_points, 1) /= 2) return
        n_vertices = size(vertices, 2)
        n_triangles = size(triangles, 2)
        target_count = size(target_points, 2)
        if (n_vertices < 3 .or. n_triangles < 1 .or. target_count < 1) return
        if (size(target_cells) /= target_count .or. &
            size(interpolation_map, 1) /= target_count .or. &
            size(interpolation_map, 2) /= n_vertices) return
        if (any(.not. ieee_is_finite(vertices)) .or. &
            any(.not. ieee_is_finite(target_points))) return
        do row = 1, target_count
            cell = target_cells(row)
            if (cell < 1 .or. cell > n_triangles) return
            if (any(triangles(:, cell) < 1) .or. &
                any(triangles(:, cell) > n_vertices)) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_triangle_map_inputs

    subroutine triangle_weights( &
            vertices, triangles, target_points, target_cells, row, vertex_1, &
            vertex_2, vertex_3, determinant, numerator_1, numerator_2, weights)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: target_points(:, :)
        integer, intent(in) :: target_cells(:)
        integer, intent(in) :: row
        integer, intent(out) :: vertex_1, vertex_2, vertex_3
        real(dp), intent(out) :: determinant, numerator_1, numerator_2
        real(dp), intent(out) :: weights(3)

        real(dp) :: edge_1(2), edge_2(2)

        vertex_1 = triangles(1, target_cells(row))
        vertex_2 = triangles(2, target_cells(row))
        vertex_3 = triangles(3, target_cells(row))
        edge_1 = vertices(:, vertex_2) - vertices(:, vertex_1)
        edge_2 = vertices(:, vertex_3) - vertices(:, vertex_1)
        determinant = cross_2d(edge_1, edge_2)
        edge_1 = vertices(:, vertex_2) - target_points(:, row)
        edge_2 = vertices(:, vertex_3) - target_points(:, row)
        numerator_1 = cross_2d(edge_1, edge_2)
        edge_1 = vertices(:, vertex_3) - target_points(:, row)
        edge_2 = vertices(:, vertex_1) - target_points(:, row)
        numerator_2 = cross_2d(edge_1, edge_2)
        weights = [numerator_1/determinant, numerator_2/determinant, &
            1.0_dp - numerator_1/determinant - numerator_2/determinant]
    end subroutine triangle_weights

    pure real(dp) function cross_2d(left, right)
        real(dp), intent(in) :: left(2), right(2)

        cross_2d = left(1)*right(2) - left(2)*right(1)
    end function cross_2d

    pure real(dp) function cross_2d_dot(left, right, left_dot, right_dot)
        real(dp), intent(in) :: left(2), right(2), left_dot(2), right_dot(2)

        cross_2d_dot = left_dot(1)*right(2) + left(1)*right_dot(2) - &
            left_dot(2)*right(1) - left(2)*right_dot(1)
    end function cross_2d_dot

    subroutine build_fci_triangle_interpolation_maps_2d( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, forward_map, backward_map, status)
        !! Build per-segment barycentric maps on a fixed triangle mesh.
        !!
        !! Point arrays have shape `(2, n_staggered, n_segment)`, cell arrays
        !! have shape `(n_staggered, n_segment)`, and map tensors have shape
        !! `(n_staggered, n_vertices, n_segment)`.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: forward_points(:, :, :)
        integer, intent(in) :: forward_cells(:, :)
        real(dp), intent(in) :: backward_points(:, :, :)
        integer, intent(in) :: backward_cells(:, :)
        real(dp), intent(out) :: forward_map(:, :, :)
        real(dp), intent(out) :: backward_map(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: segment

        forward_map = 0.0_dp
        backward_map = 0.0_dp
        call validate_triangle_batch_inputs( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, forward_map, backward_map, status)
        if (status%code /= FORTSPARSE_OK) return
        do segment = 1, size(forward_points, 3)
            call build_fci_triangle_interpolation_map_2d( &
                vertices, triangles, forward_points(:, :, segment), &
                forward_cells(:, segment), forward_map(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            call build_fci_triangle_interpolation_map_2d( &
                vertices, triangles, backward_points(:, :, segment), &
                backward_cells(:, segment), backward_map(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_triangle_interpolation_maps_2d

    subroutine build_fci_triangle_interpolation_maps_2d_jvp( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, vertices_dot, forward_points_dot, backward_points_dot, &
            forward_map_dot, backward_map_dot, status)
        !! Apply the fixed-cell JVP of the batched triangle map builder.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: forward_points(:, :, :)
        integer, intent(in) :: forward_cells(:, :)
        real(dp), intent(in) :: backward_points(:, :, :)
        integer, intent(in) :: backward_cells(:, :)
        real(dp), intent(in) :: vertices_dot(:, :)
        real(dp), intent(in) :: forward_points_dot(:, :, :)
        real(dp), intent(in) :: backward_points_dot(:, :, :)
        real(dp), intent(out) :: forward_map_dot(:, :, :)
        real(dp), intent(out) :: backward_map_dot(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: segment

        forward_map_dot = 0.0_dp
        backward_map_dot = 0.0_dp
        call validate_triangle_batch_inputs( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, forward_map_dot, backward_map_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        if (any(shape(vertices_dot) /= shape(vertices)) .or. &
            any(shape(forward_points_dot) /= shape(forward_points)) .or. &
            any(shape(backward_points_dot) /= shape(backward_points))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle batched JVP received incompatible tangents")
            return
        end if
        if (any(.not. ieee_is_finite(vertices_dot)) .or. &
            any(.not. ieee_is_finite(forward_points_dot)) .or. &
            any(.not. ieee_is_finite(backward_points_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle batched JVP received non-finite tangents")
            return
        end if
        do segment = 1, size(forward_points, 3)
            call build_fci_triangle_interpolation_map_2d_jvp( &
                vertices, triangles, forward_points(:, :, segment), &
                forward_cells(:, segment), vertices_dot, &
                forward_points_dot(:, :, segment), forward_map_dot(:, :, segment), &
                status)
            if (status%code /= FORTSPARSE_OK) return
            call build_fci_triangle_interpolation_map_2d_jvp( &
                vertices, triangles, backward_points(:, :, segment), &
                backward_cells(:, segment), vertices_dot, &
                backward_points_dot(:, :, segment), &
                backward_map_dot(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_triangle_interpolation_maps_2d_jvp

    subroutine build_fci_triangle_interpolation_maps_2d_vjp( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, forward_map_bar, backward_map_bar, vertices_bar, &
            forward_points_bar, backward_points_bar, status)
        !! Apply the fixed-cell VJP of the batched triangle map builder.
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: forward_points(:, :, :)
        integer, intent(in) :: forward_cells(:, :)
        real(dp), intent(in) :: backward_points(:, :, :)
        integer, intent(in) :: backward_cells(:, :)
        real(dp), intent(in) :: forward_map_bar(:, :, :)
        real(dp), intent(in) :: backward_map_bar(:, :, :)
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: forward_points_bar(:, :, :)
        real(dp), intent(out) :: backward_points_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: segment, n_vertices
        real(dp), allocatable :: vertices_bar_local(:, :)

        vertices_bar = 0.0_dp
        forward_points_bar = 0.0_dp
        backward_points_bar = 0.0_dp
        call validate_triangle_batch_inputs( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, forward_map_bar, backward_map_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        n_vertices = size(vertices, 2)
        if (any(shape(vertices_bar) /= shape(vertices))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle batched VJP received an incompatible vertex cotangent")
            return
        end if
        if (any(shape(forward_points_bar) /= shape(forward_points)) .or. &
            any(shape(backward_points_bar) /= shape(backward_points))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI triangle batched VJP received incompatible point cotangents")
            return
        end if
        allocate(vertices_bar_local(2, n_vertices))
        do segment = 1, size(forward_points, 3)
            call build_fci_triangle_interpolation_map_2d_vjp( &
                vertices, triangles, forward_points(:, :, segment), &
                forward_cells(:, segment), forward_map_bar(:, :, segment), &
                vertices_bar_local, forward_points_bar(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            vertices_bar = vertices_bar + vertices_bar_local
            call build_fci_triangle_interpolation_map_2d_vjp( &
                vertices, triangles, backward_points(:, :, segment), &
                backward_cells(:, segment), backward_map_bar(:, :, segment), &
                vertices_bar_local, backward_points_bar(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            vertices_bar = vertices_bar + vertices_bar_local
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_triangle_interpolation_maps_2d_vjp

    subroutine validate_triangle_batch_inputs( &
            vertices, triangles, forward_points, forward_cells, backward_points, &
            backward_cells, forward_map, backward_map, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp), intent(in) :: forward_points(:, :, :)
        integer, intent(in) :: forward_cells(:, :)
        real(dp), intent(in) :: backward_points(:, :, :)
        integer, intent(in) :: backward_cells(:, :)
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: n_vertices, n_staggered, n_segment

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI batched triangle map received incompatible arrays")
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3 .or. &
            size(forward_points, 1) /= 2 .or. size(backward_points, 1) /= 2) return
        n_vertices = size(vertices, 2)
        n_staggered = size(forward_points, 2)
        n_segment = size(forward_points, 3)
        if (n_vertices < 3 .or. n_staggered < 1 .or. n_segment < 1) return
        if (any(shape(backward_points) /= shape(forward_points)) .or. &
            any(shape(forward_cells) /= [n_staggered, n_segment]) .or. &
            any(shape(backward_cells) /= shape(forward_cells))) return
        if (any(shape(forward_map) /= [n_staggered, n_vertices, n_segment]) .or. &
            any(shape(backward_map) /= shape(forward_map))) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_triangle_batch_inputs

    subroutine build_fci_bilinear_interpolation_maps_2d( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            forward_map, backward_map, status)
        !! Build per-segment bilinear maps from traced FCI endpoints.
        !!
        !! Endpoint arrays use shape `(n_staggered, n_segment)`.  The map
        !! tensors use shape `(n_staggered, size(source_x)*size(source_y),
        !! n_segment)`, with the same x-fast source-column ordering as the
        !! single-slice builder.  Every segment is checked independently, so
        !! an out-of-domain endpoint cannot silently contaminate the operator.
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: forward_x(:, :)
        real(dp), intent(in) :: forward_y(:, :)
        real(dp), intent(in) :: backward_x(:, :)
        real(dp), intent(in) :: backward_y(:, :)
        real(dp), intent(out) :: forward_map(:, :, :)
        real(dp), intent(out) :: backward_map(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: segment, n_segment

        forward_map = 0.0_dp
        backward_map = 0.0_dp
        call validate_fci_bilinear_batch_shapes( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            forward_map, backward_map, status)
        if (status%code /= FORTSPARSE_OK) return
        n_segment = size(forward_x, 2)
        do segment = 1, n_segment
            call build_fci_bilinear_interpolation_map_2d( &
                source_x, source_y, forward_x(:, segment), forward_y(:, segment), &
                forward_map(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            call build_fci_bilinear_interpolation_map_2d( &
                source_x, source_y, backward_x(:, segment), &
                backward_y(:, segment), backward_map(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_bilinear_interpolation_maps_2d

    subroutine build_fci_bilinear_interpolation_maps_2d_jvp( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            source_x_dot, source_y_dot, forward_x_dot, forward_y_dot, &
            backward_x_dot, backward_y_dot, forward_map_dot, backward_map_dot, &
            status)
        !! Apply the fixed-topology JVP of the batched FCI map builder.
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: forward_x(:, :)
        real(dp), intent(in) :: forward_y(:, :)
        real(dp), intent(in) :: backward_x(:, :)
        real(dp), intent(in) :: backward_y(:, :)
        real(dp), intent(in) :: source_x_dot(:)
        real(dp), intent(in) :: source_y_dot(:)
        real(dp), intent(in) :: forward_x_dot(:, :)
        real(dp), intent(in) :: forward_y_dot(:, :)
        real(dp), intent(in) :: backward_x_dot(:, :)
        real(dp), intent(in) :: backward_y_dot(:, :)
        real(dp), intent(out) :: forward_map_dot(:, :, :)
        real(dp), intent(out) :: backward_map_dot(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: segment, nx, ny, n_segment

        forward_map_dot = 0.0_dp
        backward_map_dot = 0.0_dp
        call validate_fci_bilinear_batch_shapes( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            forward_map_dot, backward_map_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        nx = size(source_x)
        ny = size(source_y)
        n_segment = size(forward_x, 2)
        if (size(source_x_dot) /= nx .or. size(source_y_dot) /= ny .or. &
            any(shape(forward_x_dot) /= shape(forward_x)) .or. &
            any(shape(forward_y_dot) /= shape(forward_y)) .or. &
            any(shape(backward_x_dot) /= shape(backward_x)) .or. &
            any(shape(backward_y_dot) /= shape(backward_y))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched-map JVP received incompatible tangents")
            return
        end if
        do segment = 1, n_segment
            call build_fci_bilinear_interpolation_map_2d_jvp( &
                source_x, source_y, forward_x(:, segment), forward_y(:, segment), &
                source_x_dot, source_y_dot, forward_x_dot(:, segment), &
                forward_y_dot(:, segment), forward_map_dot(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            call build_fci_bilinear_interpolation_map_2d_jvp( &
                source_x, source_y, backward_x(:, segment), &
                backward_y(:, segment), source_x_dot, source_y_dot, &
                backward_x_dot(:, segment), backward_y_dot(:, segment), &
                backward_map_dot(:, :, segment), status)
            if (status%code /= FORTSPARSE_OK) return
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_bilinear_interpolation_maps_2d_jvp

    subroutine build_fci_bilinear_interpolation_maps_2d_vjp( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            forward_map_bar, backward_map_bar, source_x_bar, source_y_bar, &
            forward_x_bar, forward_y_bar, backward_x_bar, backward_y_bar, status)
        !! Apply the fixed-topology VJP of the batched FCI map builder.
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: forward_x(:, :)
        real(dp), intent(in) :: forward_y(:, :)
        real(dp), intent(in) :: backward_x(:, :)
        real(dp), intent(in) :: backward_y(:, :)
        real(dp), intent(in) :: forward_map_bar(:, :, :)
        real(dp), intent(in) :: backward_map_bar(:, :, :)
        real(dp), intent(out) :: source_x_bar(:)
        real(dp), intent(out) :: source_y_bar(:)
        real(dp), intent(out) :: forward_x_bar(:, :)
        real(dp), intent(out) :: forward_y_bar(:, :)
        real(dp), intent(out) :: backward_x_bar(:, :)
        real(dp), intent(out) :: backward_y_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: segment, nx, ny, n_segment
        real(dp), allocatable :: source_x_bar_local(:), source_y_bar_local(:)

        source_x_bar = 0.0_dp
        source_y_bar = 0.0_dp
        forward_x_bar = 0.0_dp
        forward_y_bar = 0.0_dp
        backward_x_bar = 0.0_dp
        backward_y_bar = 0.0_dp
        call validate_fci_bilinear_batch_shapes( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            forward_map_bar, backward_map_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        nx = size(source_x)
        ny = size(source_y)
        n_segment = size(forward_x, 2)
        if (size(source_x_bar) /= nx .or. size(source_y_bar) /= ny .or. &
            any(shape(forward_x_bar) /= shape(forward_x)) .or. &
            any(shape(forward_y_bar) /= shape(forward_y)) .or. &
            any(shape(backward_x_bar) /= shape(backward_x)) .or. &
            any(shape(backward_y_bar) /= shape(backward_y))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "FCI batched-map VJP received incompatible cotangents")
            return
        end if
        allocate( &
            source_x_bar_local(nx), source_y_bar_local(ny))
        do segment = 1, n_segment
            call build_fci_bilinear_interpolation_map_2d_vjp( &
                source_x, source_y, forward_x(:, segment), forward_y(:, segment), &
                forward_map_bar(:, :, segment), source_x_bar_local, &
                source_y_bar_local, forward_x_bar(:, segment), &
                forward_y_bar(:, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            source_x_bar = source_x_bar + source_x_bar_local
            source_y_bar = source_y_bar + source_y_bar_local
            call build_fci_bilinear_interpolation_map_2d_vjp( &
                source_x, source_y, backward_x(:, segment), &
                backward_y(:, segment), backward_map_bar(:, :, segment), &
                source_x_bar_local, source_y_bar_local, &
                backward_x_bar(:, segment), backward_y_bar(:, segment), status)
            if (status%code /= FORTSPARSE_OK) return
            source_x_bar = source_x_bar + source_x_bar_local
            source_y_bar = source_y_bar + source_y_bar_local
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fci_bilinear_interpolation_maps_2d_vjp

    subroutine validate_fci_bilinear_batch_shapes( &
            source_x, source_y, forward_x, forward_y, backward_x, backward_y, &
            forward_map, backward_map, status)
        real(dp), intent(in) :: source_x(:)
        real(dp), intent(in) :: source_y(:)
        real(dp), intent(in) :: forward_x(:, :)
        real(dp), intent(in) :: forward_y(:, :)
        real(dp), intent(in) :: backward_x(:, :)
        real(dp), intent(in) :: backward_y(:, :)
        real(dp), intent(in) :: forward_map(:, :, :)
        real(dp), intent(in) :: backward_map(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: nx, ny, n_staggered, n_segment

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "FCI batched bilinear map received incompatible arrays")
        nx = size(source_x)
        ny = size(source_y)
        n_staggered = size(forward_x, 1)
        n_segment = size(forward_x, 2)
        if (nx < 2 .or. ny < 2 .or. n_staggered < 1 .or. n_segment < 1) return
        if (any(shape(forward_y) /= shape(forward_x)) .or. &
            any(shape(backward_x) /= shape(forward_x)) .or. &
            any(shape(backward_y) /= shape(forward_x))) return
        if (any(shape(forward_map) /= [n_staggered, nx*ny, n_segment])) return
        if (any(shape(backward_map) /= [n_staggered, nx*ny, n_segment])) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_fci_bilinear_batch_shapes

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
