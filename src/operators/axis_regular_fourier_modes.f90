module fortfem_axis_regular_fourier_modes
    !! Scalar Fourier/radial metadata for a regular polar axis.
    !!
    !! A smooth scalar coefficient with poloidal label m has a radial Taylor
    !! series containing only powers p >= |m| with p = |m| (mod 2).  This
    !! contract is intentionally only a basis metadata validator: vector
    !! components may use shifted effective labels, and all equilibrium and
    !! profile choices remain caller-owned.
    use fortfem_fourier_mode_registry, only: &
        find_fourier_mode, fourier_mode_registry_t, &
        validate_fourier_mode_registry
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: AXIS_RADIAL_PARITY_EVEN = 0
    integer, parameter, public :: AXIS_RADIAL_PARITY_ODD = 1

    type, public :: axis_regular_mode_record_t
        integer :: registry_index = 0
        integer :: poloidal_mode = 0
        integer :: toroidal_mode = 0
        integer :: radial_power = 0
        integer :: required_minimum_power = 0
        integer :: required_parity = AXIS_RADIAL_PARITY_EVEN
        integer :: conjugate_registry_index = 0
        integer :: conjugate_order_index = 0
        logical :: axis_regular = .false.
    end type axis_regular_mode_record_t

    type, public :: axis_regular_mode_table_t
        type(axis_regular_mode_record_t), allocatable :: modes(:)
        integer :: field_periods = 1
        logical :: real_packed = .false.
    end type axis_regular_mode_table_t

    public :: axis_regular_mode_requirements
    public :: build_axis_regular_mode_table
    public :: validate_axis_regular_mode_table

contains

    subroutine axis_regular_mode_requirements( &
            poloidal_mode, radial_power, required_minimum_power, &
            required_parity, axis_regular, status)
        !! Report and validate the scalar axis contract for one mode.
        integer, intent(in) :: poloidal_mode, radial_power
        integer, intent(out) :: required_minimum_power, required_parity
        logical, intent(out) :: axis_regular
        type(fortsparse_status_t), intent(out) :: status

        required_minimum_power = abs(poloidal_mode)
        required_parity = modulo(required_minimum_power, 2)
        axis_regular = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "axis-regular Fourier mode has inconsistent radial metadata")
        if (radial_power < 0) return
        axis_regular = radial_power >= required_minimum_power .and. &
            modulo(radial_power, 2) == required_parity
        if (.not. axis_regular) return
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine axis_regular_mode_requirements

    subroutine build_axis_regular_mode_table(registry, table, status)
        !! Canonically order a registry and validate scalar axis metadata.
        !!
        !! Conjugate pairs are adjacent in the result whenever both modes are
        !! retained.  `conjugate_order_index` is zero for a missing pair in a
        !! complex (non-real-packed) registry and points to the self/pair row
        !! otherwise.
        type(fourier_mode_registry_t), intent(in) :: registry
        type(axis_regular_mode_table_t), intent(out) :: table
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: order(:)
        integer :: mode_count, mode, source, minimum_power, parity
        logical :: regular

        call clear_axis_regular_mode_table(table)
        if (.not. validate_fourier_mode_registry(registry, status)) return
        mode_count = size(registry%poloidal_modes)
        allocate(order(mode_count))
        do mode = 1, mode_count
            order(mode) = mode
        end do
        call sort_mode_order(registry, order)
        allocate(table%modes(mode_count))
        table%field_periods = registry%field_periods
        table%real_packed = registry%real_packed
        do mode = 1, mode_count
            source = order(mode)
            call axis_regular_mode_requirements( &
                registry%poloidal_modes(source), registry%radial_powers(source), &
                minimum_power, parity, regular, status)
            table%modes(mode)%registry_index = source
            table%modes(mode)%poloidal_mode = registry%poloidal_modes(source)
            table%modes(mode)%toroidal_mode = registry%toroidal_modes(source)
            table%modes(mode)%radial_power = registry%radial_powers(source)
            table%modes(mode)%required_minimum_power = minimum_power
            table%modes(mode)%required_parity = parity
            table%modes(mode)%axis_regular = regular
            table%modes(mode)%conjugate_registry_index = find_fourier_mode( &
                registry, -registry%poloidal_modes(source), &
                -registry%toroidal_modes(source))
            if (status%code /= FORTSPARSE_OK) then
                deallocate(order)
                call clear_axis_regular_mode_table(table)
                return
            end if
        end do
        deallocate(order)
        do mode = 1, mode_count
            do source = 1, mode_count
                if (table%modes(source)%registry_index /= &
                        table%modes(mode)%conjugate_registry_index) cycle
                table%modes(mode)%conjugate_order_index = source
                exit
            end do
        end do
        if (.not. validate_axis_regular_mode_table(table, status)) then
            call clear_axis_regular_mode_table(table)
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_axis_regular_mode_table

    logical function validate_axis_regular_mode_table(table, status) &
            result(valid)
        type(axis_regular_mode_table_t), intent(in) :: table
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, other, minimum_power, parity
        logical :: regular

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "axis-regular Fourier mode table is invalid")
        if (.not. allocated(table%modes)) return
        if (size(table%modes) < 1 .or. table%field_periods < 1) return
        do mode = 1, size(table%modes)
            if (table%modes(mode)%registry_index < 1) return
            call axis_regular_mode_requirements( &
                table%modes(mode)%poloidal_mode, table%modes(mode)%radial_power, &
                minimum_power, parity, regular, status)
            if (.not. regular .or. status%code /= FORTSPARSE_OK) return
            if (table%modes(mode)%required_minimum_power /= minimum_power .or. &
                    table%modes(mode)%required_parity /= parity .or. &
                    .not. table%modes(mode)%axis_regular) return
            if (mode > 1) then
                if (.not. mode_precedes(table%modes(mode - 1), table%modes(mode))) &
                    return
            end if
            do other = 1, mode - 1
                if (table%modes(other)%registry_index == &
                        table%modes(mode)%registry_index) return
            end do
            if (table%modes(mode)%conjugate_order_index < 0 .or. &
                    table%modes(mode)%conjugate_order_index > size(table%modes)) return
            if (table%modes(mode)%conjugate_order_index > 0) then
                other = table%modes(mode)%conjugate_order_index
                if (table%modes(other)%registry_index /= &
                        table%modes(mode)%conjugate_registry_index) return
                if (table%modes(other)%conjugate_order_index /= mode) return
                if (table%real_packed .and. table%modes(other)%radial_power /= &
                        table%modes(mode)%radial_power) return
            else if (table%real_packed) then
                return
            end if
        end do
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_axis_regular_mode_table

    subroutine sort_mode_order(registry, order)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(inout) :: order(:)

        integer :: position, previous, candidate

        do position = 2, size(order)
            candidate = order(position)
            previous = position - 1
            do while (previous >= 1)
                if (.not. registry_mode_precedes( &
                        registry, candidate, order(previous))) exit
                order(previous + 1) = order(previous)
                previous = previous - 1
            end do
            order(previous + 1) = candidate
        end do
    end subroutine sort_mode_order

    logical function registry_mode_precedes(registry, left, right) result(precedes)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: left, right
        integer :: left_m, left_n, left_orientation
        integer :: right_m, right_n, right_orientation

        call canonical_pair_key( &
            registry%poloidal_modes(left), registry%toroidal_modes(left), &
            left_m, left_n, left_orientation)
        call canonical_pair_key( &
            registry%poloidal_modes(right), registry%toroidal_modes(right), &
            right_m, right_n, right_orientation)
        precedes = left_m < right_m
        if (left_m == right_m) precedes = left_n < right_n
        if (left_m == right_m .and. left_n == right_n) then
            precedes = left_orientation < right_orientation
        end if
        if (left_m == right_m .and. left_n == right_n .and. &
                left_orientation == right_orientation) precedes = left < right
    end function registry_mode_precedes

    logical function mode_precedes(left, right) result(precedes)
        type(axis_regular_mode_record_t), intent(in) :: left, right
        integer :: left_m, left_n, left_orientation
        integer :: right_m, right_n, right_orientation

        call canonical_pair_key( &
            left%poloidal_mode, left%toroidal_mode, &
            left_m, left_n, left_orientation)
        call canonical_pair_key( &
            right%poloidal_mode, right%toroidal_mode, &
            right_m, right_n, right_orientation)
        precedes = left_m < right_m
        if (left_m == right_m) precedes = left_n < right_n
        if (left_m == right_m .and. left_n == right_n) then
            precedes = left_orientation < right_orientation
        end if
        if (left_m == right_m .and. left_n == right_n .and. &
                left_orientation == right_orientation) then
            precedes = left%registry_index < right%registry_index
        end if
    end function mode_precedes

    subroutine canonical_pair_key( &
            poloidal_mode, toroidal_mode, key_m, key_n, orientation)
        integer, intent(in) :: poloidal_mode, toroidal_mode
        integer, intent(out) :: key_m, key_n, orientation

        if (poloidal_mode > 0 .or. &
                (poloidal_mode == 0 .and. toroidal_mode >= 0)) then
            key_m = poloidal_mode
            key_n = toroidal_mode
            orientation = 0
        else
            key_m = -poloidal_mode
            key_n = -toroidal_mode
            orientation = 1
        end if
    end subroutine canonical_pair_key

    subroutine clear_axis_regular_mode_table(table)
        type(axis_regular_mode_table_t), intent(out) :: table

        if (allocated(table%modes)) deallocate(table%modes)
        table%field_periods = 1
        table%real_packed = .false.
    end subroutine clear_axis_regular_mode_table

end module fortfem_axis_regular_fourier_modes
