program test_axis_regular_fourier_modes
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        AXIS_RADIAL_PARITY_EVEN, AXIS_RADIAL_PARITY_ODD, &
        axis_regular_mode_requirements, axis_regular_mode_table_t, &
        build_axis_regular_mode_table, fourier_mode_registry_t, &
        initialize_fourier_mode_registry, validate_axis_regular_mode_table
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(fourier_mode_registry_t) :: registry, bad_registry
    type(fourier_mode_registry_t) :: reordered, complex_registry
    type(axis_regular_mode_table_t) :: table, reordered_table, complex_table
    type(fortsparse_status_t) :: status
    integer :: minimum_power, parity
    logical :: regular
    logical :: all_passed

    all_passed = .true.
    call initialize_fourier_mode_registry( &
        registry, [2, -1, 0, -2, 1], [0, 1, 0, 0, -1], 3, &
        0.1_dp, -0.2_dp, .true., radial_powers=[2, 1, 0, 2, 1], &
        status=status)
    call record(status%code == 0, &
        "real-packed registry accepts a scalar axis-regular mode set")
    call build_axis_regular_mode_table(registry, table, status)
    call record(status%code == 0 .and. &
        validate_axis_regular_mode_table(table, status), &
        "axis-regular mode table validates the registry")
    call record(table%field_periods == 3 .and. table%real_packed, &
        "axis-regular table preserves field-period and packing metadata")
    call record(size(table%modes) == 5 .and. &
        table%modes(1)%poloidal_mode == 0 .and. &
        table%modes(2)%poloidal_mode == 1 .and. &
        table%modes(3)%poloidal_mode == -1 .and. &
        table%modes(4)%poloidal_mode == 2 .and. &
        table%modes(5)%poloidal_mode == -2, &
        "axis-regular table puts each conjugate pair in deterministic order")
    call record(table%modes(2)%conjugate_order_index == 3 .and. &
        table%modes(3)%conjugate_order_index == 2 .and. &
        table%modes(4)%conjugate_order_index == 5 .and. &
        table%modes(5)%conjugate_order_index == 4, &
        "axis-regular table records conjugate-safe order indices")
    call record(table%modes(1)%required_minimum_power == 0 .and. &
        table%modes(1)%required_parity == AXIS_RADIAL_PARITY_EVEN .and. &
        table%modes(1)%axis_regular, &
        "axisymmetric mode requires an even radial power starting at zero")
    call record(table%modes(2)%required_minimum_power == 1 .and. &
        table%modes(2)%required_parity == AXIS_RADIAL_PARITY_ODD .and. &
        table%modes(2)%axis_regular, &
        "non-axisymmetric mode requires its odd scalar radial parity")

    bad_registry = registry
    bad_registry%radial_powers(1) = 3
    call build_axis_regular_mode_table(bad_registry, reordered_table, status)
    call record(status%code /= 0 .and. .not. allocated(reordered_table%modes), &
        "axis-regular table rejects inconsistent registry radial metadata")

    call axis_regular_mode_requirements(2, 4, minimum_power, parity, regular, status)
    call record(status%code == 0 .and. minimum_power == 2 .and. &
        parity == AXIS_RADIAL_PARITY_EVEN .and. regular, &
        "axis metadata reports a valid even power for m=2")
    call axis_regular_mode_requirements(2, 3, minimum_power, parity, regular, status)
    call record(status%code /= 0 .and. minimum_power == 2 .and. &
        parity == AXIS_RADIAL_PARITY_EVEN .and. .not. regular, &
        "axis metadata rejects an inconsistent parity")
    call axis_regular_mode_requirements(-1, 0, minimum_power, parity, regular, status)
    call record(status%code /= 0 .and. minimum_power == 1 .and. &
        parity == AXIS_RADIAL_PARITY_ODD .and. .not. regular, &
        "axis metadata rejects a power below |m|")

    call initialize_fourier_mode_registry( &
        reordered, [-2, 1, 0, 2, -1], [0, -1, 0, 0, 1], 3, &
        0.1_dp, -0.2_dp, .true., radial_powers=[2, 1, 0, 2, 1], &
        status=status)
    call build_axis_regular_mode_table(reordered, reordered_table, status)
    call record(status%code == 0 .and. size(reordered_table%modes) == &
        size(table%modes), "reordered input still produces a complete table")
    call record(all(reordered_table%modes%poloidal_mode == &
        table%modes%poloidal_mode) .and. &
        all(reordered_table%modes%toroidal_mode == table%modes%toroidal_mode) .and. &
        all(reordered_table%modes%radial_power == table%modes%radial_power), &
        "axis-regular ordering is independent of registry input order")

    call initialize_fourier_mode_registry( &
        complex_registry, [1], [2], 2, 0.0_dp, 0.0_dp, .false., &
        radial_powers=[1], status=status)
    call build_axis_regular_mode_table(complex_registry, complex_table, status)
    call record(status%code == 0 .and. &
        complex_table%modes(1)%conjugate_order_index == 0, &
        "complex mode tables permit an intentionally omitted conjugate")

    call check_summary("axis-regular Fourier metadata")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_axis_regular_fourier_modes
