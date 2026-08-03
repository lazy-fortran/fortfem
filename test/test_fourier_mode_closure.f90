program test_fourier_mode_closure
    use check, only: check_condition, check_summary
    use fortfem_fourier, only: &
        build_fourier_mode_closure_registry, &
        build_fourier_mode_padded_registry, &
        find_fourier_mode, fourier_mode_registry_t, &
        initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: input_mode_count = 3
    integer, parameter :: input_poloidal(input_mode_count) = [-1, 0, 1]
    integer, parameter :: input_toroidal(input_mode_count) = [0, 0, 0]
    integer, parameter :: field_periods = 4
    real(dp), parameter :: poloidal_phase = 0.21_dp
    real(dp), parameter :: toroidal_phase = -0.31_dp
    integer, parameter :: radial_powers(input_mode_count) = [1, 0, 1]
    type(fourier_mode_registry_t) :: input_registry, padded_registry
    type(fourier_mode_registry_t) :: closure_registry, one_round_registry
    type(fortsparse_status_t) :: status
    integer :: first_mode, second_mode, closure_mode
    logical :: pair_closure

    call initialize_fourier_mode_registry( &
        input_registry, input_poloidal, input_toroidal, field_periods, &
        poloidal_phase, toroidal_phase, .false., &
        radial_powers=radial_powers, status=status)
    call check_condition(status%code == 0, &
        "repeated Fourier closure accepts the input registry")
    call build_fourier_mode_padded_registry( &
        input_registry, padded_registry, status)
    if (status%code /= 0) error stop 1
    call build_fourier_mode_closure_registry( &
        input_registry, 1, one_round_registry, status)
    if (status%code /= 0) error stop 1
    call check_condition(size(one_round_registry%poloidal_modes) == &
        size(padded_registry%poloidal_modes) .and. &
        all(one_round_registry%poloidal_modes == padded_registry%poloidal_modes) .and. &
        all(one_round_registry%toroidal_modes == padded_registry%toroidal_modes), &
        "one closure round equals the padded Fourier work registry")

    call build_fourier_mode_closure_registry( &
        input_registry, 2, closure_registry, status)
    call check_condition(status%code == 0 .and. &
        size(closure_registry%poloidal_modes) > &
        size(padded_registry%poloidal_modes), &
        "two closure rounds expand the Fourier work registry")
    pair_closure = .true.
    do first_mode = 1, size(padded_registry%poloidal_modes)
        do second_mode = 1, size(padded_registry%poloidal_modes)
            closure_mode = find_fourier_mode(closure_registry, &
                padded_registry%poloidal_modes(first_mode) + &
                padded_registry%poloidal_modes(second_mode), &
                padded_registry%toroidal_modes(first_mode) + &
                padded_registry%toroidal_modes(second_mode))
            if (closure_mode == 0) pair_closure = .false.
        end do
    end do
    call check_condition(pair_closure, &
        "repeated Fourier closure retains every prior-work-set pair sum")
    call check_condition(closure_registry%field_periods == field_periods .and. &
        abs(closure_registry%poloidal_phase - poloidal_phase) < 1.0e-14_dp .and. &
        abs(closure_registry%toroidal_phase - toroidal_phase) < 1.0e-14_dp, &
        "repeated Fourier closure preserves phase and field-period metadata")

    call build_fourier_mode_closure_registry( &
        input_registry, 0, closure_registry, status)
    call check_condition(status%code /= 0, &
        "repeated Fourier closure rejects a non-positive round count")
    call check_summary("repeated Fourier mode closure")
end program test_fourier_mode_closure
