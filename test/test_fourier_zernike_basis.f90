program test_fourier_zernike_basis
    use check, only: check_condition, check_summary
    use fortfem_fourier_zernike_basis, only: &
        FOURIER_ZERNIKE_PARITY_EVEN, FOURIER_ZERNIKE_PARITY_ODD, &
        build_fourier_zernike_basis, evaluate_fourier_zernike_radial, &
        fourier_zernike_basis_t, fourier_zernike_mode_requirements, &
        validate_fourier_zernike_basis
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, FORTSPARSE_OK
    implicit none

    integer, parameter :: sample_count = 4
    real(dp), parameter :: rho(sample_count) = [ &
        0.0_dp, 0.25_dp, 0.7_dp, 1.0_dp]
    real(dp), parameter :: tolerance = 2.0e-13_dp

    type(fourier_zernike_basis_t) :: basis, repeat_basis, invalid_basis
    type(fortsparse_status_t) :: status
    real(dp) :: values(sample_count), expected(sample_count)
    real(dp) :: short_values(sample_count - 1)
    integer :: mode, pair, minimum_power, parity
    logical :: admissible, all_passed

    all_passed = .true.

    call build_fourier_zernike_basis(3, 2, 1, basis, status)
    call record(status%code == FORTSPARSE_OK .and. &
        validate_fourier_zernike_basis(basis, status), &
        "Fourier-Zernike mode table validates")
    call record(size(basis%modes) == 24, &
        "mode enumeration keeps only admissible l,m pairs")

    call record(basis%modes(1)%l == 0 .and. basis%modes(1)%m == 0 .and. &
        basis%modes(1)%n == 0 .and. basis%modes(2)%n == 1 .and. &
        basis%modes(3)%n == -1, &
        "axisymmetric zero and conjugate toroidal modes have stable order")
    call record(basis%modes(4)%l == 1 .and. basis%modes(4)%m == 1 .and. &
        basis%modes(4)%n == -1 .and. basis%modes(5)%m == -1 .and. &
        basis%modes(5)%n == 1, &
        "non-axisymmetric conjugate representatives are adjacent")

    do mode = 1, size(basis%modes)
        call record(basis%modes(mode)%minimum_power == &
            abs(basis%modes(mode)%m), "mode stores its minimum radial power")
        call record(basis%modes(mode)%parity == modulo( &
            abs(basis%modes(mode)%m), 2), "mode stores scalar radial parity")
        pair = basis%modes(mode)%conjugate_index
        call record(pair >= 1 .and. pair <= size(basis%modes), &
            "every mode has a conjugate index")
        if (pair >= 1 .and. pair <= size(basis%modes)) then
            call record(basis%modes(pair)%l == basis%modes(mode)%l .and. &
                basis%modes(pair)%m == -basis%modes(mode)%m .and. &
                basis%modes(pair)%n == -basis%modes(mode)%n, &
                "conjugate index negates both Fourier labels")
            if (pair /= mode) call record(abs(pair - mode) == 1, &
                "conjugate pairs are contiguous")
        end if
    end do

    call build_fourier_zernike_basis(3, 2, 1, repeat_basis, status)
    call record(status%code == FORTSPARSE_OK .and. &
        same_mode_labels(basis, repeat_basis), &
        "repeated enumeration is deterministic")

    expected = 2.0_dp*rho**2 - 1.0_dp
    call evaluate_fourier_zernike_radial(2, 0, rho, values, status)
    call record(status%code == FORTSPARSE_OK .and. &
        maxval(abs(values - expected)) < tolerance, &
        "R_2^0 agrees with the independent polynomial oracle")

    expected = rho**2
    call evaluate_fourier_zernike_radial(2, 2, rho, values, status)
    call record(status%code == FORTSPARSE_OK .and. &
        maxval(abs(values - expected)) < tolerance, &
        "R_2^2 agrees with the independent polynomial oracle")

    expected = 3.0_dp*rho**3 - 2.0_dp*rho
    call evaluate_fourier_zernike_radial(3, -1, rho, values, status)
    call record(status%code == FORTSPARSE_OK .and. &
        maxval(abs(values - expected)) < tolerance, &
        "negative m uses the scalar |m| Zernike radial polynomial")

    call fourier_zernike_mode_requirements(4, -2, minimum_power, parity, &
        admissible, status)
    call record(status%code == FORTSPARSE_OK .and. admissible .and. &
        minimum_power == 2 .and. parity == FOURIER_ZERNIKE_PARITY_EVEN, &
        "valid labels report minimum power and even parity")
    call fourier_zernike_mode_requirements(1, 0, minimum_power, parity, &
        admissible, status)
    call record(status%code /= FORTSPARSE_OK .and. .not. admissible .and. &
        minimum_power == 0 .and. parity == FOURIER_ZERNIKE_PARITY_EVEN, &
        "invalid l,m parity is rejected")

    call evaluate_fourier_zernike_radial(2, 1, rho, values, status)
    call record(status%code /= FORTSPARSE_OK .and. all(values == 0.0_dp), &
        "radial evaluation rejects a non-admissible l,m pair")
    call evaluate_fourier_zernike_radial(2, 0, [0.0_dp, 1.1_dp], &
        values(:2), status)
    call record(status%code /= FORTSPARSE_OK .and. all(values(:2) == 0.0_dp), &
        "radial evaluation rejects samples outside the unit interval")
    call evaluate_fourier_zernike_radial(2, 0, rho, short_values, status)
    call record(status%code /= FORTSPARSE_OK, &
        "radial evaluation rejects an incompatible output shape")

    call build_fourier_zernike_basis(-1, 2, 1, invalid_basis, status)
    call record(status%code /= FORTSPARSE_OK .and. &
        .not. allocated(invalid_basis%modes), &
        "mode enumeration rejects negative bounds")

    call check_summary("Fourier-Zernike basis")
    if (.not. all_passed) error stop 1

contains

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

    logical function same_mode_labels(left, right) result(equal)
        type(fourier_zernike_basis_t), intent(in) :: left, right
        integer :: index

        equal = .false.
        if (.not. allocated(left%modes) .or. .not. allocated(right%modes)) return
        if (size(left%modes) /= size(right%modes)) return
        do index = 1, size(left%modes)
            if (left%modes(index)%l /= right%modes(index)%l .or. &
                    left%modes(index)%m /= right%modes(index)%m .or. &
                    left%modes(index)%n /= right%modes(index)%n .or. &
                    left%modes(index)%conjugate_index /= &
                    right%modes(index)%conjugate_index) return
        end do
        equal = .true.
    end function same_mode_labels

end program test_fourier_zernike_basis
