module fortfem_toroidal_spectral_metadata
    !! Fixed-topology truncation and zero-mode diagnostics for toroidal traces.
    !!
    !! The modal window is a caller-owned rectangular policy:
    !! degree <= degree_limit and order <= order_limit are retained.  The
    !! routines report exactly what is present in the supplied list; they do
    !! not estimate an unprovided spectral tail or choose a Green kernel.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: analyze_toroidal_spectral_modes
    public :: analyze_toroidal_spectral_modes_jvp
    public :: analyze_toroidal_spectral_modes_vjp

contains

    subroutine analyze_toroidal_spectral_modes( &
            degree_indices, orders, coefficients, degree_limit, order_limit, &
            allow_zero_mode, retained_count, omitted_count, zero_mode_count, &
            total_energy, omitted_energy, status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        integer, intent(in) :: degree_limit, order_limit
        logical, intent(in) :: allow_zero_mode
        integer, intent(out) :: retained_count, omitted_count, zero_mode_count
        real(dp), intent(out) :: total_energy, omitted_energy
        integer, intent(out) :: status

        integer :: mode

        retained_count = 0
        omitted_count = 0
        zero_mode_count = 0
        total_energy = 0.0_dp
        omitted_energy = 0.0_dp
        status = 1
        if (.not. valid_inputs( &
            degree_indices, orders, coefficients, degree_limit, order_limit)) return
        do mode = 1, size(coefficients)
            total_energy = total_energy + abs(coefficients(mode))**2
            if (degree_indices(mode) == 0 .and. orders(mode) == 0) then
                zero_mode_count = zero_mode_count + 1
            end if
            if (degree_indices(mode) <= degree_limit .and. &
                orders(mode) <= order_limit) then
                retained_count = retained_count + 1
            else
                omitted_count = omitted_count + 1
                omitted_energy = omitted_energy + abs(coefficients(mode))**2
            end if
        end do
        if (.not. allow_zero_mode .and. zero_mode_count > 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine analyze_toroidal_spectral_modes

    subroutine analyze_toroidal_spectral_modes_jvp( &
            degree_indices, orders, coefficients, degree_limit, order_limit, &
            allow_zero_mode, coefficients_dot, total_energy_dot, omitted_energy_dot, &
            status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:), coefficients_dot(:)
        integer, intent(in) :: degree_limit, order_limit
        logical, intent(in) :: allow_zero_mode
        real(dp), intent(out) :: total_energy_dot, omitted_energy_dot
        integer, intent(out) :: status

        integer :: mode
        integer :: retained_count, omitted_count, zero_mode_count
        real(dp) :: total_energy, omitted_energy

        total_energy_dot = 0.0_dp
        omitted_energy_dot = 0.0_dp
        call analyze_toroidal_spectral_modes( &
            degree_indices, orders, coefficients, degree_limit, order_limit, &
            allow_zero_mode, retained_count, omitted_count, zero_mode_count, &
            total_energy, omitted_energy, status)
        if (status /= 0 .or. size(coefficients_dot) /= size(coefficients) .or. &
            .not. finite_complex(coefficients_dot)) then
            status = 1
            return
        end if
        do mode = 1, size(coefficients)
            total_energy_dot = total_energy_dot + 2.0_dp*real( &
                conjg(coefficients(mode))*coefficients_dot(mode), dp)
            if (degree_indices(mode) > degree_limit .or. &
                orders(mode) > order_limit) then
                omitted_energy_dot = omitted_energy_dot + 2.0_dp*real( &
                    conjg(coefficients(mode))*coefficients_dot(mode), dp)
            end if
        end do
        status = 0
    end subroutine analyze_toroidal_spectral_modes_jvp

    subroutine analyze_toroidal_spectral_modes_vjp( &
            degree_indices, orders, coefficients, degree_limit, order_limit, &
            allow_zero_mode, total_energy_bar, omitted_energy_bar, coefficients_bar, &
            status)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        integer, intent(in) :: degree_limit, order_limit
        logical, intent(in) :: allow_zero_mode
        real(dp), intent(in) :: total_energy_bar, omitted_energy_bar
        complex(dp), intent(out) :: coefficients_bar(:)
        integer, intent(out) :: status

        integer :: mode, retained_count, omitted_count, zero_mode_count
        real(dp) :: total_energy, omitted_energy, mode_bar

        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        call analyze_toroidal_spectral_modes( &
            degree_indices, orders, coefficients, degree_limit, order_limit, &
            allow_zero_mode, retained_count, omitted_count, zero_mode_count, &
            total_energy, omitted_energy, status)
        if (status /= 0 .or. size(coefficients_bar) /= size(coefficients) .or. &
            .not. ieee_is_finite(total_energy_bar) .or. &
            .not. ieee_is_finite(omitted_energy_bar)) then
            status = 1
            return
        end if
        do mode = 1, size(coefficients)
            mode_bar = total_energy_bar
            if (degree_indices(mode) > degree_limit .or. &
                orders(mode) > order_limit) mode_bar = mode_bar + omitted_energy_bar
            coefficients_bar(mode) = 2.0_dp*mode_bar*coefficients(mode)
        end do
        status = 0
    end subroutine analyze_toroidal_spectral_modes_vjp

    logical function valid_inputs( &
            degree_indices, orders, coefficients, degree_limit, order_limit) result(valid)
        integer, intent(in) :: degree_indices(:), orders(:)
        complex(dp), intent(in) :: coefficients(:)
        integer, intent(in) :: degree_limit, order_limit

        valid = size(degree_indices) > 0 .and. &
            size(orders) == size(degree_indices) .and. &
            size(coefficients) == size(degree_indices) .and. &
            degree_limit >= 0 .and. order_limit >= 0 .and. &
            all(degree_indices >= 0) .and. all(orders >= 0) .and. &
            finite_complex(coefficients)
    end function valid_inputs

    pure logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_toroidal_spectral_metadata
