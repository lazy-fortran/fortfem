module fortfem_fourier_zernike_basis
    !! Neutral Fourier--Zernike mode metadata and scalar radial basis.
    !!
    !! The labels are `(l,m,n)`, where `l` is the radial polynomial degree,
    !! `m` is the poloidal Fourier mode, and `n` is the toroidal Fourier mode.
    !! Only scalar axis-regular pairs with `l >= abs(m)` and
    !! `mod(l-abs(m),2) == 0` are enumerated.  No profiles, equilibrium
    !! equations, coordinate readers, or closure assumptions are included.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    integer, parameter, public :: FOURIER_ZERNIKE_PARITY_EVEN = 0
    integer, parameter, public :: FOURIER_ZERNIKE_PARITY_ODD = 1

    type, public :: fourier_zernike_mode_t
        integer :: l = 0
        integer :: m = 0
        integer :: n = 0
        integer :: minimum_power = 0
        integer :: parity = FOURIER_ZERNIKE_PARITY_EVEN
        integer :: conjugate_index = 0
    end type fourier_zernike_mode_t

    type, public :: fourier_zernike_basis_t
        type(fourier_zernike_mode_t), allocatable :: modes(:)
        integer :: l_max = 0
        integer :: m_max = 0
        integer :: n_max = 0
    end type fourier_zernike_basis_t

    public :: build_fourier_zernike_basis
    public :: validate_fourier_zernike_basis
    public :: fourier_zernike_mode_requirements
    public :: evaluate_fourier_zernike_radial

contains

    subroutine build_fourier_zernike_basis( &
            l_max, m_max, n_max, basis, status)
        !! Enumerate a complete, deterministic, conjugate-safe mode set.
        integer, intent(in) :: l_max, m_max, n_max
        type(fourier_zernike_basis_t), intent(out) :: basis
        type(fortsparse_status_t), intent(out) :: status

        integer :: count, l, m, n, mode, bound

        call clear_basis(basis)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier-Zernike bounds must be non-negative")
        if (l_max < 0 .or. m_max < 0 .or. n_max < 0) return

        count = count_modes(l_max, m_max, n_max)
        if (count < 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier-Zernike mode set is empty")
            return
        end if

        allocate(basis%modes(count))
        basis%l_max = l_max
        basis%m_max = m_max
        basis%n_max = n_max
        mode = 0
        do l = 0, l_max
            bound = min(l, m_max)
            do m = -bound, bound
                if (mod(l - abs(m), 2) /= 0) cycle
                do n = -n_max, n_max
                    mode = mode + 1
                    basis%modes(mode)%l = l
                    basis%modes(mode)%m = m
                    basis%modes(mode)%n = n
                    basis%modes(mode)%minimum_power = abs(m)
                    basis%modes(mode)%parity = modulo(abs(m), 2)
                end do
            end do
        end do

        call sort_modes(basis%modes)
        call assign_conjugate_indices(basis%modes)
        if (.not. validate_fourier_zernike_basis(basis, status)) then
            call clear_basis(basis)
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fourier_zernike_basis

    logical function validate_fourier_zernike_basis(basis, status) &
            result(valid)
        !! Validate labels, metadata, ordering, and conjugate pair links.
        type(fourier_zernike_basis_t), intent(in) :: basis
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, pair, minimum_power, parity
        logical :: admissible

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier-Zernike basis is invalid")
        if (.not. allocated(basis%modes)) return
        if (size(basis%modes) < 1) return
        if (basis%l_max < 0 .or. basis%m_max < 0 .or. basis%n_max < 0) return

        do mode = 1, size(basis%modes)
            if (basis%modes(mode)%l < 0 .or. &
                    basis%modes(mode)%l > basis%l_max) return
            if (abs(basis%modes(mode)%m) > basis%m_max .or. &
                    abs(basis%modes(mode)%m) > basis%modes(mode)%l) return
            if (abs(basis%modes(mode)%n) > basis%n_max) return
            call fourier_zernike_mode_requirements( &
                basis%modes(mode)%l, basis%modes(mode)%m, minimum_power, &
                parity, admissible, status)
            if (.not. admissible .or. status%code /= FORTSPARSE_OK) return
            if (basis%modes(mode)%minimum_power /= minimum_power .or. &
                    basis%modes(mode)%parity /= parity) return
            if (mode > 1) then
                if (.not. mode_precedes(basis%modes(mode - 1), &
                        basis%modes(mode))) return
            end if
            pair = basis%modes(mode)%conjugate_index
            if (pair < 1 .or. pair > size(basis%modes)) return
            if (basis%modes(pair)%l /= basis%modes(mode)%l .or. &
                    basis%modes(pair)%m /= -basis%modes(mode)%m .or. &
                    basis%modes(pair)%n /= -basis%modes(mode)%n) return
            if (basis%modes(pair)%conjugate_index /= mode) return
            if (pair /= mode .and. abs(pair - mode) /= 1) return
        end do
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_fourier_zernike_basis

    subroutine fourier_zernike_mode_requirements( &
            l, m, minimum_power, parity, admissible, status)
        !! Return scalar axis-regularity metadata for one `(l,m)` label.
        integer, intent(in) :: l, m
        integer, intent(out) :: minimum_power, parity
        logical, intent(out) :: admissible
        type(fortsparse_status_t), intent(out) :: status

        minimum_power = abs(m)
        parity = modulo(abs(m), 2)
        admissible = l >= 0 .and. abs(m) <= l .and. &
            modulo(l - abs(m), 2) == 0
        if (admissible) then
            call status_set(status, FORTSPARSE_OK, "")
        else
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier-Zernike label violates scalar axis regularity")
        end if
    end subroutine fourier_zernike_mode_requirements

    subroutine evaluate_fourier_zernike_radial( &
            l, m, rho, values, status)
        !! Evaluate the scalar Zernike radial polynomial `R_l^|m|(rho)`.
        !!
        !! The polynomial is normalized by `R_l^|m|(1)=1`; samples must be in
        !! the closed unit interval.  The implementation uses the finite sum
        !!
        !!   sum_s (-1)^s (l-s)! rho^(l-2s)
        !!       / (s! ((l+|m|)/2-s)! ((l-|m|)/2-s)!).
        integer, intent(in) :: l, m
        real(dp), intent(in) :: rho(:)
        real(dp), intent(out) :: values(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: s, radial_order, minimum_power, parity
        logical :: admissible
        real(dp) :: coefficient

        values = 0.0_dp
        if (size(values) /= size(rho)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier-Zernike radial output shape is incompatible")
            return
        end if
        if (.not. finite_real(rho) .or. any(rho < 0.0_dp) .or. &
                any(rho > 1.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier-Zernike radial samples must lie in [0,1]")
            return
        end if
        call fourier_zernike_mode_requirements( &
            l, m, minimum_power, parity, admissible, status)
        if (.not. admissible .or. status%code /= FORTSPARSE_OK) return
        radial_order = (l - abs(m))/2
        do s = 0, radial_order
            coefficient = zernike_coefficient(l, abs(m), s)
            values = values + coefficient*rho**(l - 2*s)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_zernike_radial

    integer function count_modes(l_max, m_max, n_max) result(count)
        integer, intent(in) :: l_max, m_max, n_max
        integer :: l, m, bound

        count = 0
        do l = 0, l_max
            bound = min(l, m_max)
            do m = -bound, bound
                if (mod(l - abs(m), 2) == 0) then
                    count = count + 2*n_max + 1
                end if
            end do
        end do
    end function count_modes

    real(dp) function zernike_coefficient(l, m, s) result(coefficient)
        integer, intent(in) :: l, m, s

        coefficient = factorial_real(l - s)/ &
            (factorial_real(s)*factorial_real((l + m)/2 - s)* &
            factorial_real((l - m)/2 - s))
        if (mod(s, 2) /= 0) coefficient = -coefficient
    end function zernike_coefficient

    pure real(dp) function factorial_real(number) result(value)
        integer, intent(in) :: number
        integer :: factor

        value = 1.0_dp
        do factor = 2, number
            value = value*real(factor, dp)
        end do
    end function factorial_real

    subroutine assign_conjugate_indices(modes)
        type(fourier_zernike_mode_t), intent(inout) :: modes(:)
        integer :: mode, other

        do mode = 1, size(modes)
            do other = 1, size(modes)
                if (modes(other)%l /= modes(mode)%l) cycle
                if (modes(other)%m /= -modes(mode)%m) cycle
                if (modes(other)%n /= -modes(mode)%n) cycle
                modes(mode)%conjugate_index = other
                exit
            end do
        end do
    end subroutine assign_conjugate_indices

    subroutine sort_modes(modes)
        type(fourier_zernike_mode_t), intent(inout) :: modes(:)
        type(fourier_zernike_mode_t) :: candidate
        integer :: position, previous

        do position = 2, size(modes)
            candidate = modes(position)
            previous = position - 1
            do while (previous >= 1)
                if (.not. mode_precedes(candidate, modes(previous))) exit
                modes(previous + 1) = modes(previous)
                previous = previous - 1
            end do
            modes(previous + 1) = candidate
        end do
    end subroutine sort_modes

    logical function mode_precedes(left, right) result(precedes)
        type(fourier_zernike_mode_t), intent(in) :: left, right
        integer :: left_m, left_n, left_orientation
        integer :: right_m, right_n, right_orientation

        call canonical_pair_key(left%m, left%n, left_m, left_n, left_orientation)
        call canonical_pair_key(right%m, right%n, right_m, right_n, right_orientation)
        precedes = left%l < right%l
        if (left%l == right%l) precedes = left_m < right_m
        if (left%l == right%l .and. left_m == right_m) then
            precedes = left_n < right_n
        end if
        if (left%l == right%l .and. left_m == right_m .and. &
                left_n == right_n) precedes = left_orientation < right_orientation
    end function mode_precedes

    subroutine canonical_pair_key(m, n, key_m, key_n, orientation)
        integer, intent(in) :: m, n
        integer, intent(out) :: key_m, key_n, orientation

        if (m > 0 .or. (m == 0 .and. n >= 0)) then
            key_m = m
            key_n = n
            orientation = 0
        else
            key_m = -m
            key_n = -n
            orientation = 1
        end if
    end subroutine canonical_pair_key

    pure logical function finite_real(values) result(finite)
        real(dp), intent(in) :: values(:)

        finite = all(ieee_is_finite(values))
    end function finite_real

    subroutine clear_basis(basis)
        type(fourier_zernike_basis_t), intent(out) :: basis

        if (allocated(basis%modes)) deallocate(basis%modes)
        basis%l_max = 0
        basis%m_max = 0
        basis%n_max = 0
    end subroutine clear_basis

end module fortfem_fourier_zernike_basis
