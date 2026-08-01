module fortfem_fourier_mode_registry
    !! Fixed-topology metadata and basis evaluation for Fourier-FEM modes.
    !!
    !! A toroidal index is measured per field period.  The phase convention is
    !! therefore
    !!
    !!   m (theta + theta_phase)
    !!     + n * field_periods * (phi + phi_phase).
    !!
    !! The registry is deliberately neutral: it does not choose a plasma
    !! coordinate system, a profile, or a nonlinear closure.  It records the
    !! normalization, conjugate packing, triad lookup, and caller-selected
    !! radial regularity power needed by Fourier-FEM and IGA clients.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: fourier_mode_registry_t
        integer, allocatable :: poloidal_modes(:)
        integer, allocatable :: toroidal_modes(:)
        integer, allocatable :: radial_powers(:)
        real(dp), allocatable :: normalization(:)
        integer :: field_periods = 1
        real(dp) :: poloidal_phase = 0.0_dp
        real(dp) :: toroidal_phase = 0.0_dp
        logical :: real_packed = .false.
    contains
        procedure, private :: assign_fourier_mode_registry
        generic :: assignment(=) => assign_fourier_mode_registry
    end type fourier_mode_registry_t

    public :: evaluate_fourier_mode
    public :: evaluate_fourier_mode_jvp
    public :: evaluate_fourier_mode_vjp
    public :: find_fourier_mode
    public :: fourier_mode_conjugate_index
    public :: fourier_mode_triad_closed
    public :: build_fourier_mode_triad_map
    public :: build_fourier_mode_padded_registry
    public :: build_fourier_mode_closure_registry
    public :: initialize_fourier_mode_registry
    public :: validate_fourier_mode_registry

contains

    subroutine assign_fourier_mode_registry(lhs, rhs)
        class(fourier_mode_registry_t), intent(out) :: lhs
        type(fourier_mode_registry_t), intent(in) :: rhs

        lhs%field_periods = rhs%field_periods
        lhs%poloidal_phase = rhs%poloidal_phase
        lhs%toroidal_phase = rhs%toroidal_phase
        lhs%real_packed = rhs%real_packed
        if (allocated(rhs%poloidal_modes)) then
            allocate(lhs%poloidal_modes(size(rhs%poloidal_modes)))
            lhs%poloidal_modes = rhs%poloidal_modes
        end if
        if (allocated(rhs%toroidal_modes)) then
            allocate(lhs%toroidal_modes(size(rhs%toroidal_modes)))
            lhs%toroidal_modes = rhs%toroidal_modes
        end if
        if (allocated(rhs%radial_powers)) then
            allocate(lhs%radial_powers(size(rhs%radial_powers)))
            lhs%radial_powers = rhs%radial_powers
        end if
        if (allocated(rhs%normalization)) then
            allocate(lhs%normalization(size(rhs%normalization)))
            lhs%normalization = rhs%normalization
        end if
    end subroutine assign_fourier_mode_registry

    subroutine initialize_fourier_mode_registry( &
            registry, poloidal_modes, toroidal_modes, field_periods, &
            poloidal_phase, toroidal_phase, real_packed, radial_powers, &
            normalization, status)
        type(fourier_mode_registry_t), intent(out) :: registry
        integer, intent(in) :: poloidal_modes(:), toroidal_modes(:)
        integer, intent(in) :: field_periods
        real(dp), intent(in) :: poloidal_phase, toroidal_phase
        logical, intent(in) :: real_packed
        integer, intent(in), optional :: radial_powers(:)
        real(dp), intent(in), optional :: normalization(:)
        type(fortsparse_status_t), intent(out) :: status

        call clear_registry(registry)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier mode registry received incompatible optional arrays")
        if (size(poloidal_modes) /= size(toroidal_modes)) return
        if (present(radial_powers) .and. &
            size(radial_powers) /= size(poloidal_modes)) return
        if (present(normalization) .and. &
            size(normalization) /= size(poloidal_modes)) return
        registry%field_periods = field_periods
        registry%poloidal_phase = poloidal_phase
        registry%toroidal_phase = toroidal_phase
        registry%real_packed = real_packed
        if (size(poloidal_modes) > 0) then
            allocate(registry%poloidal_modes(size(poloidal_modes)), &
                registry%toroidal_modes(size(toroidal_modes)), &
                registry%radial_powers(size(poloidal_modes)), &
                registry%normalization(size(poloidal_modes)))
            registry%poloidal_modes = poloidal_modes
            registry%toroidal_modes = toroidal_modes
            registry%radial_powers = abs(poloidal_modes)
            registry%normalization = 1.0_dp
            if (present(radial_powers)) then
                if (size(radial_powers) == size(poloidal_modes)) then
                    registry%radial_powers = radial_powers
                end if
            end if
            if (present(normalization)) then
                if (size(normalization) == size(poloidal_modes)) then
                    registry%normalization = normalization
                end if
            end if
        end if
        if (.not. validate_fourier_mode_registry(registry, status)) return
    end subroutine initialize_fourier_mode_registry

    logical function validate_fourier_mode_registry(registry, status) &
            result(valid)
        type(fourier_mode_registry_t), intent(in) :: registry
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, other, mode_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier mode registry has incompatible metadata")
        if (.not. allocated(registry%poloidal_modes) .or. &
            .not. allocated(registry%toroidal_modes) .or. &
            .not. allocated(registry%radial_powers) .or. &
            .not. allocated(registry%normalization)) return
        mode_count = size(registry%poloidal_modes)
        if (mode_count < 1 .or. size(registry%toroidal_modes) /= mode_count .or. &
            size(registry%radial_powers) /= mode_count .or. &
            size(registry%normalization) /= mode_count .or. &
            registry%field_periods < 1 .or. &
            .not. ieee_is_finite(registry%poloidal_phase) .or. &
            .not. ieee_is_finite(registry%toroidal_phase)) return
        if (any(registry%radial_powers < 0) .or. &
            any(.not. ieee_is_finite(registry%normalization)) .or. &
            any(registry%normalization <= 0.0_dp)) return
        do mode = 1, mode_count
            do other = mode + 1, mode_count
                if (registry%poloidal_modes(mode) == &
                    registry%poloidal_modes(other) .and. &
                    registry%toroidal_modes(mode) == &
                    registry%toroidal_modes(other)) return
            end do
        end do
        if (registry%real_packed) then
            do mode = 1, mode_count
                other = find_fourier_mode(registry, &
                    -registry%poloidal_modes(mode), &
                    -registry%toroidal_modes(mode))
                if (other == 0) return
                if (abs(registry%normalization(other) - &
                    registry%normalization(mode)) > 1.0e-12_dp) return
            end do
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_fourier_mode_registry

    integer function find_fourier_mode(registry, poloidal_mode, toroidal_mode)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: poloidal_mode, toroidal_mode
        integer :: mode

        find_fourier_mode = 0
        if (.not. allocated(registry%poloidal_modes) .or. &
            .not. allocated(registry%toroidal_modes)) return
        do mode = 1, size(registry%poloidal_modes)
            if (registry%poloidal_modes(mode) /= poloidal_mode .or. &
                registry%toroidal_modes(mode) /= toroidal_mode) cycle
            find_fourier_mode = mode
            return
        end do
    end function find_fourier_mode

    integer function fourier_mode_conjugate_index( &
            registry, mode_index, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: mode_index
        type(fortsparse_status_t), intent(out) :: status

        fourier_mode_conjugate_index = 0
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier conjugate lookup received an invalid mode")
        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (mode_index < 1 .or. mode_index > size(registry%poloidal_modes)) return
        fourier_mode_conjugate_index = find_fourier_mode(registry, &
            -registry%poloidal_modes(mode_index), &
            -registry%toroidal_modes(mode_index))
        if (fourier_mode_conjugate_index == 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier conjugate mode is not retained")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end function fourier_mode_conjugate_index

    logical function fourier_mode_triad_closed( &
            registry, first_index, second_index, output_index, status) &
            result(closed)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: first_index, second_index, output_index
        type(fortsparse_status_t), intent(out) :: status
        integer :: expected_poloidal, expected_toroidal

        closed = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier triad lookup received an invalid mode")
        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (first_index < 1 .or. first_index > size(registry%poloidal_modes) .or. &
            second_index < 1 .or. second_index > size(registry%poloidal_modes) .or. &
            output_index < 1 .or. output_index > size(registry%poloidal_modes)) return
        expected_poloidal = registry%poloidal_modes(first_index) + &
            registry%poloidal_modes(second_index)
        expected_toroidal = registry%toroidal_modes(first_index) + &
            registry%toroidal_modes(second_index)
        closed = registry%poloidal_modes(output_index) == expected_poloidal .and. &
            registry%toroidal_modes(output_index) == expected_toroidal
        call status_set(status, FORTSPARSE_OK, "")
    end function fourier_mode_triad_closed

    subroutine build_fourier_mode_triad_map( &
            registry, triad_map, missing_count, status)
        !! Map every retained input pair to its retained output mode.
        !!
        !! `triad_map(first,second)` is the output index for the sum of the
        !! two input mode labels, or zero when that sum is not retained.
        !! `missing_count` is therefore an explicit de-aliasing diagnostic;
        !! callers may reject the registry, pad it, or intentionally project
        !! the omitted sums according to their nonlinear policy.
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, allocatable, intent(out) :: triad_map(:, :)
        integer, intent(out) :: missing_count
        type(fortsparse_status_t), intent(out) :: status
        integer :: first_mode, second_mode, mode_count

        if (allocated(triad_map)) deallocate(triad_map)
        missing_count = 0
        if (.not. validate_fourier_mode_registry(registry, status)) then
            allocate(triad_map(0, 0))
            return
        end if
        mode_count = size(registry%poloidal_modes)
        allocate(triad_map(mode_count, mode_count))
        triad_map = 0
        do first_mode = 1, mode_count
            do second_mode = 1, mode_count
                triad_map(first_mode, second_mode) = find_fourier_mode( &
                    registry, registry%poloidal_modes(first_mode) + &
                    registry%poloidal_modes(second_mode), &
                    registry%toroidal_modes(first_mode) + &
                    registry%toroidal_modes(second_mode))
                if (triad_map(first_mode, second_mode) == 0) then
                    missing_count = missing_count + 1
                end if
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fourier_mode_triad_map

    subroutine build_fourier_mode_padded_registry( &
            registry, padded_registry, status)
        !! Build a deterministic work registry closed under one bilinear sum.
        !!
        !! Original modes are retained first, followed by unique sums of all
        !! ordered input pairs. Original radial powers and normalizations are
        !! preserved; new modes use the regularity default `abs(m)` and unit
        !! normalization. This is a one-product padding primitive, not a
        !! promise of closure under arbitrarily many nonlinear products.
        type(fourier_mode_registry_t), intent(in) :: registry
        type(fourier_mode_registry_t), intent(out) :: padded_registry
        type(fortsparse_status_t), intent(out) :: status
        integer, allocatable :: poloidal_buffer(:), toroidal_buffer(:)
        integer, allocatable :: radial_buffer(:)
        real(dp), allocatable :: normalization_buffer(:)
        integer :: mode_count, maximum_count, mode_total
        integer :: first_mode, second_mode

        call clear_registry(padded_registry)
        if (.not. validate_fourier_mode_registry(registry, status)) return
        mode_count = size(registry%poloidal_modes)
        maximum_count = mode_count + mode_count*mode_count
        allocate(poloidal_buffer(maximum_count), toroidal_buffer(maximum_count), &
            radial_buffer(maximum_count), normalization_buffer(maximum_count))
        mode_total = 0
        do first_mode = 1, mode_count
            call append_padded_mode( &
                registry, registry%poloidal_modes(first_mode), &
                registry%toroidal_modes(first_mode), mode_total, poloidal_buffer, &
                toroidal_buffer, radial_buffer, normalization_buffer)
        end do
        do first_mode = 1, mode_count
            do second_mode = 1, mode_count
                call append_padded_mode( &
                    registry, registry%poloidal_modes(first_mode) + &
                    registry%poloidal_modes(second_mode), &
                    registry%toroidal_modes(first_mode) + &
                    registry%toroidal_modes(second_mode), mode_total, &
                    poloidal_buffer, toroidal_buffer, radial_buffer, &
                    normalization_buffer)
            end do
        end do
        call initialize_fourier_mode_registry( &
            padded_registry, poloidal_buffer(:mode_total), &
            toroidal_buffer(:mode_total), registry%field_periods, &
            registry%poloidal_phase, registry%toroidal_phase, registry%real_packed, &
            radial_powers=radial_buffer(:mode_total), &
            normalization=normalization_buffer(:mode_total), status=status)
    end subroutine build_fourier_mode_padded_registry

    subroutine build_fourier_mode_closure_registry( &
            registry, closure_rounds, closure_registry, status)
        !! Apply the one-product padding primitive for a bounded number of rounds.
        !!
        !! One round is exactly `build_fourier_mode_padded_registry`. Each later
        !! round closes the preceding work set under another bilinear sum. This
        !! makes repeated-product growth explicit and caller-controlled; it is
        !! not an unbounded algebraic closure.
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: closure_rounds
        type(fourier_mode_registry_t), intent(out) :: closure_registry
        type(fortsparse_status_t), intent(out) :: status
        type(fourier_mode_registry_t) :: current_registry, next_registry
        integer :: round

        call clear_registry(closure_registry)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier closure requires a positive round count")
        if (closure_rounds < 1) return
        if (.not. validate_fourier_mode_registry(registry, status)) return
        current_registry = registry
        do round = 1, closure_rounds
            call build_fourier_mode_padded_registry( &
                current_registry, next_registry, status)
            if (status%code /= FORTSPARSE_OK) then
                call clear_registry(closure_registry)
                return
            end if
            current_registry = next_registry
        end do
        closure_registry = current_registry
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine build_fourier_mode_closure_registry

    subroutine append_padded_mode( &
            registry, poloidal_mode, toroidal_mode, mode_total, poloidal_buffer, &
            toroidal_buffer, radial_buffer, normalization_buffer)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: poloidal_mode, toroidal_mode
        integer, intent(inout) :: mode_total
        integer, intent(inout) :: poloidal_buffer(:), toroidal_buffer(:)
        integer, intent(inout) :: radial_buffer(:)
        real(dp), intent(inout) :: normalization_buffer(:)
        integer :: existing_mode, source_mode

        do existing_mode = 1, mode_total
            if (poloidal_buffer(existing_mode) /= poloidal_mode .or. &
                toroidal_buffer(existing_mode) /= toroidal_mode) cycle
            return
        end do
        mode_total = mode_total + 1
        poloidal_buffer(mode_total) = poloidal_mode
        toroidal_buffer(mode_total) = toroidal_mode
        source_mode = find_fourier_mode(registry, poloidal_mode, toroidal_mode)
        if (source_mode > 0) then
            radial_buffer(mode_total) = registry%radial_powers(source_mode)
            normalization_buffer(mode_total) = registry%normalization(source_mode)
        else
            radial_buffer(mode_total) = abs(poloidal_mode)
            normalization_buffer(mode_total) = 1.0_dp
        end if
    end subroutine append_padded_mode

    subroutine evaluate_fourier_mode( &
            registry, mode_index, radius, theta, phi, value, radial_derivative, &
            theta_derivative, phi_derivative, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: mode_index
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(out) :: value, radial_derivative, theta_derivative
        complex(dp), intent(out) :: phi_derivative
        type(fortsparse_status_t), intent(out) :: status

        real(dp) :: radial_factor, radial_factor_derivative, phase
        complex(dp) :: phase_factor

        value = cmplx(0.0_dp, 0.0_dp, dp)
        radial_derivative = value
        theta_derivative = value
        phi_derivative = value
        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (mode_index < 1 .or. mode_index > size(registry%poloidal_modes) .or. &
            .not. ieee_is_finite(radius) .or. radius < 0.0_dp .or. &
            .not. ieee_is_finite(theta) .or. .not. ieee_is_finite(phi)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode evaluation received invalid coordinates")
            return
        end if
        call radial_factor_and_derivative( &
            radius, registry%radial_powers(mode_index), radial_factor, &
            radial_factor_derivative)
        phase = real(registry%poloidal_modes(mode_index), dp)* &
            (theta + registry%poloidal_phase) + &
            real(registry%toroidal_modes(mode_index), dp)* &
            real(registry%field_periods, dp)*(phi + registry%toroidal_phase)
        phase_factor = exp(cmplx(0.0_dp, phase, dp))
        value = registry%normalization(mode_index)*radial_factor*phase_factor
        radial_derivative = registry%normalization(mode_index)* &
            radial_factor_derivative*phase_factor
        theta_derivative = cmplx(0.0_dp, &
            real(registry%poloidal_modes(mode_index), dp), dp)*value
        phi_derivative = cmplx(0.0_dp, real(registry%toroidal_modes(mode_index)* &
            registry%field_periods, dp), dp)*value
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode

    subroutine evaluate_fourier_mode_jvp( &
            registry, mode_index, radius, theta, phi, radius_dot, theta_dot, &
            phi_dot, value_dot, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: mode_index
        real(dp), intent(in) :: radius, theta, phi, radius_dot, theta_dot, phi_dot
        complex(dp), intent(out) :: value_dot
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: value, radial_derivative, theta_derivative, phi_derivative

        value_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call evaluate_fourier_mode( &
            registry, mode_index, radius, theta, phi, value, radial_derivative, &
            theta_derivative, phi_derivative, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(radius_dot) .or. &
            .not. ieee_is_finite(theta_dot) .or. .not. ieee_is_finite(phi_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode JVP received non-finite increments")
            return
        end if
        value_dot = radial_derivative*radius_dot + theta_derivative*theta_dot + &
            phi_derivative*phi_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_jvp

    subroutine evaluate_fourier_mode_vjp( &
            registry, mode_index, radius, theta, phi, value_bar, radius_bar, &
            theta_bar, phi_bar, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: mode_index
        real(dp), intent(in) :: radius, theta, phi
        complex(dp), intent(in) :: value_bar
        real(dp), intent(out) :: radius_bar, theta_bar, phi_bar
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: value, radial_derivative, theta_derivative, phi_derivative

        radius_bar = 0.0_dp
        theta_bar = 0.0_dp
        phi_bar = 0.0_dp
        call evaluate_fourier_mode( &
            registry, mode_index, radius, theta, phi, value, radial_derivative, &
            theta_derivative, phi_derivative, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. ieee_is_finite(real(value_bar, dp)) .or. &
            .not. ieee_is_finite(aimag(value_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode VJP received a non-finite cotangent")
            return
        end if
        radius_bar = real(conjg(value_bar)*radial_derivative, dp)
        theta_bar = real(conjg(value_bar)*theta_derivative, dp)
        phi_bar = real(conjg(value_bar)*phi_derivative, dp)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_fourier_mode_vjp

    subroutine radial_factor_and_derivative( &
            radius, radial_power, factor, derivative)
        real(dp), intent(in) :: radius
        integer, intent(in) :: radial_power
        real(dp), intent(out) :: factor, derivative

        if (radius == 0.0_dp) then
            select case (radial_power)
            case (0)
                factor = 1.0_dp
                derivative = 0.0_dp
            case (1)
                factor = 0.0_dp
                derivative = 1.0_dp
            case default
                factor = 0.0_dp
                derivative = 0.0_dp
            end select
        else
            factor = radius**radial_power
            if (radial_power == 0) then
                derivative = 0.0_dp
            else
                derivative = real(radial_power, dp)*radius**(radial_power - 1)
            end if
        end if
    end subroutine radial_factor_and_derivative

    subroutine clear_registry(registry)
        type(fourier_mode_registry_t), intent(out) :: registry

        if (allocated(registry%poloidal_modes)) deallocate(registry%poloidal_modes)
        if (allocated(registry%toroidal_modes)) deallocate(registry%toroidal_modes)
        if (allocated(registry%radial_powers)) deallocate(registry%radial_powers)
        if (allocated(registry%normalization)) deallocate(registry%normalization)
        registry%field_periods = 1
        registry%poloidal_phase = 0.0_dp
        registry%toroidal_phase = 0.0_dp
        registry%real_packed = .false.
    end subroutine clear_registry

end module fortfem_fourier_mode_registry
