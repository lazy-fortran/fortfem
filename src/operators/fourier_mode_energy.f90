module fortfem_fourier_mode_energy
    !! Weighted modal energy and real-packed Fourier symmetry diagnostics.
    !!
    !! For physical modal coefficients c(component, point, mode), positive
    !! point weights w and positive modal weights alpha, the neutral energy is
    !!
    !!   E_mode = 1/2 sum_point,component alpha* w |c|^2.
    !!
    !! The caller owns the metric, quadrature and normalization represented by
    !! those weights.  This block only supplies the fixed-topology algebra and
    !! its real-coordinate derivative actions.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fourier_mode_registry, only: &
        find_fourier_mode, fourier_mode_registry_t, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fourier_mode_energy
    public :: assemble_fourier_mode_energy_jvp
    public :: assemble_fourier_mode_energy_vjp
    public :: fourier_coefficients_conjugate_symmetric

contains

    logical function fourier_coefficients_conjugate_symmetric( &
            registry, coefficients, tolerance, residual, status) result(symmetric)
        !! Check conjugate coefficient packing on a fixed Fourier topology.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :, :)
        real(dp), intent(in) :: tolerance
        real(dp), intent(out) :: residual
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, conjugate
        real(dp) :: difference

        symmetric = .false.
        residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Fourier coefficient symmetry received invalid data")
        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (size(coefficients, 1) < 1 .or. size(coefficients, 2) < 1 .or. &
            size(coefficients, 3) /= size(registry%poloidal_modes) .or. &
            .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) return
        if (.not. finite_complex_3d(coefficients)) return
        if (.not. registry%real_packed) then
            symmetric = .true.
            call status_set(status, FORTSPARSE_OK, "")
            return
        end if
        do mode = 1, size(registry%poloidal_modes)
            conjugate = find_fourier_mode(registry, &
                -registry%poloidal_modes(mode), -registry%toroidal_modes(mode))
            if (conjugate == 0) return
            difference = maxval(abs(coefficients(:, :, conjugate) - &
                conjg(coefficients(:, :, mode))))
            residual = max(residual, difference)
        end do
        symmetric = residual <= tolerance
        call status_set(status, FORTSPARSE_OK, "")
    end function fourier_coefficients_conjugate_symmetric

    subroutine assemble_fourier_mode_energy( &
            registry, coefficients, point_weights, mode_weights, mode_energy, &
            total_energy, status)
        !! Assemble weighted energy per retained Fourier mode and in total.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :, :)
        real(dp), intent(in) :: point_weights(:), mode_weights(:)
        real(dp), intent(out) :: mode_energy(:), total_energy
        type(fortsparse_status_t), intent(out) :: status

        integer :: component, point, mode

        mode_energy = 0.0_dp
        total_energy = 0.0_dp
        if (.not. validate_energy_inputs( &
            registry, coefficients, point_weights, mode_weights, mode_energy, &
            status)) return
        do mode = 1, size(mode_energy)
            do point = 1, size(point_weights)
                do component = 1, size(coefficients, 1)
                    mode_energy(mode) = mode_energy(mode) + 0.5_dp* &
                        mode_weights(mode)*point_weights(point)*real( &
                        conjg(coefficients(component, point, mode))* &
                        coefficients(component, point, mode), dp)
                end do
            end do
        end do
        total_energy = sum(mode_energy)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fourier_mode_energy

    subroutine assemble_fourier_mode_energy_jvp( &
            registry, coefficients, coefficients_dot, point_weights, &
            point_weights_dot, mode_weights, mode_weights_dot, mode_energy_dot, &
            total_energy_dot, status)
        !! Apply the tangent of the weighted modal energy.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :, :), coefficients_dot(:, :, :)
        real(dp), intent(in) :: point_weights(:), point_weights_dot(:)
        real(dp), intent(in) :: mode_weights(:), mode_weights_dot(:)
        real(dp), intent(out) :: mode_energy_dot(:), total_energy_dot
        type(fortsparse_status_t), intent(out) :: status

        integer :: component, point, mode
        real(dp) :: squared_norm, directional_norm

        mode_energy_dot = 0.0_dp
        total_energy_dot = 0.0_dp
        if (.not. validate_energy_inputs( &
            registry, coefficients, point_weights, mode_weights, mode_energy_dot, &
            status)) return
        if (.not. validate_energy_direction( &
            coefficients, coefficients_dot, point_weights_dot, mode_weights_dot, &
            size(point_weights), size(mode_weights))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier modal energy JVP received incompatible increments")
            return
        end if
        do mode = 1, size(mode_energy_dot)
            do point = 1, size(point_weights)
                do component = 1, size(coefficients, 1)
                    squared_norm = real(conjg(coefficients(component, point, mode))* &
                        coefficients(component, point, mode), dp)
                    directional_norm = real(conjg( &
                        coefficients(component, point, mode))*coefficients_dot( &
                        component, point, mode), dp)
                    mode_energy_dot(mode) = mode_energy_dot(mode) + 0.5_dp* &
                        point_weights(point)*mode_weights_dot(mode)*squared_norm + &
                        0.5_dp*point_weights_dot(point)*mode_weights(mode)* &
                        squared_norm + point_weights(point)*mode_weights(mode)* &
                        directional_norm
                end do
            end do
        end do
        total_energy_dot = sum(mode_energy_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fourier_mode_energy_jvp

    subroutine assemble_fourier_mode_energy_vjp( &
            registry, coefficients, point_weights, mode_weights, mode_energy_bar, &
            total_energy_bar, coefficients_bar, point_weights_bar, mode_weights_bar, &
            status)
        !! Apply the real adjoint of the weighted modal energy.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :, :)
        real(dp), intent(in) :: point_weights(:), mode_weights(:)
        real(dp), intent(in) :: mode_energy_bar(:), total_energy_bar
        complex(dp), intent(out) :: coefficients_bar(:, :, :)
        real(dp), intent(out) :: point_weights_bar(:), mode_weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: component, point, mode
        real(dp) :: squared_norm, energy_bar

        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        point_weights_bar = 0.0_dp
        mode_weights_bar = 0.0_dp
        if (.not. validate_energy_inputs( &
            registry, coefficients, point_weights, mode_weights, mode_energy_bar, &
            status)) return
        if (size(coefficients_bar, 1) /= size(coefficients, 1) .or. &
            size(coefficients_bar, 2) /= size(coefficients, 2) .or. &
            size(coefficients_bar, 3) /= size(coefficients, 3) .or. &
            size(point_weights_bar) /= size(point_weights) .or. &
            size(mode_weights_bar) /= size(mode_weights) .or. &
            .not. ieee_is_finite(total_energy_bar) .or. &
            .not. all(ieee_is_finite(mode_energy_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier modal energy VJP received incompatible cotangents")
            return
        end if
        do mode = 1, size(mode_weights)
            energy_bar = mode_energy_bar(mode) + total_energy_bar
            do point = 1, size(point_weights)
                do component = 1, size(coefficients, 1)
                    squared_norm = real(conjg(coefficients(component, point, mode))* &
                        coefficients(component, point, mode), dp)
                    coefficients_bar(component, point, mode) = &
                        coefficients_bar(component, point, mode) + energy_bar* &
                        point_weights(point)*mode_weights(mode)* &
                        coefficients(component, point, mode)
                    point_weights_bar(point) = point_weights_bar(point) + 0.5_dp* &
                        energy_bar*mode_weights(mode)*squared_norm
                    mode_weights_bar(mode) = mode_weights_bar(mode) + 0.5_dp* &
                        energy_bar*point_weights(point)*squared_norm
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fourier_mode_energy_vjp

    logical function validate_energy_inputs( &
            registry, coefficients, point_weights, mode_weights, mode_energy, &
            status) result(valid)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :, :)
        real(dp), intent(in) :: point_weights(:), mode_weights(:)
        real(dp), intent(in) :: mode_energy(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: mode_count, point_count

        valid = .false.
        if (.not. validate_fourier_mode_registry(registry, status)) return
        mode_count = size(registry%poloidal_modes)
        point_count = size(coefficients, 2)
        if (size(coefficients, 1) < 1 .or. point_count < 1 .or. &
            size(coefficients, 3) /= mode_count .or. &
            size(point_weights) /= point_count .or. &
            size(mode_weights) /= mode_count .or. size(mode_energy) /= mode_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier modal energy arrays have incompatible shapes")
            return
        end if
        if (.not. finite_complex_3d(coefficients) .or. &
            .not. all(ieee_is_finite(point_weights)) .or. &
            .not. all(ieee_is_finite(mode_weights)) .or. &
            any(point_weights <= 0.0_dp) .or. any(mode_weights <= 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier modal energy received non-positive or non-finite data")
            return
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_energy_inputs

    logical function validate_energy_direction( &
            coefficients, coefficients_dot, point_weights_dot, mode_weights_dot, &
            point_count, mode_count) result(valid)
        complex(dp), intent(in) :: coefficients(:, :, :), coefficients_dot(:, :, :)
        real(dp), intent(in) :: point_weights_dot(:), mode_weights_dot(:)
        integer, intent(in) :: point_count, mode_count

        valid = all(shape(coefficients_dot) == shape(coefficients)) .and. &
            size(point_weights_dot) == point_count .and. &
            size(mode_weights_dot) == mode_count .and. &
            finite_complex_3d(coefficients_dot) .and. &
            all(ieee_is_finite(point_weights_dot)) .and. &
            all(ieee_is_finite(mode_weights_dot))
    end function validate_energy_direction

    pure logical function finite_complex_3d(values) result(valid)
        complex(dp), intent(in) :: values(:, :, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_3d

end module fortfem_fourier_mode_energy
