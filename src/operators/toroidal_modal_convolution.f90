module fortfem_toroidal_modal_convolution
    !! Matrix-free retained-mode convolution for periodic Green operators.
    !!
    !! For a fixed Fourier registry, the action is
    !!
    !!   y_k = sum_l K_(k-l) x_l,
    !!
    !! where a difference is included only when its mode is retained by the
    !! registry.  The kernel coefficients are caller-owned modal data; no
    !! regularization, Green normalization, or zero-mode convention is hidden
    !! here.  This makes the primitive suitable for NESTOR-like periodic
    !! surface operators while keeping application-specific geometry outside
    !! FortFEM.
    !!
    !! The implementation deliberately scans the mode labels instead of
    !! allocating an O(N^2) difference map.  It therefore has bounded auxiliary
    !! memory for large truncations; callers can provide a specialized map or
    !! accelerated backend when the O(N^3) label lookup becomes limiting.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fourier_mode_registry, only: &
        find_fourier_mode, fourier_mode_registry_t, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: apply_toroidal_modal_convolution
    public :: apply_toroidal_modal_convolution_jvp
    public :: apply_toroidal_modal_convolution_vjp

contains

    subroutine apply_toroidal_modal_convolution( &
            registry, kernel, source, potential, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: kernel(:), source(:)
        complex(dp), intent(out) :: potential(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: kernel_mode, output_mode, source_mode

        potential = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "toroidal modal convolution received incompatible arrays")
        if (.not. valid_inputs(registry, kernel, source, potential)) return
        do output_mode = 1, size(potential)
            do source_mode = 1, size(source)
                kernel_mode = find_fourier_mode(registry, &
                    registry%poloidal_modes(output_mode) - &
                    registry%poloidal_modes(source_mode), &
                    registry%toroidal_modes(output_mode) - &
                    registry%toroidal_modes(source_mode))
                if (kernel_mode == 0) cycle
                potential(output_mode) = potential(output_mode) + &
                    kernel(kernel_mode)*source(source_mode)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_toroidal_modal_convolution

    subroutine apply_toroidal_modal_convolution_jvp( &
            registry, kernel, source, kernel_dot, source_dot, potential_dot, &
            status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: kernel(:), source(:)
        complex(dp), intent(in) :: kernel_dot(:), source_dot(:)
        complex(dp), intent(out) :: potential_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: kernel_mode, output_mode, source_mode

        potential_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "toroidal modal convolution JVP received incompatible arrays")
        if (.not. valid_inputs(registry, kernel, source, potential_dot)) return
        if (size(kernel_dot) /= size(kernel) .or. size(source_dot) /= size(source) .or. &
            .not. finite_complex(kernel_dot) .or. .not. finite_complex(source_dot)) return
        do output_mode = 1, size(source)
            do source_mode = 1, size(source)
                kernel_mode = find_fourier_mode(registry, &
                    registry%poloidal_modes(output_mode) - &
                    registry%poloidal_modes(source_mode), &
                    registry%toroidal_modes(output_mode) - &
                    registry%toroidal_modes(source_mode))
                if (kernel_mode == 0) cycle
                potential_dot(output_mode) = potential_dot(output_mode) + &
                    kernel_dot(kernel_mode)*source(source_mode) + &
                    kernel(kernel_mode)*source_dot(source_mode)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_toroidal_modal_convolution_jvp

    subroutine apply_toroidal_modal_convolution_vjp( &
            registry, kernel, source, potential_bar, kernel_bar, source_bar, &
            status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: kernel(:), source(:), potential_bar(:)
        complex(dp), intent(out) :: kernel_bar(:), source_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: kernel_mode, output_mode, source_mode

        kernel_bar = cmplx(0.0_dp, 0.0_dp, dp)
        source_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "toroidal modal convolution VJP received incompatible arrays")
        if (.not. valid_inputs(registry, kernel, source, potential_bar)) return
        if (size(kernel_bar) /= size(kernel) .or. size(source_bar) /= size(source)) return
        do output_mode = 1, size(source)
            do source_mode = 1, size(source)
                kernel_mode = find_fourier_mode(registry, &
                    registry%poloidal_modes(output_mode) - &
                    registry%poloidal_modes(source_mode), &
                    registry%toroidal_modes(output_mode) - &
                    registry%toroidal_modes(source_mode))
                if (kernel_mode == 0) cycle
                kernel_bar(kernel_mode) = kernel_bar(kernel_mode) + &
                    conjg(source(source_mode))*potential_bar(output_mode)
                source_bar(source_mode) = source_bar(source_mode) + &
                    conjg(kernel(kernel_mode))*potential_bar(output_mode)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_toroidal_modal_convolution_vjp

    logical function valid_inputs( &
            registry, kernel, source, output) result(valid)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: kernel(:), source(:), output(:)
        type(fortsparse_status_t) :: validation_status

        valid = .false.
        if (.not. validate_fourier_mode_registry(registry, validation_status)) return
        if (size(kernel) /= size(registry%poloidal_modes) .or. &
            size(source) /= size(kernel) .or. size(output) /= size(source)) return
        if (.not. finite_complex(kernel) .or. .not. finite_complex(source) .or. &
            .not. finite_complex(output)) return
        valid = .true.
    end function valid_inputs

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_toroidal_modal_convolution
