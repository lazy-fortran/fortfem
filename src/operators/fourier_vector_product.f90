module fortfem_fourier_vector_product
    !! Pointwise bilinear Fourier-mode coupling for vector or tensor fields.
    !!
    !! For retained modes n, the product at mode k contains every pair p+q=k:
    !!
    !!   C[a,b,c] L[b,x,p] R[c,x,q].
    !!
    !! The coupling tensor is caller-owned, so the same primitive covers
    !! scalar, vector, H(curl), H(div), and nonlinear tensor-valued clients.
    !! Truncated pairs are omitted explicitly; no hidden de-aliasing rule is
    !! imposed by this neutral algebraic layer.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fourier_mode_registry, only: &
        find_fourier_mode, fourier_mode_registry_t, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_fourier_vector_product
    public :: assemble_fourier_vector_product_jvp
    public :: assemble_fourier_vector_product_vjp

contains

    subroutine assemble_fourier_vector_product( &
            registry, coupling, left_values, right_values, product_values, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coupling(:, :, :)
        complex(dp), intent(in) :: left_values(:, :, :), right_values(:, :, :)
        complex(dp), intent(out) :: product_values(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: first_mode, second_mode, output_mode, point
        integer :: output_component, left_component, right_component
        complex(dp) :: term

        product_values = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_product_inputs( &
            registry, coupling, left_values, right_values, product_values, &
            status)) return
        do first_mode = 1, size(registry%poloidal_modes)
            do second_mode = 1, size(registry%poloidal_modes)
                output_mode = find_fourier_mode(registry, &
                    registry%poloidal_modes(first_mode) + &
                    registry%poloidal_modes(second_mode), &
                    registry%toroidal_modes(first_mode) + &
                    registry%toroidal_modes(second_mode))
                if (output_mode == 0) cycle
                do point = 1, size(left_values, 2)
                    do output_component = 1, size(product_values, 1)
                        term = cmplx(0.0_dp, 0.0_dp, dp)
                        do left_component = 1, size(left_values, 1)
                            do right_component = 1, size(right_values, 1)
                                term = term + coupling(output_component, &
                                    left_component, right_component)* &
                                    left_values(left_component, point, first_mode)* &
                                    right_values(right_component, point, second_mode)
                            end do
                        end do
                        product_values(output_component, point, output_mode) = &
                            product_values(output_component, point, output_mode) + term
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fourier_vector_product

    subroutine assemble_fourier_vector_product_jvp( &
            registry, coupling, left_values, right_values, coupling_dot, &
            left_values_dot, right_values_dot, product_values_dot, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coupling(:, :, :)
        complex(dp), intent(in) :: left_values(:, :, :), right_values(:, :, :)
        complex(dp), intent(in) :: coupling_dot(:, :, :)
        complex(dp), intent(in) :: left_values_dot(:, :, :), right_values_dot(:, :, :)
        complex(dp), intent(out) :: product_values_dot(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: first_mode, second_mode, output_mode, point
        integer :: output_component, left_component, right_component
        complex(dp) :: term_dot

        product_values_dot = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_product_inputs( &
            registry, coupling, left_values, right_values, product_values_dot, &
            status)) return
        if (.not. validate_product_direction( &
            coupling_dot, left_values_dot, right_values_dot, product_values_dot, &
            size(coupling, 1), size(left_values, 1), size(right_values, 1), &
            size(left_values, 2), size(registry%poloidal_modes))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier vector-product JVP has incompatible increments")
            return
        end if
        do first_mode = 1, size(registry%poloidal_modes)
            do second_mode = 1, size(registry%poloidal_modes)
                output_mode = find_fourier_mode(registry, &
                    registry%poloidal_modes(first_mode) + &
                    registry%poloidal_modes(second_mode), &
                    registry%toroidal_modes(first_mode) + &
                    registry%toroidal_modes(second_mode))
                if (output_mode == 0) cycle
                do point = 1, size(left_values, 2)
                    do output_component = 1, size(product_values_dot, 1)
                        term_dot = cmplx(0.0_dp, 0.0_dp, dp)
                        do left_component = 1, size(left_values, 1)
                            do right_component = 1, size(right_values, 1)
                                term_dot = term_dot + &
                                    coupling_dot(output_component, left_component, &
                                    right_component)*left_values( &
                                    left_component, point, first_mode)* &
                                    right_values( &
                                    right_component, point, second_mode) + &
                                    coupling(output_component, left_component, &
                                    right_component)*left_values_dot( &
                                    left_component, point, first_mode)* &
                                    right_values( &
                                    right_component, point, second_mode) + &
                                    coupling(output_component, left_component, &
                                    right_component)*left_values( &
                                    left_component, point, first_mode)* &
                                    right_values_dot( &
                                    right_component, point, second_mode)
                            end do
                        end do
                        product_values_dot(output_component, point, output_mode) = &
                            product_values_dot(output_component, point, output_mode) + &
                            term_dot
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fourier_vector_product_jvp

    subroutine assemble_fourier_vector_product_vjp( &
            registry, coupling, left_values, right_values, product_values_bar, &
            left_values_bar, right_values_bar, coupling_bar, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coupling(:, :, :)
        complex(dp), intent(in) :: left_values(:, :, :), right_values(:, :, :)
        complex(dp), intent(in) :: product_values_bar(:, :, :)
        complex(dp), intent(out) :: left_values_bar(:, :, :), right_values_bar(:, :, :)
        complex(dp), intent(out) :: coupling_bar(:, :, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: first_mode, second_mode, output_mode, point
        integer :: output_component, left_component, right_component
        complex(dp) :: cotangent

        left_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        right_values_bar = cmplx(0.0_dp, 0.0_dp, dp)
        coupling_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_product_inputs( &
            registry, coupling, left_values, right_values, product_values_bar, &
            status)) return
        if (.not. validate_product_adjoint( &
            left_values_bar, right_values_bar, coupling_bar, product_values_bar, &
            size(coupling, 1), size(left_values, 1), size(right_values, 1), &
            size(left_values, 2), size(registry%poloidal_modes))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier vector-product VJP has incompatible cotangents")
            return
        end if
        do first_mode = 1, size(registry%poloidal_modes)
            do second_mode = 1, size(registry%poloidal_modes)
                output_mode = find_fourier_mode(registry, &
                    registry%poloidal_modes(first_mode) + &
                    registry%poloidal_modes(second_mode), &
                    registry%toroidal_modes(first_mode) + &
                    registry%toroidal_modes(second_mode))
                if (output_mode == 0) cycle
                do point = 1, size(left_values, 2)
                    do output_component = 1, size(product_values_bar, 1)
                        cotangent = product_values_bar( &
                            output_component, point, output_mode)
                        do left_component = 1, size(left_values, 1)
                            do right_component = 1, size(right_values, 1)
                                left_values_bar(left_component, point, first_mode) = &
                                    left_values_bar( &
                                    left_component, point, first_mode) + &
                                    conjg(coupling(output_component, left_component, &
                                    right_component)*right_values( &
                                    right_component, point, second_mode))*cotangent
                                right_values_bar( &
                                    right_component, point, second_mode) = &
                                    right_values_bar( &
                                    right_component, point, second_mode) + &
                                    conjg(coupling(output_component, left_component, &
                                    right_component)*left_values( &
                                    left_component, point, first_mode))* &
                                    cotangent
                                coupling_bar(output_component, left_component, &
                                    right_component) = &
                                    coupling_bar(output_component, left_component, &
                                    right_component) + &
                                    conjg(left_values( &
                                    left_component, point, first_mode)* &
                                    right_values( &
                                    right_component, point, second_mode))* &
                                    cotangent
                            end do
                        end do
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_fourier_vector_product_vjp

    logical function validate_product_inputs( &
            registry, coupling, left_values, right_values, product_values, status) &
            result(valid)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coupling(:, :, :)
        complex(dp), intent(in) :: left_values(:, :, :), right_values(:, :, :)
        complex(dp), intent(in) :: product_values(:, :, :)
        type(fortsparse_status_t), intent(out) :: status
        integer :: output_count, left_count, right_count, point_count, mode_count

        valid = .false.
        if (.not. validate_fourier_mode_registry(registry, status)) return
        output_count = size(coupling, 1)
        left_count = size(coupling, 2)
        right_count = size(coupling, 3)
        point_count = size(left_values, 2)
        mode_count = size(registry%poloidal_modes)
        if (output_count < 1 .or. left_count < 1 .or. right_count < 1 .or. &
            point_count < 1 .or. size(left_values, 1) /= left_count .or. &
            size(right_values, 1) /= right_count .or. &
            size(left_values, 3) /= mode_count .or. &
            size(right_values, 2) /= point_count .or. &
            size(right_values, 3) /= mode_count .or. &
            size(product_values, 1) /= output_count .or. &
            size(product_values, 2) /= point_count .or. &
            size(product_values, 3) /= mode_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier vector-product arrays have incompatible shapes")
            return
        end if
        if (.not. finite_complex_3d(coupling) .or. &
            .not. finite_complex_3d(left_values) .or. &
            .not. finite_complex_3d(right_values) .or. &
            .not. finite_complex_3d(product_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier vector-product arrays contain non-finite data")
            return
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_product_inputs

    logical function validate_product_direction( &
            coupling_dot, left_values_dot, right_values_dot, product_values_dot, &
            output_count, left_count, right_count, point_count, mode_count) &
            result(valid)
        complex(dp), intent(in) :: coupling_dot(:, :, :)
        complex(dp), intent(in) :: left_values_dot(:, :, :), right_values_dot(:, :, :)
        complex(dp), intent(in) :: product_values_dot(:, :, :)
        integer, intent(in) :: output_count, left_count, right_count
        integer, intent(in) :: point_count, mode_count

        valid = size(coupling_dot, 1) == output_count .and. &
            size(coupling_dot, 2) == left_count .and. &
            size(coupling_dot, 3) == right_count .and. &
            size(left_values_dot, 1) == left_count .and. &
            size(left_values_dot, 2) == point_count .and. &
            size(left_values_dot, 3) == mode_count .and. &
            size(right_values_dot, 1) == right_count .and. &
            size(right_values_dot, 2) == point_count .and. &
            size(right_values_dot, 3) == mode_count .and. &
            size(product_values_dot, 1) == output_count .and. &
            size(product_values_dot, 2) == point_count .and. &
            size(product_values_dot, 3) == mode_count .and. &
            finite_complex_3d(coupling_dot) .and. &
            finite_complex_3d(left_values_dot) .and. &
            finite_complex_3d(right_values_dot)
    end function validate_product_direction

    logical function validate_product_adjoint( &
            left_values_bar, right_values_bar, coupling_bar, product_values_bar, &
            output_count, left_count, right_count, point_count, mode_count) &
            result(valid)
        complex(dp), intent(in) :: left_values_bar(:, :, :), right_values_bar(:, :, :)
        complex(dp), intent(in) :: coupling_bar(:, :, :), product_values_bar(:, :, :)
        integer, intent(in) :: output_count, left_count, right_count
        integer, intent(in) :: point_count, mode_count

        valid = size(left_values_bar, 1) == left_count .and. &
            size(left_values_bar, 2) == point_count .and. &
            size(left_values_bar, 3) == mode_count .and. &
            size(right_values_bar, 1) == right_count .and. &
            size(right_values_bar, 2) == point_count .and. &
            size(right_values_bar, 3) == mode_count .and. &
            size(coupling_bar, 1) == output_count .and. &
            size(coupling_bar, 2) == left_count .and. &
            size(coupling_bar, 3) == right_count .and. &
            size(product_values_bar, 1) == output_count .and. &
            size(product_values_bar, 2) == point_count .and. &
            size(product_values_bar, 3) == mode_count .and. &
            finite_complex_3d(product_values_bar)
    end function validate_product_adjoint

    pure logical function finite_complex_3d(values) result(valid)
        complex(dp), intent(in) :: values(:, :, :)
        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_3d

end module fortfem_fourier_vector_product
