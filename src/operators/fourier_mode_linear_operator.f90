module fortfem_fourier_mode_linear_operator
    !! Mode-diagonal complex linear maps and their exact derivatives.
    !!
    !! The map is deliberately agnostic to the field interpretation.  For each
    !! retained Fourier mode m, a caller-owned matrix A(:,:,m) maps an input
    !! component vector to an output component vector:
    !!
    !!     y(:,m) = A(:,:,m) c(:,m).
    !!
    !! This is the small linear building block shared by scalar, H(curl),
    !! H(div), and tensor-valued Fourier-FEM compositions.  No physical model,
    !! mode registry, or basis convention is imposed here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: apply_fourier_mode_linear_operator
    public :: apply_fourier_mode_linear_operator_jvp
    public :: apply_fourier_mode_linear_operator_vjp

contains

    subroutine apply_fourier_mode_linear_operator( &
            operator, coefficients, output, status)
        !! Apply a caller-owned matrix independently at every Fourier mode.
        complex(dp), intent(in) :: operator(:, :, :)
        complex(dp), intent(in) :: coefficients(:, :)
        complex(dp), intent(out) :: output(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: input_component, output_component, mode

        output = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_operator_inputs( &
                operator, coefficients, output, status)) return
        do mode = 1, size(coefficients, 2)
            do output_component = 1, size(operator, 1)
                do input_component = 1, size(operator, 2)
                    output(output_component, mode) = &
                        output(output_component, mode) + &
                        operator(output_component, input_component, mode)* &
                        coefficients(input_component, mode)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fourier_mode_linear_operator

    subroutine apply_fourier_mode_linear_operator_jvp( &
            operator, coefficients, operator_dot, coefficients_dot, &
            output_dot, status)
        !! Apply the exact product-rule tangent of the mode-wise map.
        complex(dp), intent(in) :: operator(:, :, :)
        complex(dp), intent(in) :: coefficients(:, :)
        complex(dp), intent(in) :: operator_dot(:, :, :)
        complex(dp), intent(in) :: coefficients_dot(:, :)
        complex(dp), intent(out) :: output_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: input_component, output_component, mode

        output_dot = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_operator_inputs( &
                operator, coefficients, output_dot, status)) return
        if (.not. validate_operator_direction( &
                operator_dot, coefficients_dot, output_dot, status)) return
        do mode = 1, size(coefficients, 2)
            do output_component = 1, size(operator, 1)
                do input_component = 1, size(operator, 2)
                    output_dot(output_component, mode) = &
                        output_dot(output_component, mode) + &
                        operator_dot(output_component, input_component, mode)* &
                        coefficients(input_component, mode) + &
                        operator(output_component, input_component, mode)* &
                        coefficients_dot(input_component, mode)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fourier_mode_linear_operator_jvp

    subroutine apply_fourier_mode_linear_operator_vjp( &
            operator, coefficients, output_bar, operator_bar, &
            coefficients_bar, status)
        !! Apply the real-adjoint of the complex mode-wise map.
        !!
        !! The pairing is Re(sum(conjg(a_bar)*a_dot)).  Consequently the
        !! coefficient cotangent is A^H output_bar and the matrix cotangent is
        !! output_bar*conjg(coefficients), mode by mode.
        complex(dp), intent(in) :: operator(:, :, :)
        complex(dp), intent(in) :: coefficients(:, :)
        complex(dp), intent(in) :: output_bar(:, :)
        complex(dp), intent(out) :: operator_bar(:, :, :)
        complex(dp), intent(out) :: coefficients_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: input_component, output_component, mode

        operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_operator_inputs( &
                operator, coefficients, output_bar, status)) return
        if (.not. all_finite_vector(output_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode linear-operator cotangent is not finite")
            return
        end if
        if (size(operator_bar, 1) /= size(operator, 1) .or. &
                size(operator_bar, 2) /= size(operator, 2) .or. &
                size(operator_bar, 3) /= size(operator, 3) .or. &
                size(coefficients_bar, 1) /= size(coefficients, 1) .or. &
                size(coefficients_bar, 2) /= size(coefficients, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode linear-operator VJP output has wrong shape")
            return
        end if
        do mode = 1, size(coefficients, 2)
            do input_component = 1, size(operator, 2)
                do output_component = 1, size(operator, 1)
                    coefficients_bar(input_component, mode) = &
                        coefficients_bar(input_component, mode) + &
                        conjg(operator(output_component, input_component, mode))* &
                        output_bar(output_component, mode)
                    operator_bar(output_component, input_component, mode) = &
                        output_bar(output_component, mode)* &
                        conjg(coefficients(input_component, mode))
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine apply_fourier_mode_linear_operator_vjp

    logical function validate_operator_inputs(operator, coefficients, output, status)
        complex(dp), intent(in) :: operator(:, :, :)
        complex(dp), intent(in) :: coefficients(:, :)
        complex(dp), intent(in) :: output(:, :)
        type(fortsparse_status_t), intent(out) :: status

        validate_operator_inputs = .false.
        if (size(operator, 1) <= 0 .or. size(operator, 2) <= 0 .or. &
                size(operator, 3) <= 0 .or. &
                size(coefficients, 1) /= size(operator, 2) .or. &
                size(coefficients, 2) /= size(operator, 3) .or. &
                size(output, 1) /= size(operator, 1) .or. &
                size(output, 2) /= size(operator, 3)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode linear-operator arrays have incompatible shapes")
            return
        end if
        if (.not. all_finite_matrix(operator) .or. &
                .not. all_finite_vector(coefficients)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode linear-operator inputs are not finite")
            return
        end if
        validate_operator_inputs = .true.
    end function validate_operator_inputs

    logical function validate_operator_direction( &
            operator_dot, coefficients_dot, output_dot, status)
        complex(dp), intent(in) :: operator_dot(:, :, :)
        complex(dp), intent(in) :: coefficients_dot(:, :)
        complex(dp), intent(in) :: output_dot(:, :)
        type(fortsparse_status_t), intent(inout) :: status

        validate_operator_direction = .false.
        if (size(operator_dot, 1) /= size(output_dot, 1) .or. &
                size(operator_dot, 2) /= size(coefficients_dot, 1) .or. &
                size(operator_dot, 3) /= size(output_dot, 2) .or. &
                size(coefficients_dot, 2) /= size(output_dot, 2)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode linear-operator JVP has incompatible shapes")
            return
        end if
        if (.not. all_finite_matrix(operator_dot) .or. &
                .not. all_finite_vector(coefficients_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "Fourier mode linear-operator increments are not finite")
            return
        end if
        validate_operator_direction = .true.
    end function validate_operator_direction

    logical function all_finite_matrix(values)
        complex(dp), intent(in) :: values(:, :, :)

        all_finite_matrix = all(scalar_is_finite(values))
    end function all_finite_matrix

    logical function all_finite_vector(values)
        complex(dp), intent(in) :: values(:, :)

        all_finite_vector = all(scalar_is_finite(values))
    end function all_finite_vector

    pure elemental logical function scalar_is_finite(value)
        complex(dp), intent(in) :: value

        scalar_is_finite = ieee_is_finite(real(value, dp)) .and. &
            ieee_is_finite(aimag(value))
    end function scalar_is_finite

end module fortfem_fourier_mode_linear_operator
