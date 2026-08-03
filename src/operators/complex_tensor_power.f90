module fortfem_complex_tensor_power
    !! Complex quadratic power for caller-owned anisotropic tensors.
    !!
    !! The value is ``p = x^H K x``.  This is a neutral frequency-domain
    !! contraction: callers provide the tensor, normalization, geometry, and
    !! constitutive law.  The VJP follows the real-part complex convention.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: evaluate_complex_tensor_power
    public :: evaluate_complex_tensor_power_jvp
    public :: evaluate_complex_tensor_power_vjp

contains

    subroutine evaluate_complex_tensor_power(tensor, vector, power, status)
        complex(dp), intent(in) :: tensor(:, :), vector(:)
        complex(dp), intent(out) :: power
        type(fortsparse_status_t), intent(out) :: status

        integer :: i, j

        power = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex tensor power received incompatible arrays")
        if (.not. valid_inputs(tensor, vector)) return
        do i = 1, size(vector)
            do j = 1, size(vector)
                power = power + conjg(vector(i))*tensor(i, j)*vector(j)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_complex_tensor_power

    subroutine evaluate_complex_tensor_power_jvp( &
            tensor, vector, tensor_dot, vector_dot, power_dot, status)
        complex(dp), intent(in) :: tensor(:, :), vector(:)
        complex(dp), intent(in) :: tensor_dot(:, :), vector_dot(:)
        complex(dp), intent(out) :: power_dot
        type(fortsparse_status_t), intent(out) :: status

        integer :: i, j

        power_dot = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex tensor power JVP received incompatible arrays")
        if (.not. valid_inputs(tensor, vector) .or. &
            .not. valid_inputs(tensor_dot, vector_dot)) return
        do i = 1, size(vector)
            do j = 1, size(vector)
                power_dot = power_dot + conjg(vector_dot(i))* &
                    tensor(i, j)*vector(j) + &
                    conjg(vector(i))*tensor_dot(i, j)*vector(j) + &
                    conjg(vector(i))*tensor(i, j)*vector_dot(j)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_complex_tensor_power_jvp

    subroutine evaluate_complex_tensor_power_vjp( &
            tensor, vector, power_bar, tensor_bar, vector_bar, status)
        complex(dp), intent(in) :: tensor(:, :), vector(:), power_bar
        complex(dp), intent(out) :: tensor_bar(:, :), vector_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        complex(dp) :: tensor_vector(size(vector)), tensor_adjoint_vector(size(vector))
        integer :: i, j

        tensor_bar = cmplx(0.0_dp, 0.0_dp, dp)
        vector_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "complex tensor power VJP received incompatible arrays")
        if (.not. ieee_is_finite(real(power_bar, dp)) .or. &
            .not. ieee_is_finite(aimag(power_bar))) return
        if (.not. valid_inputs(tensor, vector)) return
        if (size(tensor_bar, 1) /= size(vector) .or. &
            size(tensor_bar, 2) /= size(vector) .or. &
            size(vector_bar) /= size(vector)) return
        tensor_vector = cmplx(0.0_dp, 0.0_dp, dp)
        tensor_adjoint_vector = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, size(vector)
            do j = 1, size(vector)
                tensor_vector(i) = tensor_vector(i) + tensor(i, j)*vector(j)
                tensor_adjoint_vector(i) = tensor_adjoint_vector(i) + &
                    conjg(tensor(j, i))*vector(j)
                tensor_bar(i, j) = power_bar*vector(i)*conjg(vector(j))
            end do
        end do
        vector_bar = conjg(power_bar)*tensor_vector + &
            power_bar*tensor_adjoint_vector
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_complex_tensor_power_vjp

    pure logical function valid_inputs(tensor, vector) result(valid)
        complex(dp), intent(in) :: tensor(:, :), vector(:)

        valid = size(vector) > 0 .and. size(tensor, 1) == size(vector) .and. &
            size(tensor, 2) == size(vector) .and. &
            all(ieee_is_finite(real(tensor, dp))) .and. &
            all(ieee_is_finite(aimag(tensor))) .and. &
            all(ieee_is_finite(real(vector, dp))) .and. &
            all(ieee_is_finite(aimag(vector)))
    end function valid_inputs

end module fortfem_complex_tensor_power
