module fortfem_tensor_power_split
    !! Neutral quadratic power split for a tensor constitutive block.
    !!
    !! For a vector x and a (possibly nonsymmetric) tensor K, this contract
    !! reports x^T Sym(K) x, x^T Skew(K) x, and x^T K x separately.  The skew
    !! contribution is identically zero for a real vector, which makes the
    !! split useful as an energy/dissipation oracle for Hall or gyrotropic
    !! terms without imposing a plasma closure.  All products are pointwise;
    !! callers own quadrature, geometry, and field-aligned pullbacks.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_generated_tensor_power_split, only: &
        generated_tensor_power_split
    use fortfem_generated_tensor_power_split_jvp, only: &
        generated_tensor_power_split_jvp
    use fortfem_generated_tensor_power_split_vjp, only: &
        generated_tensor_power_split_vjp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_tensor_power_split
    public :: evaluate_tensor_power_split_jvp
    public :: evaluate_tensor_power_split_vjp

contains

    pure subroutine evaluate_tensor_power_split( &
            tensor, vector, symmetric_power, skew_power, total_power, status)
        real(dp), intent(in) :: tensor(:, :), vector(:)
        real(dp), intent(out) :: symmetric_power, skew_power, total_power
        type(fortsparse_status_t), intent(out) :: status

        symmetric_power = 0.0_dp
        skew_power = 0.0_dp
        total_power = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "tensor power split received incompatible arrays")
        if (.not. valid_inputs(tensor, vector)) return
        call generated_tensor_power_split( &
            tensor(1, 1), tensor(1, 2), tensor(1, 3), &
            tensor(2, 1), tensor(2, 2), tensor(2, 3), &
            tensor(3, 1), tensor(3, 2), tensor(3, 3), &
            vector(1), vector(2), vector(3), symmetric_power, skew_power, &
            total_power)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_tensor_power_split

    pure subroutine evaluate_tensor_power_split_jvp( &
            tensor, vector, tensor_dot, vector_dot, symmetric_power_dot, &
            skew_power_dot, total_power_dot, status)
        real(dp), intent(in) :: tensor(:, :), vector(:)
        real(dp), intent(in) :: tensor_dot(:, :), vector_dot(:)
        real(dp), intent(out) :: symmetric_power_dot, skew_power_dot
        real(dp), intent(out) :: total_power_dot
        type(fortsparse_status_t), intent(out) :: status

        symmetric_power_dot = 0.0_dp
        skew_power_dot = 0.0_dp
        total_power_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "tensor power split JVP received incompatible arrays")
        if (.not. valid_inputs(tensor, vector) .or. &
            .not. valid_inputs(tensor_dot, vector_dot)) return
        call generated_tensor_power_split_jvp( &
            tensor(1, 1), tensor(1, 2), tensor(1, 3), &
            tensor(2, 1), tensor(2, 2), tensor(2, 3), &
            tensor(3, 1), tensor(3, 2), tensor(3, 3), &
            vector(1), vector(2), vector(3), tensor_dot(1, 1), &
            tensor_dot(1, 2), tensor_dot(1, 3), tensor_dot(2, 1), &
            tensor_dot(2, 2), tensor_dot(2, 3), tensor_dot(3, 1), &
            tensor_dot(3, 2), tensor_dot(3, 3), vector_dot(1), vector_dot(2), &
            vector_dot(3), symmetric_power_dot, skew_power_dot, total_power_dot)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_tensor_power_split_jvp

    pure subroutine evaluate_tensor_power_split_vjp( &
            tensor, vector, symmetric_power_bar, skew_power_bar, total_power_bar, &
            tensor_bar, vector_bar, status)
        real(dp), intent(in) :: tensor(:, :), vector(:)
        real(dp), intent(in) :: symmetric_power_bar, skew_power_bar
        real(dp), intent(in) :: total_power_bar
        real(dp), intent(out) :: tensor_bar(:, :), vector_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        tensor_bar = 0.0_dp
        vector_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "tensor power split VJP received incompatible arrays")
        if (.not. valid_inputs(tensor, vector) .or. &
            size(tensor_bar, 1) /= 3 .or. size(tensor_bar, 2) /= 3 .or. &
            size(vector_bar) /= 3 .or. &
            .not. ieee_is_finite(symmetric_power_bar) .or. &
            .not. ieee_is_finite(skew_power_bar) .or. &
            .not. ieee_is_finite(total_power_bar)) return
        call generated_tensor_power_split_vjp( &
            tensor(1, 1), tensor(1, 2), tensor(1, 3), &
            tensor(2, 1), tensor(2, 2), tensor(2, 3), &
            tensor(3, 1), tensor(3, 2), tensor(3, 3), &
            vector(1), vector(2), vector(3), symmetric_power_bar, &
            skew_power_bar, total_power_bar, tensor_bar(1, 1), &
            tensor_bar(1, 2), tensor_bar(1, 3), tensor_bar(2, 1), &
            tensor_bar(2, 2), tensor_bar(2, 3), tensor_bar(3, 1), &
            tensor_bar(3, 2), tensor_bar(3, 3), vector_bar(1), vector_bar(2), &
            vector_bar(3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_tensor_power_split_vjp

    pure logical function valid_inputs(tensor, vector) result(valid)
        real(dp), intent(in) :: tensor(:, :), vector(:)

        valid = size(tensor, 1) == 3 .and. size(tensor, 2) == 3 .and. &
            size(vector) == 3 .and. all(ieee_is_finite(tensor)) .and. &
            all(ieee_is_finite(vector))
    end function valid_inputs

end module fortfem_tensor_power_split
