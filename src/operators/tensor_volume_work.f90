module fortfem_tensor_volume_work
    !! Neutral tensor stress work/load contraction over volume quadrature.
    !!
    !! The residual is
    !!
    !!   r_i = sum_q w_q sum_{a,b} grad(v_i)_{a,b}(q) P_{a,b}(q).
    !!
    !! `P` may be a CGL pressure tensor, Maxwell stress, elastic stress, or
    !! any caller-owned anisotropic tensor.  This module deliberately does
    !! not choose a constitutive law or identify a physical field.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: assemble_tensor_volume_work
    public :: assemble_tensor_volume_work_jvp
    public :: assemble_tensor_volume_work_vjp

contains

    subroutine assemble_tensor_volume_work( &
            test_gradient, stress, weights, residual, status)
        real(dp), intent(in) :: test_gradient(:, :, :, :), stress(:, :, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: residual(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, test_function, row, column

        residual = 0.0_dp
        if (.not. validate_inputs( &
            test_gradient, stress, weights, residual, status)) return
        do quadrature = 1, size(weights)
            do test_function = 1, size(test_gradient, 2)
                do row = 1, 3
                    do column = 1, 3
                        residual(test_function) = residual(test_function) + &
                            weights(quadrature)*test_gradient( &
                            quadrature, test_function, row, column)* &
                            stress(quadrature, row, column)
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_tensor_volume_work

    subroutine assemble_tensor_volume_work_jvp( &
            test_gradient, stress, weights, test_gradient_dot, stress_dot, &
            weights_dot, residual_dot, status)
        real(dp), intent(in) :: test_gradient(:, :, :, :), stress(:, :, :)
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: test_gradient_dot(:, :, :, :), stress_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        real(dp), intent(out) :: residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, test_function, row, column

        residual_dot = 0.0_dp
        if (.not. validate_inputs( &
            test_gradient, stress, weights, residual_dot, status)) return
        if (.not. validate_direction( &
            test_gradient_dot, stress_dot, weights_dot, size(weights), &
            size(test_gradient, 2))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "tensor volume-work JVP has incompatible increments")
            return
        end if
        do quadrature = 1, size(weights)
            do test_function = 1, size(test_gradient, 2)
                do row = 1, 3
                    do column = 1, 3
                        residual_dot(test_function) = &
                            residual_dot(test_function) + &
                            weights_dot(quadrature)*test_gradient( &
                            quadrature, test_function, row, column)* &
                            stress(quadrature, row, column) + &
                            weights(quadrature)*test_gradient_dot( &
                            quadrature, test_function, row, column)* &
                            stress(quadrature, row, column) + &
                            weights(quadrature)*test_gradient( &
                            quadrature, test_function, row, column)* &
                            stress_dot(quadrature, row, column)
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_tensor_volume_work_jvp

    subroutine assemble_tensor_volume_work_vjp( &
            test_gradient, stress, weights, residual_bar, test_gradient_bar, &
            stress_bar, weights_bar, status)
        real(dp), intent(in) :: test_gradient(:, :, :, :), stress(:, :, :)
        real(dp), intent(in) :: weights(:), residual_bar(:)
        real(dp), intent(out) :: test_gradient_bar(:, :, :, :), stress_bar(:, :, :)
        real(dp), intent(out) :: weights_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: quadrature, test_function, row, column

        test_gradient_bar = 0.0_dp
        stress_bar = 0.0_dp
        weights_bar = 0.0_dp
        if (.not. validate_inputs( &
            test_gradient, stress, weights, residual_bar, status)) return
        if (.not. validate_adjoint( &
            test_gradient_bar, stress_bar, weights_bar, residual_bar, &
            size(weights), size(test_gradient, 2))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "tensor volume-work VJP has incompatible cotangents")
            return
        end if
        do quadrature = 1, size(weights)
            do test_function = 1, size(test_gradient, 2)
                do row = 1, 3
                    do column = 1, 3
                        test_gradient_bar(quadrature, test_function, row, column) = &
                            test_gradient_bar( &
                            quadrature, test_function, row, column) + &
                            weights(quadrature)*residual_bar(test_function)* &
                            stress(quadrature, row, column)
                        stress_bar(quadrature, row, column) = &
                            stress_bar(quadrature, row, column) + &
                            weights(quadrature)*residual_bar(test_function)* &
                            test_gradient(quadrature, test_function, row, column)
                        weights_bar(quadrature) = weights_bar(quadrature) + &
                            residual_bar(test_function)*test_gradient( &
                            quadrature, test_function, row, column)* &
                            stress(quadrature, row, column)
                    end do
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_tensor_volume_work_vjp

    logical function validate_inputs( &
            test_gradient, stress, weights, residual, status) result(valid)
        real(dp), intent(in) :: test_gradient(:, :, :, :), stress(:, :, :)
        real(dp), intent(in) :: weights(:), residual(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: quadrature_count, test_count

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "tensor volume work has incompatible arrays")
        quadrature_count = size(test_gradient, 1)
        test_count = size(test_gradient, 2)
        if (quadrature_count < 1 .or. test_count < 1 .or. &
            size(test_gradient, 3) /= 3 .or. size(test_gradient, 4) /= 3 .or. &
            size(stress, 1) /= quadrature_count .or. &
            size(stress, 2) /= 3 .or. size(stress, 3) /= 3 .or. &
            size(weights) /= quadrature_count .or. size(residual) /= test_count) return
        if (any(.not. ieee_is_finite(test_gradient)) .or. &
            any(.not. ieee_is_finite(stress)) .or. &
            any(.not. ieee_is_finite(weights)) .or. &
            any(.not. ieee_is_finite(residual)) .or. &
            any(weights <= 0.0_dp)) return
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_inputs

    logical function validate_direction( &
            test_gradient_dot, stress_dot, weights_dot, quadrature_count, &
            test_count) result(valid)
        real(dp), intent(in) :: test_gradient_dot(:, :, :, :), stress_dot(:, :, :)
        real(dp), intent(in) :: weights_dot(:)
        integer, intent(in) :: quadrature_count, test_count

        valid = size(test_gradient_dot, 1) == quadrature_count .and. &
            size(test_gradient_dot, 2) == test_count .and. &
            size(test_gradient_dot, 3) == 3 .and. &
            size(test_gradient_dot, 4) == 3 .and. &
            size(stress_dot, 1) == quadrature_count .and. &
            size(stress_dot, 2) == 3 .and. size(stress_dot, 3) == 3 .and. &
            size(weights_dot) == quadrature_count .and. &
            all(ieee_is_finite(test_gradient_dot)) .and. &
            all(ieee_is_finite(stress_dot)) .and. &
            all(ieee_is_finite(weights_dot))
    end function validate_direction

    logical function validate_adjoint( &
            test_gradient_bar, stress_bar, weights_bar, residual_bar, &
            quadrature_count, test_count) result(valid)
        real(dp), intent(in) :: test_gradient_bar(:, :, :, :), stress_bar(:, :, :)
        real(dp), intent(in) :: weights_bar(:), residual_bar(:)
        integer, intent(in) :: quadrature_count, test_count

        valid = size(test_gradient_bar, 1) == quadrature_count .and. &
            size(test_gradient_bar, 2) == test_count .and. &
            size(test_gradient_bar, 3) == 3 .and. &
            size(test_gradient_bar, 4) == 3 .and. &
            size(stress_bar, 1) == quadrature_count .and. &
            size(stress_bar, 2) == 3 .and. size(stress_bar, 3) == 3 .and. &
            size(weights_bar) == quadrature_count .and. &
            size(residual_bar) == test_count
    end function validate_adjoint

end module fortfem_tensor_volume_work
