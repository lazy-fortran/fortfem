program test_complex_tensor_power
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    use fortfem_feec, only: &
        evaluate_complex_tensor_power, evaluate_complex_tensor_power_jvp, &
        evaluate_complex_tensor_power_vjp
    implicit none

    complex(dp), parameter :: tensor(3, 3) = reshape([ &
        cmplx(1.2_dp, 0.2_dp, dp), cmplx(-0.3_dp, 0.4_dp, dp), &
        cmplx(0.1_dp, -0.2_dp, dp), cmplx(0.5_dp, -0.1_dp, dp), &
        cmplx(1.7_dp, 0.3_dp, dp), cmplx(-0.4_dp, 0.2_dp, dp), &
        cmplx(-0.2_dp, 0.5_dp, dp), cmplx(0.6_dp, -0.3_dp, dp), &
        cmplx(0.9_dp, -0.4_dp, dp)], [3, 3])
    complex(dp), parameter :: vector(3) = [ &
        cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.6_dp, dp), &
        cmplx(0.3_dp, 0.5_dp, dp)]
    complex(dp), parameter :: tensor_dot(3, 3) = reshape([ &
        cmplx(0.1_dp, -0.05_dp, dp), cmplx(-0.2_dp, 0.03_dp, dp), &
        cmplx(0.04_dp, 0.08_dp, dp), cmplx(0.07_dp, 0.02_dp, dp), &
        cmplx(-0.06_dp, 0.09_dp, dp), cmplx(0.03_dp, -0.04_dp, dp), &
        cmplx(-0.05_dp, 0.01_dp, dp), cmplx(0.02_dp, 0.06_dp, dp), &
        cmplx(0.08_dp, -0.07_dp, dp)], [3, 3])
    complex(dp), parameter :: vector_dot(3) = [ &
        cmplx(-0.03_dp, 0.02_dp, dp), cmplx(0.05_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, -0.04_dp, dp)]
    complex(dp), parameter :: power_bar = cmplx(0.6_dp, -0.8_dp, dp)
    real(dp), parameter :: epsilon = 1.0e-7_dp

    complex(dp) :: power, power_dot, power_plus, power_minus
    complex(dp) :: tensor_bar(3, 3), vector_bar(3), tensor_bad(3, 3)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status, status_plus, status_minus

    call evaluate_complex_tensor_power(tensor, vector, power, status)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(power - independent_power(tensor, vector)) < 1.0e-14_dp, &
        "complex tensor power matches independent nested-loop oracle")

    call evaluate_complex_tensor_power_jvp( &
        tensor, vector, tensor_dot, vector_dot, power_dot, status)
    call evaluate_complex_tensor_power( &
        tensor + epsilon*tensor_dot, vector + epsilon*vector_dot, &
        power_plus, status_plus)
    call evaluate_complex_tensor_power( &
        tensor - epsilon*tensor_dot, vector - epsilon*vector_dot, &
        power_minus, status_minus)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        status_plus%code == FORTSPARSE_OK .and. &
        status_minus%code == FORTSPARSE_OK .and. &
        abs(power_dot - (power_plus - power_minus)/(2.0_dp*epsilon)) < 2.0e-8_dp, &
        "complex tensor power JVP matches central reassembly")

    call evaluate_complex_tensor_power_vjp( &
        tensor, vector, power_bar, tensor_bar, vector_bar, status)
    lhs = real(conjg(power_bar)*power_dot, dp)
    rhs = real(sum(conjg(tensor_bar)*tensor_dot) + &
        sum(conjg(vector_bar)*vector_dot), dp)
    call check_condition(status%code == FORTSPARSE_OK .and. &
        abs(lhs - rhs) < 2.0e-14_dp, &
        "complex tensor power VJP satisfies real-part adjoint identity")

    tensor_bad = tensor
    tensor_bad(2, 1) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0.0_dp, dp)
    call evaluate_complex_tensor_power(tensor_bad, vector, power, status)
    call check_condition(status%code /= FORTSPARSE_OK .and. &
        power == cmplx(0.0_dp, 0.0_dp, dp), &
        "complex tensor power rejects non-finite tensor entries")

    call check_summary("complex anisotropic tensor power")

contains

    pure function independent_power(matrix, value) result(result)
        complex(dp), intent(in) :: matrix(:, :), value(:)
        complex(dp) :: result
        integer :: i, j

        result = cmplx(0.0_dp, 0.0_dp, dp)
        do i = 1, size(value)
            do j = 1, size(value)
                result = result + conjg(value(i))*matrix(i, j)*value(j)
            end do
        end do
    end function independent_power

end program test_complex_tensor_power
