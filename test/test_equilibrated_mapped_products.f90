program test_equilibrated_mapped_products
    ! Independent polynomial oracle for whole-cell determinant hulls. The
    ! Bernstein determinant is evaluated separately from det(evaluate(J)).
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_generated_equilibrated_determinant_controls, only: &
        generated_equilibrated_determinant_controls
    use fortfem_generated_equilibrated_relative_jacobian_norm, only: &
        generated_equilibrated_relative_jacobian_norm
    implicit none
    real(dp) :: jac(2, 2, 6), controls(15), point(3), weights(6), matrix(2, 2)
    real(dp) :: direct, polynomial, coefficient, max_error, norm, invt(2, 2), angle
    integer :: i, j, a, b, c, n

    do n = 1, 6
        jac(:, :, n) = reshape([1.0_dp + n/20.0_dp, sin(real(n, dp))/10, &
            cos(real(n, dp))/8, 0.8_dp + n/30.0_dp], [2, 2])
    end do
    call generated_equilibrated_determinant_controls(jac, controls)
    max_error = 0
    do i = 0, 10
        do j = 0, 10 - i
            point = [real(i, dp)/10, real(j, dp)/10, real(10-i-j, dp)/10]
            weights = [point**2, 2*point(1)*point(2), &
                2*point(2)*point(3), 2*point(1)*point(3)]
            matrix = 0
            do n = 1, 6
                matrix = matrix + weights(n)*jac(:, :, n)
            end do
            direct = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
            polynomial = 0
            n = 0
            do a = 0, 4
                do b = 0, 4 - a
                    c = 4 - a - b
                    n = n + 1
                    coefficient = real(factorial(4), dp) &
                        /(factorial(a)*factorial(b)*factorial(c))
                    polynomial = polynomial + controls(n)*coefficient &
                        *point(1)**a*point(2)**b*point(3)**c
                end do
            end do
            max_error = max(max_error, abs(polynomial - direct))
            call check_condition(direct >= minval(controls) - 1.0e-14_dp &
                .and. direct <= maxval(controls) + 1.0e-14_dp, &
                'determinant lies in whole-cell Bernstein hull')
        end do
    end do
    call check_condition(max_error < 2.0e-14_dp, 'determinant polynomial identity')

    ! A rotation times diag(1,2) has spectral norm exactly two. A nonsymmetric
    ! inverse transpose catches transposition mistakes in J J0^(-1).
    angle = 0.37_dp
    matrix = reshape([cos(angle), sin(angle), -2*sin(angle), 2*cos(angle)], [2, 2])
    invt = reshape([1.0_dp, 0.0_dp, 0.3_dp, 1.0_dp], [2, 2])
    jac(:, :, 1) = matmul(matrix, reshape([1.0_dp, -0.3_dp, 0.0_dp, 1.0_dp], &
        [2, 2]))
    call generated_equilibrated_relative_jacobian_norm(jac(:, :, 1), invt, norm)
    call check_condition(abs(norm - 2.0_dp) < 1.0e-14_dp, 'spectral norm oracle')
    call check_summary('equilibrated_mapped_products')
contains
    integer function factorial(k) result(value)
        integer, intent(in) :: k
        integer :: m
        value = 1
        do m = 2, k
            value = value*m
        end do
    end function factorial
end program test_equilibrated_mapped_products
