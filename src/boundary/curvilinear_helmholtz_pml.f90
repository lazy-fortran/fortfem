module fortfem_curvilinear_helmholtz_pml
    !! Full complex 3-by-3 coordinate-stretch tensors for curved PML layers.
    !!
    !! For a nonsingular complex stretch S, the scalar H1 tensors are
    !!
    !!   G = det(S) S^{-1} S^{-T},   M = det(S),
    !!
    !! while the covariant H(curl) curl tensor is
    !!
    !!   C = det(S)^{-1} S^T S.
    !!
    !! Transposes are intentionally not conjugate-transposes: the complex
    !! coordinate transform is holomorphic.  The API accepts a full matrix so
    !! curvilinear, IGA, and geometry-generated layers do not have to be
    !! approximated by diagonal Cartesian stretches.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: curvilinear_scalar_helmholtz_pml_coefficients
    public :: curvilinear_scalar_helmholtz_pml_coefficients_jvp
    public :: curvilinear_scalar_helmholtz_pml_coefficients_vjp
    public :: curvilinear_curl_curl_pml_coefficients
    public :: curvilinear_curl_curl_pml_coefficients_jvp
    public :: curvilinear_curl_curl_pml_coefficients_vjp

    interface finite_complex
        module procedure finite_complex_scalar
        module procedure finite_complex_matrix
    end interface finite_complex

contains

    subroutine curvilinear_scalar_helmholtz_pml_coefficients( &
            stretch, gradient_coefficient, mass_coefficient, status)
        complex(dp), intent(in) :: stretch(3, 3)
        complex(dp), intent(out) :: gradient_coefficient(3, 3)
        complex(dp), intent(out) :: mass_coefficient
        integer, intent(out) :: status

        complex(dp) :: inverse(3, 3), determinant

        gradient_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        call inverse_and_determinant(stretch, inverse, determinant, status)
        if (status /= 0) return
        gradient_coefficient = determinant*matmul( &
            inverse, transpose(inverse))
        mass_coefficient = determinant
        status = 0
    end subroutine curvilinear_scalar_helmholtz_pml_coefficients

    subroutine curvilinear_scalar_helmholtz_pml_coefficients_jvp( &
            stretch, stretch_dot, gradient_coefficient_dot, &
            mass_coefficient_dot, status)
        complex(dp), intent(in) :: stretch(3, 3), stretch_dot(3, 3)
        complex(dp), intent(out) :: gradient_coefficient_dot(3, 3)
        complex(dp), intent(out) :: mass_coefficient_dot
        integer, intent(out) :: status

        complex(dp) :: inverse(3, 3), inverse_dot(3, 3), determinant
        complex(dp) :: determinant_dot

        gradient_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. finite_complex(stretch_dot)) return
        call inverse_and_determinant(stretch, inverse, determinant, status)
        if (status /= 0) return
        call determinant_and_cofactor_dot( &
            stretch, stretch_dot, determinant_dot)
        inverse_dot = -matmul(inverse, matmul(stretch_dot, inverse))
        gradient_coefficient_dot = determinant_dot*matmul( &
            inverse, transpose(inverse)) + determinant*matmul( &
            inverse_dot, transpose(inverse)) + determinant*matmul( &
            inverse, transpose(inverse_dot))
        mass_coefficient_dot = determinant_dot
        status = 0
    end subroutine curvilinear_scalar_helmholtz_pml_coefficients_jvp

    subroutine curvilinear_scalar_helmholtz_pml_coefficients_vjp( &
            stretch, gradient_coefficient_bar, mass_coefficient_bar, &
            stretch_bar, status)
        complex(dp), intent(in) :: stretch(3, 3), gradient_coefficient_bar(3, 3)
        complex(dp), intent(in) :: mass_coefficient_bar
        complex(dp), intent(out) :: stretch_bar(3, 3)
        integer, intent(out) :: status

        complex(dp) :: inverse(3, 3), determinant, p(3, 3), p_bar(3, 3)
        complex(dp) :: inverse_bar(3, 3), determinant_bar, cofactor(3, 3)

        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. finite_complex(gradient_coefficient_bar) .or. &
            .not. finite_complex(mass_coefficient_bar)) return
        call inverse_and_determinant(stretch, inverse, determinant, status)
        if (status /= 0) return
        p = matmul(inverse, transpose(inverse))
        determinant_bar = mass_coefficient_bar + sum( &
            gradient_coefficient_bar*conjg(p))
        p_bar = conjg(determinant)*gradient_coefficient_bar
        inverse_bar = matmul(p_bar, conjg(inverse)) + transpose(matmul( &
            conjg(transpose(inverse)), p_bar))
        stretch_bar = -matmul(conjg(transpose(inverse)), matmul( &
            inverse_bar, conjg(transpose(inverse))))
        call cofactor_matrix(stretch, cofactor)
        stretch_bar = stretch_bar + determinant_bar*conjg(cofactor)
        status = 0
    end subroutine curvilinear_scalar_helmholtz_pml_coefficients_vjp

    subroutine curvilinear_curl_curl_pml_coefficients( &
            stretch, curl_coefficient, mass_coefficient, status)
        complex(dp), intent(in) :: stretch(3, 3)
        complex(dp), intent(out) :: curl_coefficient(3, 3)
        complex(dp), intent(out) :: mass_coefficient(3, 3)
        integer, intent(out) :: status

        complex(dp) :: inverse(3, 3), determinant

        curl_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        call inverse_and_determinant(stretch, inverse, determinant, status)
        if (status /= 0) return
        curl_coefficient = matmul(transpose(stretch), stretch)/determinant
        mass_coefficient = determinant*matmul( &
            inverse, transpose(inverse))
        status = 0
    end subroutine curvilinear_curl_curl_pml_coefficients

    subroutine curvilinear_curl_curl_pml_coefficients_jvp( &
            stretch, stretch_dot, curl_coefficient_dot, &
            mass_coefficient_dot, status)
        complex(dp), intent(in) :: stretch(3, 3), stretch_dot(3, 3)
        complex(dp), intent(out) :: curl_coefficient_dot(3, 3)
        complex(dp), intent(out) :: mass_coefficient_dot(3, 3)
        integer, intent(out) :: status

        complex(dp) :: inverse(3, 3), inverse_dot(3, 3), determinant
        complex(dp) :: determinant_dot, tensor(3, 3), tensor_dot(3, 3)

        curl_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        mass_coefficient_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. finite_complex(stretch_dot)) return
        call inverse_and_determinant(stretch, inverse, determinant, status)
        if (status /= 0) return
        call determinant_and_cofactor_dot( &
            stretch, stretch_dot, determinant_dot)
        inverse_dot = -matmul(inverse, matmul(stretch_dot, inverse))
        tensor = matmul(transpose(stretch), stretch)
        tensor_dot = matmul(transpose(stretch_dot), stretch) + &
            matmul(transpose(stretch), stretch_dot)
        curl_coefficient_dot = tensor_dot/determinant - &
            tensor*determinant_dot/(determinant*determinant)
        mass_coefficient_dot = determinant_dot*matmul( &
            inverse, transpose(inverse)) + determinant*matmul( &
            inverse_dot, transpose(inverse)) + determinant*matmul( &
            inverse, transpose(inverse_dot))
        status = 0
    end subroutine curvilinear_curl_curl_pml_coefficients_jvp

    subroutine curvilinear_curl_curl_pml_coefficients_vjp( &
            stretch, curl_coefficient_bar, mass_coefficient_bar, &
            stretch_bar, status)
        complex(dp), intent(in) :: stretch(3, 3), curl_coefficient_bar(3, 3)
        complex(dp), intent(in) :: mass_coefficient_bar(3, 3)
        complex(dp), intent(out) :: stretch_bar(3, 3)
        integer, intent(out) :: status

        complex(dp) :: inverse(3, 3), determinant
        complex(dp) :: tensor(3, 3), tensor_bar(3, 3), tensor_a_bar(3, 3)
        complex(dp) :: determinant_bar, q_bar, q, cofactor(3, 3)
        complex(dp) :: zero_scalar

        stretch_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        zero_scalar = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. finite_complex(curl_coefficient_bar) .or. &
            .not. finite_complex(mass_coefficient_bar)) return
        call inverse_and_determinant(stretch, inverse, determinant, status)
        if (status /= 0) return
        call curvilinear_scalar_helmholtz_pml_coefficients_vjp( &
            stretch, mass_coefficient_bar, zero_scalar, stretch_bar, status)
        if (status /= 0) return

        q = 1.0_dp/determinant
        tensor = matmul(transpose(stretch), stretch)
        tensor_bar = conjg(q)*curl_coefficient_bar
        q_bar = sum(curl_coefficient_bar*conjg(tensor))
        determinant_bar = -q_bar*conjg(1.0_dp/(determinant*determinant))
        tensor_a_bar = matmul(tensor_bar, transpose(conjg(stretch)))
        stretch_bar = stretch_bar + transpose(tensor_a_bar) + &
            matmul(conjg(stretch), tensor_bar)
        call cofactor_matrix(stretch, cofactor)
        stretch_bar = stretch_bar + determinant_bar*conjg(cofactor)
        status = 0
    end subroutine curvilinear_curl_curl_pml_coefficients_vjp

    subroutine inverse_and_determinant(matrix, inverse, determinant, status)
        complex(dp), intent(in) :: matrix(3, 3)
        complex(dp), intent(out) :: inverse(3, 3), determinant
        integer, intent(out) :: status

        determinant = determinant3(matrix)
        inverse = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        if (.not. finite_complex(matrix) .or. .not. finite_complex(determinant)) return
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(matrix)))**3) return
        call cofactor_matrix(matrix, inverse)
        inverse = transpose(inverse)/determinant
        status = 0
    end subroutine inverse_and_determinant

    subroutine determinant_and_cofactor_dot(matrix, direction, determinant_dot)
        complex(dp), intent(in) :: matrix(3, 3), direction(3, 3)
        complex(dp), intent(out) :: determinant_dot
        complex(dp) :: cofactor(3, 3)

        call cofactor_matrix(matrix, cofactor)
        determinant_dot = sum(cofactor*direction)
    end subroutine determinant_and_cofactor_dot

    pure function determinant3(matrix) result(determinant)
        complex(dp), intent(in) :: matrix(3, 3)
        complex(dp) :: determinant

        determinant = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2)) - matrix(1, 2)*( &
            matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))
    end function determinant3

    pure subroutine cofactor_matrix(matrix, cofactor)
        complex(dp), intent(in) :: matrix(3, 3)
        complex(dp), intent(out) :: cofactor(3, 3)

        cofactor(1, 1) = matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)
        cofactor(1, 2) = -(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1))
        cofactor(1, 3) = matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1)
        cofactor(2, 1) = -(matrix(1, 2)*matrix(3, 3) - matrix(1, 3)*matrix(3, 2))
        cofactor(2, 2) = matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1)
        cofactor(2, 3) = -(matrix(1, 1)*matrix(3, 2) - matrix(1, 2)*matrix(3, 1))
        cofactor(3, 1) = matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2)
        cofactor(3, 2) = -(matrix(1, 1)*matrix(2, 3) - matrix(1, 3)*matrix(2, 1))
        cofactor(3, 3) = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
    end subroutine cofactor_matrix

    logical function finite_complex_scalar(value) result(valid)
        complex(dp), intent(in) :: value

        valid = ieee_is_finite(real(value, dp)) .and. &
            ieee_is_finite(aimag(value))
    end function finite_complex_scalar

    logical function finite_complex_matrix(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex_matrix

end module fortfem_curvilinear_helmholtz_pml
