program test_geometry_jet_2d
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use fortfem_geometry_jet_2d, only: transform_scalar_jet_2d, &
        geometry_jet_success, geometry_jet_invalid_input, &
        geometry_jet_nonpositive_jacobian
    use basis_p2_2d_module, only: basis_p2_2d_t
    use check, only: check_condition, check_summary
    implicit none

    call test_affine_parity()
    call test_nonlinear_inverse()
    call test_invalid_samples()
    call check_summary('Sampled two-dimensional geometry jets')

contains

    subroutine test_affine_parity()
        type(basis_p2_2d_t) :: basis
        real(dp) :: vertices(2, 3), point(2), jacobian(2, 2), det_affine
        real(dp) :: map_hessian(2, 2, 2), reference_gradient(2)
        real(dp) :: reference_hessian(2, 2), inverse_jacobian(2, 2)
        real(dp) :: gradient(2), hessian(2, 2), determinant, coefficients(6)
        real(dp) :: expected_gradient(2), expected_hessian(2, 2), node(2)
        real(dp) :: xi, eta
        integer :: i, status

        vertices(:, 1) = [2.0_dp, -1.0_dp]
        vertices(:, 2) = [4.0_dp, -0.75_dp]
        vertices(:, 3) = [2.5_dp, 0.5_dp]
        xi = 0.25_dp
        eta = 0.125_dp
        call basis%compute_jacobian(vertices, jacobian, det_affine)
        call basis%transform_to_physical(xi, eta, vertices, point(1), point(2))
        do i = 1, 6
            call basis%transform_to_physical(basis%nodes(1, i), &
                basis%nodes(2, i), vertices, node(1), node(2))
            coefficients(i) = node(1)**2 + 3*node(1)*node(2) &
                + 2*node(2)**2 + 5*node(1) - 7*node(2) + 11
        end do
        reference_gradient = 0.0_dp
        reference_hessian = 0.0_dp
        do i = 1, 6
            reference_gradient = reference_gradient + &
                coefficients(i)*basis%grad(i, xi, eta)
            reference_hessian = reference_hessian + &
                coefficients(i)*basis%hessian(i, xi, eta)
        end do
        map_hessian = 0.0_dp
        call transform_scalar_jet_2d(point, jacobian, map_hessian, &
            reference_gradient, reference_hessian, determinant, &
            inverse_jacobian, gradient, hessian, status)
        expected_gradient = [2*point(1) + 3*point(2) + 5, &
            3*point(1) + 4*point(2) - 7]
        expected_hessian = reshape([2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [2, 2])
        call check_condition(status == geometry_jet_success, &
            'Affine geometry sample succeeds')
        call check_condition(abs(determinant - det_affine) < 1.0e-14_dp, &
            'Determinant agrees with existing affine basis geometry')
        call check_condition(maxval(abs(gradient - expected_gradient)) &
            < 1.0e-12_dp, 'Affine quadratic physical gradient')
        call check_condition(maxval(abs(hessian - expected_hessian)) &
            < 1.0e-12_dp, 'Affine quadratic physical Hessian')
        call check_condition(maxval(abs(matmul(jacobian, inverse_jacobian) &
            - reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) &
            < 1.0e-14_dp, 'Inverse uses the nonsymmetric Jacobian convention')
    end subroutine test_affine_parity

    subroutine test_nonlinear_inverse()
        real(dp) :: point(2), jacobian(2, 2), map_hessian(2, 2, 2)
        real(dp) :: reference_gradient(2), reference_hessian(2, 2)
        real(dp) :: inverse_jacobian(2, 2), gradient(2), hessian(2, 2)
        real(dp) :: determinant, expected_gradient(2), expected_hessian(2, 2)
        real(dp) :: wrong_hessian(2, 2), xi, eta, z
        integer :: sample, status

        ! T=(2+xi, eta+xi**2/2), inverse=(R-2, Z-(R-2)**2/2).
        ! The reference quadratic is xi**2+xi*eta+2*eta**2+3*xi-eta+1.
        ! Expand this with the explicit inverse and differentiate in R,Z.
        do sample = 1, 4
            xi = real(sample, dp)/8
            eta = 0.125_dp
            point = [2 + xi, eta + xi**2/2]
            z = point(2)
            jacobian = reshape([1.0_dp, xi, 0.0_dp, 1.0_dp], [2, 2])
            map_hessian = 0.0_dp
            map_hessian(1, 1, 2) = 1.0_dp
            reference_gradient = [2*xi + eta + 3, xi + 4*eta - 1]
            reference_hessian = &
                reshape([2.0_dp, 1.0_dp, 1.0_dp, 4.0_dp], [2, 2])
            call transform_scalar_jet_2d(point, jacobian, map_hessian, &
                reference_gradient, reference_hessian, determinant, &
                inverse_jacobian, gradient, hessian, status)
            expected_gradient = [z - 4*z*xi + 2*xi**3 - 1.5_dp*xi**2 &
                + 3*xi + 3, 4*z + xi - 2*xi**2 - 1]
            expected_hessian(1, 1) = -4*z + 6*xi**2 - 3*xi + 3
            expected_hessian(1, 2) = 1 - 4*xi
            expected_hessian(2, 1) = 1 - 4*xi
            expected_hessian(2, 2) = 4
            call check_condition(status == geometry_jet_success, &
                'Nonlinear sampled geometry succeeds')
            call check_condition(maxval(abs(gradient - expected_gradient)) &
                < 1.0e-14_dp, 'Explicit inverse physical gradient')
            call check_condition(maxval(abs(hessian - expected_hessian)) &
                < 1.0e-14_dp, 'Explicit inverse physical Hessian')
            wrong_hessian = matmul(transpose(inverse_jacobian), &
                matmul(reference_hessian, inverse_jacobian))
            if (sample < 4) then
                call check_condition(maxval(abs(wrong_hessian - hessian)) &
                    > 0.1_dp, 'Missing Hessian connection fails independent oracle')
            end if
            ! T_Z is a genuine reference P2 field representing physical Z.
            reference_gradient = [xi, 1.0_dp]
            reference_hessian = 0.0_dp
            reference_hessian(1, 1) = 1.0_dp
            call transform_scalar_jet_2d(point, jacobian, map_hessian, &
                reference_gradient, reference_hessian, determinant, &
                inverse_jacobian, gradient, hessian, status)
            call check_condition(maxval(abs(gradient - [0.0_dp, 1.0_dp])) &
                < 1.0e-14_dp, &
                'Represented nonlinear coordinate has physical unit gradient')
            call check_condition(maxval(abs(hessian)) < 1.0e-14_dp, &
                'Represented nonlinear coordinate has zero physical Hessian')
        end do
    end subroutine test_nonlinear_inverse

    subroutine test_invalid_samples()
        real(dp) :: point(2), jacobian(2, 2), map_hessian(2, 2, 2)
        real(dp) :: reference_gradient(2), reference_hessian(2, 2)
        real(dp) :: determinant, inverse_jacobian(2, 2), gradient(2), hessian(2, 2)
        real(dp) :: nan, infinity
        integer :: invalid_case, status, expected_status

        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        infinity = ieee_value(0.0_dp, ieee_positive_inf)
        do invalid_case = 1, 7
            point = [2.0_dp, 0.0_dp]
            jacobian = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
            map_hessian = 0.0_dp
            reference_gradient = [1.0_dp, 2.0_dp]
            reference_hessian = 0.0_dp
            expected_status = geometry_jet_invalid_input
            select case (invalid_case)
            case (1)
                jacobian(2, 2) = 0.0_dp
                expected_status = geometry_jet_nonpositive_jacobian
            case (2)
                jacobian(2, 2) = -1.0_dp
                expected_status = geometry_jet_nonpositive_jacobian
            case (3)
                point(1) = nan
            case (4)
                jacobian(1, 2) = nan
            case (5)
                map_hessian(1, 1, 2) = infinity
            case (6)
                reference_gradient(2) = nan
            case (7)
                reference_hessian(2, 2) = infinity
            end select
            call transform_scalar_jet_2d(point, jacobian, map_hessian, &
                reference_gradient, reference_hessian, determinant, &
                inverse_jacobian, gradient, hessian, status)
            call check_condition(status == expected_status, &
                'Invalid geometry/scalar sample rejected with declared status')
            call check_condition(determinant == 0.0_dp &
                .and. maxval(abs(inverse_jacobian)) == 0.0_dp &
                .and. maxval(abs(gradient)) == 0.0_dp &
                .and. maxval(abs(hessian)) == 0.0_dp, &
                'Invalid sample returns cleared finite outputs')
        end do
    end subroutine test_invalid_samples

end program test_geometry_jet_2d
