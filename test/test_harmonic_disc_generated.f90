program test_harmonic_disc_generated
    !! The generated local kernels of fortfem_harmonic_disc against the
    !! handwritten arithmetic they replaced: complex Cauchy terms, Nystrom
    !! self terms, the Jacobian determinant of the geometry jet, and the
    !! Bernstein interpolation constants against the former runtime dgesv;
    !! the boundary weight against the curvature of a circle and of an
    !! ellipse in complex arithmetic.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_generated_harmonic_disc_cauchy_terms, only: &
        harmonic_disc_cauchy_terms
    use fortfem_generated_harmonic_disc_nystrom_offdiag, only: &
        harmonic_disc_nystrom_offdiag
    use fortfem_generated_harmonic_disc_nystrom_self, only: &
        harmonic_disc_nystrom_self
    use fortfem_generated_harmonic_disc_interior_limit_self, only: &
        harmonic_disc_interior_limit_self
    use fortfem_generated_harmonic_disc_boundary_weight, only: &
        harmonic_disc_boundary_weight
    use fortfem_generated_spline_jacobian_determinant, only: &
        spline_jacobian_determinant
    use fortfem_generated_spline_bernstein_interpolation, only: &
        spline_bernstein_interpolation, spline_bernstein_max_degree
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

    call cauchy_terms()
    call nystrom_terms()
    call jacobian_determinant()
    call boundary_weight()
    call bernstein_constants()
    call check_summary("harmonic disc generated kernels")

contains

    subroutine cauchy_terms()
        complex(dp) :: y, z, w, q, p
        real(dp) :: qr, qi, pr, pi_, err
        integer :: i

        err = 0.0_dp
        do i = 1, 20
            y = cmplx(0.3_dp*i, 1.0_dp - 0.1_dp*i, dp)
            z = cmplx(0.1_dp, 0.2_dp*i - 1.0_dp, dp)
            w = cmplx(0.7_dp - 0.05_dp*i, 0.4_dp + 0.03_dp*i, dp)
            call harmonic_disc_cauchy_terms(real(y, dp), aimag(y), real(z, dp), &
                aimag(z), real(w, dp), aimag(w), qr, qi, pr, pi_)
            q = w/(y - z)
            p = w/(y - z)**2
            err = max(err, abs(cmplx(qr, qi, dp) - q)/abs(q), &
                abs(cmplx(pr, pi_, dp) - p)/abs(p))
        end do
        call check_condition(err < 1.0e-14_dp, &
            "Cauchy terms w/(y-z) and w/(y-z)^2 match complex arithmetic")
    end subroutine cauchy_terms

    subroutine boundary_weight()
        !! |y'| (1 + (scale kappa)^2)^(1/3): on a circle of radius a
        !! traversed at speed a b, kappa = 1/a; on the ellipse (2 cos t,
        !! sin t), kappa = 2/|y'|^3.
        complex(dp) :: dy, ddy
        real(dp) :: weight, err, t, a, b, kappa, scale
        integer :: i

        err = 0.0_dp
        do i = 1, 12
            t = 0.5_dp*i
            a = 0.2_dp + 0.3_dp*i
            b = 1.0_dp + 0.1_dp*i
            scale = 0.4_dp*i
            dy = a*b*iu*exp(iu*t)
            ddy = -a*b*b*exp(iu*t)
            call harmonic_disc_boundary_weight(real(dy, dp), aimag(dy), real(ddy, dp), &
                aimag(ddy), scale, weight)
            err = max(err, abs(weight/(a*b*(1.0_dp + (scale/a)**2)**(1.0_dp/3.0_dp)) - &
                1.0_dp))
            dy = cmplx(-2.0_dp*sin(t), cos(t), dp)
            ddy = cmplx(-2.0_dp*cos(t), -sin(t), dp)
            kappa = 2.0_dp/abs(dy)**3
            call harmonic_disc_boundary_weight(real(dy, dp), aimag(dy), real(ddy, dp), &
                aimag(ddy), scale, weight)
            err = max(err, abs(weight/(abs(dy)*(1.0_dp + (scale*kappa)**2)** &
                (1.0_dp/3.0_dp)) - 1.0_dp))
        end do
        call harmonic_disc_boundary_weight(0.6_dp, -0.8_dp, 3.0_dp, 5.0_dp, 0.0_dp, weight)
        call check_condition(err < 1.0e-14_dp .and. abs(weight - 1.0_dp) < 1.0e-15_dp, &
            "boundary weight |y'| (1 + (scale kappa)^2)^(1/3) on a circle and an ellipse")
    end subroutine boundary_weight

    subroutine nystrom_terms()
        complex(dp) :: w, yj, yk, dy, ddy, c, e
        real(dp) :: cr, ci, h, err, dmu
        integer :: i

        err = 0.0_dp
        h = 2.0_dp*pi/64.0_dp
        do i = 1, 10
            w = cmplx(0.5_dp + 0.1_dp*i, -0.3_dp*i, dp)
            yj = cmplx(1.0_dp*i, 0.2_dp, dp)
            yk = cmplx(-0.4_dp, 0.1_dp*i, dp)
            dy = cmplx(0.3_dp*i, 1.0_dp, dp)
            ddy = cmplx(-1.0_dp, 0.2_dp*i, dp)
            dmu = 0.37_dp*i - 1.0_dp
            call harmonic_disc_nystrom_offdiag(real(w, dp), aimag(w), real(yj, dp), &
                aimag(yj), real(yk, dp), aimag(yk), cr, ci)
            c = w/(yj - yk)/(2.0_dp*pi*iu)
            err = max(err, abs(cmplx(cr, ci, dp) - c)/abs(c))
            call harmonic_disc_nystrom_self(real(ddy, dp), aimag(ddy), real(dy, dp), &
                aimag(dy), h, cr, ci)
            c = ddy/(2.0_dp*dy)*h/(2.0_dp*pi*iu)
            err = max(err, abs(cmplx(cr, ci, dp) - c)/abs(c))
            call harmonic_disc_interior_limit_self(dmu, h, cr, ci)
            e = dmu*h/(2.0_dp*pi*iu)
            err = max(err, abs(cmplx(cr, ci, dp) - e)/abs(e))
        end do
        call check_condition(err < 1.0e-14_dp, &
            "Nystrom off-diagonal, self and interior-limit terms match")
    end subroutine nystrom_terms

    subroutine jacobian_determinant()
        real(dp) :: jet(2, 6), det, expected
        integer :: i

        jet = reshape([(0.1_dp*i*i - 0.3_dp*i, i = 1, 12)], [2, 6])
        call spline_jacobian_determinant(jet(1, 2), jet(2, 2), jet(1, 3), &
            jet(2, 3), det)
        expected = jet(1, 2)*jet(2, 3) - jet(1, 3)*jet(2, 2)
        call check_condition(abs(det - expected) <= 1.0e-15_dp*abs(expected), &
            "Jacobian determinant of the geometry jet")
    end subroutine jacobian_determinant

    subroutine bernstein_constants()
        !! The former runtime path: Chebyshev nodes, Bernstein matrix, dgesv
        !! against the identity.
        real(dp), allocatable :: nodes(:), inverse(:, :), matrix(:, :)
        real(dp), allocatable :: reference(:, :), a(:, :)
        integer, allocatable :: pivot(:)
        real(dp) :: node_err, inverse_err, residual
        integer :: degree, m, i, j, status

        call check_condition(spline_bernstein_max_degree >= 3, &
            "generated Bernstein constants cover the tested degrees")
        do degree = 1, spline_bernstein_max_degree
            m = 2*degree
            call spline_bernstein_interpolation(degree, nodes, inverse, status)
            call check_condition(status == 0 .and. size(nodes) == m .and. &
                all(shape(inverse) == m), "degree "//digit(degree)//": shapes")
            allocate(matrix(m, m), reference(m, m), a(m, m), pivot(m))
            node_err = 0.0_dp
            do i = 1, m
                node_err = max(node_err, abs(nodes(i) - (0.5_dp - 0.5_dp* &
                    cos(pi*(real(i, dp) - 0.5_dp)/real(m, dp)))))
                do j = 1, m
                    matrix(i, j) = binomial(m - 1, j - 1)*nodes(i)**(j - 1)* &
                        (1.0_dp - nodes(i))**(m - j)
                end do
            end do
            reference = 0.0_dp
            do i = 1, m
                reference(i, i) = 1.0_dp
            end do
            a = matrix
            call dgesv(m, m, a, m, pivot, reference, m, status)
            inverse_err = maxval(abs(inverse - reference))/maxval(abs(reference))
            residual = maxval(abs(matmul(inverse, matrix))- eye(m))
            call check_condition(node_err < 1.0e-15_dp, &
                "degree "//digit(degree)//": nodes are the Chebyshev nodes")
            call check_condition(inverse_err < 1.0e-8_dp .and. residual < 1.0e-8_dp, &
                "degree "//digit(degree)//": inverse reproduces the dgesv inverse")
            deallocate(matrix, reference, a, pivot)
        end do
        call spline_bernstein_interpolation(spline_bernstein_max_degree + 1, nodes, &
            inverse, status)
        call check_condition(status /= 0, "unsupported degree is reported")
    end subroutine bernstein_constants

    function eye(m) result(e)
        integer, intent(in) :: m
        real(dp) :: e(m, m)
        integer :: i

        e = 0.0_dp
        do i = 1, m
            e(i, i) = 1.0_dp
        end do
    end function eye

    function digit(i) result(s)
        integer, intent(in) :: i
        character(2) :: s

        write(s, "(i2)") i
    end function digit

    pure real(dp) function binomial(n, k) result(c)
        integer, intent(in) :: n, k
        integer :: i

        c = 1.0_dp
        do i = 1, k
            c = c*real(n - k + i, dp)/real(i, dp)
        end do
    end function binomial

end program test_harmonic_disc_generated
