program test_equilibrated_flux_3d
    !! Patchwise equilibrated-flux bound in 3D: (1) an exact polynomial flux
    !! grad u (psi_a grad u in P_4^3, inside RT_4) on a conforming
    !! 8-tetrahedron mesh of the reference tetrahedron is recovered, so the bound vanishes to rounding (exercises
    !! the face orientation and permutation transforms between neighbouring
    !! tetrahedra); (2) for A = 1/(1+x), u = sin(pi x) sin(pi y) sin(pi z) and
    !! the non-Galerkin P1 interpolant on Kuhn cube meshes, the bound is a
    !! true upper bound of the energy error and decreases with h.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_equilibrated_flux_3d, only: equilibrated_flux_estimate_3d, &
        equilibrated_flux_estimate_3d_t, EQUILIBRATION_OK
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: vertices(:, :), xi(:), eta(:), zeta(:), w(:)
    integer, allocatable :: tetrahedra(:, :)
    type(equilibrated_flux_estimate_3d_t) :: est
    real(dp) :: error, previous, ratio, efficiency
    integer :: status, n, k, ncell
    !> Measured efficiency 2.25/1.79 (k=0) and 1.61/1.38 (k=1) on n = 2/4.
    real(dp), parameter :: cap(0:1) = [2.5_dp, 1.8_dp]

    call tetra_duffy_quadrature(8, xi, eta, zeta, w, status)
    call check_condition(status == 0, "quadrature")
    call bey_mesh(vertices, tetrahedra)
    call polynomial_case(4, est, status)
    call check_condition(status == EQUILIBRATION_OK, "RT_4 equilibration succeeds")
    ! Rounding of the RT_4 saddle-point solves; a wrong face orientation or
    ! permutation between neighbouring tetrahedra leaves an O(1e-4) defect.
    call check_condition(est%bound < 1.0e-8_dp, &
        "RT_4 equilibration reproduces an exact polynomial flux")
    print '(a, 4es10.2)', "polynomial flux bound, flux, osc, div ", est%bound, &
        est%flux_defect, est%oscillation, est%divergence_defect

    call tetra_duffy_quadrature(6, xi, eta, zeta, w, status)
    call check_condition(status == 0, "quadrature")
    do k = 0, 1
        previous = huge(1.0_dp)
        do n = 2, 4, 2
            call cube_mesh(n, vertices, tetrahedra)
            call interpolant_case(k, est, error, status)
            call check_condition(status == EQUILIBRATION_OK, "interpolant case")
            call check_condition(est%bound >= error, &
                "equilibrated bound is above the energy error")
            call check_condition(est%bound < previous, "bound decreases with h")
            call check_condition(est%bound/error < cap(k), &
                "efficiency index below the measured cap")
            ncell = size(tetrahedra, 2)
            efficiency = est%bound/error
            ratio = 0.0_dp
            if (previous < huge(1.0_dp)) ratio = previous/est%bound
            if (ratio > 0.0_dp) call check_condition(ratio > 1.8_dp, &
                "bound converges at the energy rate (ratio > 1.8 per halving)")
            previous = est%bound
            print '(a, i2, a, i3, a, i5, 4(a, es10.2))', "k", k, " n", n, &
                " cells", ncell, " error", error, " bound", est%bound, &
                " efficiency", efficiency, " ratio", ratio
        end do
    end do
    call check_summary("equilibrated_flux_3d")

contains

    subroutine bey_mesh(points, cells)
        !! Reference tetrahedron split into 8 tetrahedra (4 corner tetrahedra
        !! and the octahedron cut along the M12-M34 diagonal).
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        integer :: c(4, 8), t

        allocate (points(3, 10))
        points(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        points(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        points(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        points(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        points(:, 5) = 0.5_dp*(points(:, 1) + points(:, 2))
        points(:, 6) = 0.5_dp*(points(:, 1) + points(:, 3))
        points(:, 7) = 0.5_dp*(points(:, 1) + points(:, 4))
        points(:, 8) = 0.5_dp*(points(:, 2) + points(:, 3))
        points(:, 9) = 0.5_dp*(points(:, 2) + points(:, 4))
        points(:, 10) = 0.5_dp*(points(:, 3) + points(:, 4))
        c(:, 1) = [1, 5, 6, 7]
        c(:, 2) = [2, 5, 8, 9]
        c(:, 3) = [3, 6, 8, 10]
        c(:, 4) = [4, 7, 9, 10]
        c(:, 5) = [5, 10, 6, 8]
        c(:, 6) = [5, 10, 8, 9]
        c(:, 7) = [5, 10, 9, 7]
        c(:, 8) = [5, 10, 7, 6]
        allocate (cells(4, 8))
        do t = 1, 8
            cells(:, t) = c(:, t)
        end do
        call orient_positively(points, cells)
    end subroutine bey_mesh

    subroutine cube_mesh(n, points, cells)
        !! Unit cube, n^3 cubes, each cut into 6 tetrahedra (Kuhn): the path
        !! from the lower corner to the upper corner along each axis permutation.
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        integer, parameter :: perms(3, 6) = reshape([1, 2, 3, 1, 3, 2, 2, 1, 3, &
            2, 3, 1, 3, 1, 2, 3, 2, 1], [3, 6])
        integer :: i, j, k, p, c, l, d(3)
        integer :: v(0:3)

        allocate (points(3, (n + 1)**3), cells(4, 6*n**3))
        do k = 0, n
            do j = 0, n
                do i = 0, n
                    points(:, i + (n + 1)*(j + (n + 1)*k) + 1) = &
                        [real(i, dp), real(j, dp), real(k, dp)]/n
                end do
            end do
        end do
        c = 0
        do k = 0, n - 1
            do j = 0, n - 1
                do i = 0, n - 1
                    do p = 1, 6
                        d = [i, j, k]
                        v(0) = d(1) + (n + 1)*(d(2) + (n + 1)*d(3)) + 1
                        do l = 1, 3
                            d(perms(l, p)) = d(perms(l, p)) + 1
                            v(l) = d(1) + (n + 1)*(d(2) + (n + 1)*d(3)) + 1
                        end do
                        c = c + 1
                        cells(:, c) = v
                    end do
                end do
            end do
        end do
        call orient_positively(points, cells)
    end subroutine cube_mesh

    subroutine orient_positively(points, cells)
        real(dp), intent(in) :: points(:, :)
        integer, intent(inout) :: cells(:, :)
        integer :: t, tmp
        real(dp) :: jac(3, 3)

        do t = 1, size(cells, 2)
            jac = jacobian(points(:, cells(:, t)))
            if (det3(jac) < 0.0_dp) then
                tmp = cells(3, t)
                cells(3, t) = cells(4, t)
                cells(4, t) = tmp
            end if
        end do
    end subroutine orient_positively

    pure function jacobian(p) result(jac)
        real(dp), intent(in) :: p(3, 4)
        real(dp) :: jac(3, 3)

        jac(:, 1) = p(:, 2) - p(:, 1)
        jac(:, 2) = p(:, 3) - p(:, 1)
        jac(:, 3) = p(:, 4) - p(:, 1)
    end function jacobian

    pure real(dp) function det3(a) result(det)
        real(dp), intent(in) :: a(3, 3)

        det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
            - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
            + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
    end function det3

    pure function gradient_of_hat_values(jac, values) result(grad)
        !! Gradient of the affine interpolant with vertex values (values(2..4)
        !! - values(1)) on the cell with Jacobian jac, by Cramer's rule on
        !! J^T g = d.
        real(dp), intent(in) :: jac(3, 3), values(4)
        real(dp) :: grad(3), d(3), a(3, 3), b(3, 3)
        integer :: i

        d = [values(2) - values(1), values(3) - values(1), values(4) - values(1)]
        a = transpose(jac)
        do i = 1, 3
            b = a
            b(:, i) = d
            grad(i) = det3(b)/det3(a)
        end do
    end function gradient_of_hat_values

    subroutine sample(points, cells, t, q, x)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: cells(:, :), t, q
        real(dp), intent(out) :: x(3)
        real(dp) :: p(3, 4)

        p = points(:, cells(:, t))
        x = p(:, 1) + xi(q)*(p(:, 2) - p(:, 1)) + eta(q)*(p(:, 3) - p(:, 1)) &
            + zeta(q)*(p(:, 4) - p(:, 1))
    end subroutine sample

    subroutine polynomial_case(degree, res, stat)
        !! A = 1, u = xyz(1-x-y-z) vanishes on the reference tetrahedron, its
        !! gradient is a cubic polynomial, so psi_a grad u is in P_4^3 (inside
        !! RT_4) and f = -Laplace u = 2(xy + yz + zx) in P_2; u_h = u exactly.
        integer, intent(in) :: degree
        type(equilibrated_flux_estimate_3d_t), intent(out) :: res
        integer, intent(out) :: stat
        real(dp), allocatable :: a(:, :), f(:, :), flux(:, :, :), lower(:)
        real(dp) :: x(3)
        integer :: t, q, nt

        nt = size(tetrahedra, 2)
        allocate (a(size(w), nt), f(size(w), nt), flux(3, size(w), nt), lower(nt))
        a = 1.0_dp
        lower = 1.0_dp
        do t = 1, nt
            do q = 1, size(w)
                call sample(vertices, tetrahedra, t, q, x)
                f(q, t) = 2*(x(1)*x(2) + x(2)*x(3) + x(3)*x(1))
                flux(:, q, t) = [x(2)*x(3)*(1 - 2*x(1) - x(2) - x(3)), &
                    x(1)*x(3)*(1 - x(1) - 2*x(2) - x(3)), &
                    x(1)*x(2)*(1 - x(1) - x(2) - 2*x(3))]
            end do
        end do
        call equilibrated_flux_estimate_3d(vertices, tetrahedra, degree, xi, eta, &
            zeta, w, a, f, flux, lower, 1.0_dp, res, stat)
    end subroutine polynomial_case

    subroutine interpolant_case(degree, res, energy_error, stat)
        !! A = 1/(1+x) in [1/2, 1]; u = sin(pi x) sin(pi y) sin(pi z) and the
        !! P1 interpolant u_h. C_F = sqrt(2/3)/pi: the cube's Friedrichs constant
        !! 1/(pi sqrt 3) combined with ||grad v|| <= sqrt(2) |||v|||.
        integer, intent(in) :: degree
        type(equilibrated_flux_estimate_3d_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        real(dp), allocatable :: a(:, :), f(:, :), flux(:, :, :), lower(:)
        real(dp) :: x(3), p(3, 4), jac(3, 3), vals(4), grad_h(3), grad(3)
        real(dp) :: lap, coef, dcoef, det
        integer :: t, q, nt

        nt = size(tetrahedra, 2)
        allocate (a(size(w), nt), f(size(w), nt), flux(3, size(w), nt), lower(nt))
        energy_error = 0.0_dp
        do t = 1, nt
            p = vertices(:, tetrahedra(:, t))
            lower(t) = 1.0_dp/(1.0_dp + maxval(p(1, :)))
            jac = jacobian(p)
            det = det3(jac)
            vals = sin(pi*p(1, :))*sin(pi*p(2, :))*sin(pi*p(3, :))
            grad_h = gradient_of_hat_values(jac, vals)
            do q = 1, size(w)
                call sample(vertices, tetrahedra, t, q, x)
                coef = 1.0_dp/(1.0_dp + x(1))
                dcoef = -coef*coef
                a(q, t) = coef
                grad = pi*[cos(pi*x(1))*sin(pi*x(2))*sin(pi*x(3)), &
                    sin(pi*x(1))*cos(pi*x(2))*sin(pi*x(3)), &
                    sin(pi*x(1))*sin(pi*x(2))*cos(pi*x(3))]
                lap = -3*pi*pi*sin(pi*x(1))*sin(pi*x(2))*sin(pi*x(3))
                f(q, t) = -(dcoef*grad(1) + coef*lap)
                flux(:, q, t) = coef*grad_h
                energy_error = energy_error + w(q)*abs(det)*coef &
                    *sum((grad - grad_h)**2)
            end do
        end do
        energy_error = sqrt(energy_error)
        call equilibrated_flux_estimate_3d(vertices, tetrahedra, degree, xi, eta, &
            zeta, w, a, f, flux, lower, sqrt(2.0_dp/3.0_dp)/pi, res, stat)
    end subroutine interpolant_case

end program test_equilibrated_flux_3d
