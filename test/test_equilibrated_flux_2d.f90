program test_equilibrated_flux_2d
    !! Patchwise equilibrated-flux bound: (1) an exact polynomial flux inside
    !! RT_4 is recovered, so the bound vanishes to rounding (exercises edge
    !! orientation, moment signs and the interior-patch mean constraint);
    !! (2) for A = 1/(1+x) and the non-Galerkin P1 interpolant of a smooth
    !! solution the bound is a true upper bound of the energy error.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_t, EQUILIBRATION_OK
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: vertices(:, :), xi(:), eta(:), w(:)
    integer, allocatable :: triangles(:, :)
    type(equilibrated_flux_estimate_t) :: est
    real(dp) :: error, previous
    integer :: status, n, k

    call square_mesh(4, vertices, triangles)
    call triangle_duffy_quadrature(12, xi, eta, w, status)
    call check_condition(status == 0, "quadrature")
    call polynomial_case(4, est, status)
    call check_condition(status == EQUILIBRATION_OK, "RT_4 equilibration succeeds")
    ! Rounding of the RT_4 saddle-point solves leaves ||div sigma + Pi f|| ~ 1e-10
    ! against ||f|| ~ 0.5; a wrong edge sign or moment would give O(1e-2).
    call check_condition(est%bound < 1.0e-9_dp, &
        "RT_4 equilibration reproduces an exact polynomial flux")
    print '(a, 4es10.2)', "polynomial flux bound, flux, osc, div ", est%bound, &
        est%flux_defect, est%oscillation, est%divergence_defect

    do k = 0, 2
        previous = huge(1.0_dp)
        do n = 4, 16, 12
            call square_mesh(n, vertices, triangles)
            call interpolant_case(k, est, error, status)
            call check_condition(status == EQUILIBRATION_OK, "interpolant case")
            call check_condition(est%bound >= error, &
                "equilibrated bound is above the energy error")
            call check_condition(est%bound < previous, "bound decreases with h")
            previous = est%bound
            print '(a, i2, a, i3, 4(a, es10.2))', "k", k, " n", n, " error", &
                error, " bound", est%bound, " flux", est%flux_defect, &
                " divergence", est%divergence_defect
        end do
    end do
    call check_summary("equilibrated_flux_2d")

contains

    subroutine square_mesh(n, points, cells)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: points(:, :)
        integer, allocatable, intent(out) :: cells(:, :)
        integer :: i, j, v, c

        allocate (points(2, (n + 1)**2), cells(3, 2*n*n))
        do j = 0, n
            do i = 0, n
                points(:, j*(n + 1) + i + 1) = [real(i, dp), real(j, dp)]/n
            end do
        end do
        c = 0
        do j = 0, n - 1
            do i = 0, n - 1
                v = j*(n + 1) + i + 1
                ! Alternate diagonals so both edge orientations occur.
                if (mod(i + j, 2) == 0) then
                    cells(:, c + 1) = [v, v + 1, v + n + 2]
                    cells(:, c + 2) = [v, v + n + 2, v + n + 1]
                else
                    cells(:, c + 1) = [v, v + 1, v + n + 1]
                    cells(:, c + 2) = [v + 1, v + n + 2, v + n + 1]
                end if
                c = c + 2
            end do
        end do
    end subroutine square_mesh

    subroutine sample(cell, q, x)
        integer, intent(in) :: cell, q
        real(dp), intent(out) :: x(2)
        real(dp) :: p(2, 3)

        p = vertices(:, triangles(:, cell))
        x = p(:, 1) + xi(q)*(p(:, 2) - p(:, 1)) + eta(q)*(p(:, 3) - p(:, 1))
    end subroutine sample

    subroutine polynomial_case(degree, res, stat)
        !! A = 1, u = x(1-x)y(1-y), u_h = u (degree-4 exact flux).
        integer, intent(in) :: degree
        type(equilibrated_flux_estimate_t), intent(out) :: res
        integer, intent(out) :: stat
        real(dp), allocatable :: a(:, :), f(:, :), flux(:, :, :), lower(:)
        real(dp) :: x(2)
        integer :: t, q, nt

        nt = size(triangles, 2)
        allocate (a(size(w), nt), f(size(w), nt), flux(2, size(w), nt), lower(nt))
        a = 1.0_dp
        lower = 1.0_dp
        do t = 1, nt
            do q = 1, size(w)
                call sample(t, q, x)
                f(q, t) = 2*x(2)*(1 - x(2)) + 2*x(1)*(1 - x(1))
                flux(:, q, t) = [(1 - 2*x(1))*x(2)*(1 - x(2)), &
                    x(1)*(1 - x(1))*(1 - 2*x(2))]
            end do
        end do
        call equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, w, a, &
            f, flux, lower, 1.0_dp/pi, res, stat)
    end subroutine polynomial_case

    subroutine interpolant_case(degree, res, energy_error, stat)
        !! A = 1/(1+x), u = sin(pi x) sin(pi y), u_h its P1 interpolant.
        !! C_F: A >= 1/2 and the square's Friedrichs constant 1/(sqrt(2) pi).
        integer, intent(in) :: degree
        type(equilibrated_flux_estimate_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        real(dp), allocatable :: a(:, :), f(:, :), flux(:, :, :), lower(:)
        real(dp) :: x(2), p(2, 3), jac(2, 2), grad_h(2), grad(2), vals(3), det
        real(dp) :: lap, coef
        integer :: t, q, nt

        nt = size(triangles, 2)
        allocate (a(size(w), nt), f(size(w), nt), flux(2, size(w), nt), lower(nt))
        energy_error = 0.0_dp
        do t = 1, nt
            p = vertices(:, triangles(:, t))
            lower(t) = 1.0_dp/(1.0_dp + maxval(p(1, :)))
            jac(:, 1) = p(:, 2) - p(:, 1)
            jac(:, 2) = p(:, 3) - p(:, 1)
            det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
            vals = sin(pi*p(1, :))*sin(pi*p(2, :))
            ! grad u_h = J^(-T) [v2 - v1, v3 - v1]
            grad_h = [jac(2, 2)*(vals(2) - vals(1)) - jac(2, 1)*(vals(3) - vals(1)), &
                -jac(1, 2)*(vals(2) - vals(1)) + jac(1, 1)*(vals(3) - vals(1))]/det
            do q = 1, size(w)
                call sample(t, q, x)
                coef = 1.0_dp/(1.0_dp + x(1))
                a(q, t) = coef
                grad = pi*[cos(pi*x(1))*sin(pi*x(2)), sin(pi*x(1))*cos(pi*x(2))]
                lap = -2*pi*pi*sin(pi*x(1))*sin(pi*x(2))
                f(q, t) = -(-coef*coef*grad(1) + coef*lap)
                flux(:, q, t) = coef*grad_h
                energy_error = energy_error + w(q)*det*coef*sum((grad - grad_h)**2)
            end do
        end do
        energy_error = sqrt(energy_error)
        call equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, w, a, &
            f, flux, lower, 1.0_dp/pi, res, stat)
    end subroutine interpolant_case

end program test_equilibrated_flux_2d
