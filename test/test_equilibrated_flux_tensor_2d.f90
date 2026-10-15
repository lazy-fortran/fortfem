program test_equilibrated_flux_tensor_2d
    !! Equilibrated-flux bound with a symmetric positive definite tensor
    !! coefficient: (1) an exact polynomial flux A grad u inside RT_4 is
    !! recovered (bound at rounding level) for a constant anisotropic A;
    !! (2) A = a I reproduces the scalar estimator; (3) for the P1 Galerkin
    !! solution of a manufactured problem with rotated anisotropic A (ratio
    !! 1, 1e2, 1e4) the bound is above the energy error, decreases with h and
    !! has a bounded, non-increasing efficiency index; (4) a lower bound above the pointwise
    !! smaller eigenvalue is rejected.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_equilibrated_flux_2d, only: equilibrated_flux_estimate, &
        equilibrated_flux_estimate_tensor, equilibrated_flux_estimate_t, &
        EQUILIBRATION_OK, EQUILIBRATION_BAD_INPUT
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    !> Friedrichs constant of the unit square for the Laplacian.
    real(dp), parameter :: friedrichs = 1.0_dp/(sqrt(2.0_dp)*pi)
    real(dp), allocatable :: vertices(:, :), xi(:), eta(:), w(:)
    integer, allocatable :: triangles(:, :)
    type(equilibrated_flux_estimate_t) :: est, ref_est
    real(dp) :: a(2, 2), error, previous, index_value, ratios(3), caps(3)
    real(dp) :: previous_index
    integer :: status, n, k, r

    call triangle_duffy_quadrature(12, xi, eta, w, status)
    call check_condition(status == 0, "quadrature")

    call square_mesh(4, vertices, triangles)
    a = reshape([2.0_dp, 0.6_dp, 0.6_dp, 0.5_dp], [2, 2])
    call polynomial_case(4, a, est, status)
    call check_condition(status == EQUILIBRATION_OK, "tensor RT_4 equilibration")
    call check_condition(est%bound < 1.0e-9_dp, &
        "tensor RT_4 equilibration reproduces an exact polynomial flux")
    print '(a, 4es10.2)', "polynomial tensor flux bound, flux, osc, div ", &
        est%bound, est%flux_defect, est%oscillation, est%divergence_defect

    call square_mesh(8, vertices, triangles)
    do k = 0, 2
        call galerkin_case(k, 1.0_dp, 0.0_dp, .true., ref_est, error, status)
        call check_condition(status == EQUILIBRATION_OK, "scalar reference")
        call galerkin_case(k, 1.0_dp, 0.0_dp, .false., est, error, status)
        call check_condition(status == EQUILIBRATION_OK, "isotropic tensor")
        call check_condition(abs(est%bound - ref_est%bound) &
            <= 1.0e-12_dp*ref_est%bound, "A = a I reproduces the scalar bound")
    end do

    ! RT_k with k >= 1 for the P1 solution (RT_0 cannot represent the
    ! hat-weighted P1 flux; its index grows like the square root of the
    ! ratio through the A^(-1)-weighted defect).
    ratios = [1.0_dp, 1.0e2_dp, 1.0e4_dp]
    caps = [1.1_dp, 2.0_dp, 16.0_dp]
    do r = 1, size(ratios)
        do k = 1, 2
            previous = huge(1.0_dp)
            previous_index = huge(1.0_dp)
            do n = 8, 32, 8
                if (n == 24) cycle
                call square_mesh(n, vertices, triangles)
                call galerkin_case(k, ratios(r), pi/6, .false., est, error, status)
                call check_condition(status == EQUILIBRATION_OK, "anisotropic case")
                call check_condition(est%bound >= error, &
                    "anisotropic bound is above the energy error")
                call check_condition(est%bound < previous, "bound decreases with h")
                previous = est%bound
                index_value = est%bound/error
                ! The efficiency is not uniform in the anisotropy ratio (the
                ! hat-function localisation of the A^(-1)-weighted defect);
                ! measured: 1.05, 1.5-1.8, 8-15 at n = 8; it decreases with h.
                call check_condition(index_value < caps(r), &
                    "anisotropic efficiency index is bounded")
                call check_condition(index_value <= previous_index*(1 + 1.0e-3_dp), &
                    "anisotropic efficiency index does not grow with h")
                previous_index = index_value
                print '(a, es8.1, a, i2, a, i3, 3(a, es10.3), a, f7.3)', "ratio", &
                    ratios(r), " k", k, " n", n, " error", error, " bound", &
                    est%bound, " flux", est%flux_defect, " index", index_value
            end do
        end do
    end do

    call square_mesh(4, vertices, triangles)
    call rejected_lower_bound(status)
    call check_condition(status == EQUILIBRATION_BAD_INPUT, &
        "a lower bound above the smaller eigenvalue is rejected")
    call check_summary("equilibrated_flux_tensor_2d")

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

    pure real(dp) function smallest_eigenvalue(m) result(value)
        real(dp), intent(in) :: m(2, 2)

        value = 0.5_dp*(m(1, 1) + m(2, 2)) &
            - sqrt(0.25_dp*(m(1, 1) - m(2, 2))**2 + m(1, 2)**2)
    end function smallest_eigenvalue

    subroutine polynomial_case(degree, coef, res, stat)
        !! Constant A, u = x(1-x)y(1-y), exact flux A grad u (degree 3).
        integer, intent(in) :: degree
        real(dp), intent(in) :: coef(2, 2)
        type(equilibrated_flux_estimate_t), intent(out) :: res
        integer, intent(out) :: stat
        real(dp), allocatable :: tensor(:, :, :, :), f(:, :), flux(:, :, :), lower(:)
        real(dp) :: x(2), grad(2), hess(2, 2)
        integer :: t, q, nt

        nt = size(triangles, 2)
        allocate (tensor(2, 2, size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt))
        lower = smallest_eigenvalue(coef)
        do t = 1, nt
            do q = 1, size(w)
                call sample(t, q, x)
                tensor(:, :, q, t) = coef
                grad = [(1 - 2*x(1))*x(2)*(1 - x(2)), x(1)*(1 - x(1))*(1 - 2*x(2))]
                hess(1, 1) = -2*x(2)*(1 - x(2))
                hess(2, 2) = -2*x(1)*(1 - x(1))
                hess(1, 2) = (1 - 2*x(1))*(1 - 2*x(2))
                hess(2, 1) = hess(1, 2)
                f(q, t) = -sum(coef*hess)
                flux(:, q, t) = matmul(coef, grad)
            end do
        end do
        call equilibrated_flux_estimate_tensor(vertices, triangles, degree, xi, eta, &
            w, tensor, f, flux, lower, friedrichs/sqrt(lower(1)), res, stat)
    end subroutine polynomial_case

    subroutine galerkin_case(degree, ratio, angle, scalar, res, energy_error, stat)
        !! A = Q diag(1, 1/ratio) Q^T (Q the rotation by angle), u = sin(pi x)
        !! sin(pi y), u_h the P1 Galerkin solution with homogeneous Dirichlet
        !! data. With scalar (ratio 1 only) the scalar estimator is called.
        integer, intent(in) :: degree
        real(dp), intent(in) :: ratio, angle
        logical, intent(in) :: scalar
        type(equilibrated_flux_estimate_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        real(dp), allocatable :: tensor(:, :, :, :), f(:, :), flux(:, :, :)
        real(dp), allocatable :: lower(:), coefficient(:, :), u(:)
        real(dp) :: coef(2, 2), rot(2, 2), x(2), grad(2), hess(2, 2), e(2)
        real(dp) :: gradients(2, 3), grad_h(2), det
        integer :: t, q, nt

        rot = reshape([cos(angle), sin(angle), -sin(angle), cos(angle)], [2, 2])
        coef = matmul(rot, matmul(reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/ratio], &
            [2, 2]), transpose(rot)))
        coef(2, 1) = coef(1, 2)
        call p1_galerkin(coef, u)
        nt = size(triangles, 2)
        allocate (tensor(2, 2, size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt), coefficient(size(w), nt))
        lower = smallest_eigenvalue(coef)
        energy_error = 0.0_dp
        do t = 1, nt
            call p1_gradients(t, gradients, det)
            grad_h = matmul(gradients, u(triangles(:, t)))
            do q = 1, size(w)
                call sample(t, q, x)
                tensor(:, :, q, t) = coef
                coefficient(q, t) = coef(1, 1)
                grad = pi*[cos(pi*x(1))*sin(pi*x(2)), sin(pi*x(1))*cos(pi*x(2))]
                hess(1, 1) = -pi*pi*sin(pi*x(1))*sin(pi*x(2))
                hess(2, 2) = hess(1, 1)
                hess(1, 2) = pi*pi*cos(pi*x(1))*cos(pi*x(2))
                hess(2, 1) = hess(1, 2)
                f(q, t) = -sum(coef*hess)
                flux(:, q, t) = matmul(coef, grad_h)
                e = grad - grad_h
                energy_error = energy_error + w(q)*det*dot_product(e, matmul(coef, e))
            end do
        end do
        energy_error = sqrt(energy_error)
        if (scalar) then
            call equilibrated_flux_estimate(vertices, triangles, degree, xi, eta, &
                w, coefficient, f, flux, lower, friedrichs/sqrt(lower(1)), res, stat)
        else
            call equilibrated_flux_estimate_tensor(vertices, triangles, degree, xi, &
                eta, w, tensor, f, flux, lower, friedrichs/sqrt(lower(1)), res, stat)
        end if
    end subroutine galerkin_case

    subroutine p1_gradients(t, gradients, det)
        !! Physical gradients of the three P1 hat functions of cell t.
        integer, intent(in) :: t
        real(dp), intent(out) :: gradients(2, 3), det
        real(dp) :: p(2, 3), jac(2, 2)

        p = vertices(:, triangles(:, t))
        jac(:, 1) = p(:, 2) - p(:, 1)
        jac(:, 2) = p(:, 3) - p(:, 1)
        det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
        gradients(:, 2) = [jac(2, 2), -jac(1, 2)]/det
        gradients(:, 3) = [-jac(2, 1), jac(1, 1)]/det
        gradients(:, 1) = -gradients(:, 2) - gradients(:, 3)
    end subroutine p1_gradients

    subroutine p1_galerkin(coef, u)
        !! Dense P1 Galerkin solve of -div(coef grad u) = f, u = 0 on the
        !! boundary of the unit square; f from u = sin(pi x) sin(pi y).
        real(dp), intent(in) :: coef(2, 2)
        real(dp), allocatable, intent(out) :: u(:)
        real(dp), allocatable :: stiff(:, :), rhs(:)
        integer, allocatable :: ipiv(:)
        real(dp) :: gradients(2, 3), det, x(2), hat(3), source
        logical, allocatable :: boundary(:)
        integer :: nv, t, q, i, j, info
        interface
            subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
                import :: dp
                integer, intent(in) :: n, nrhs, lda, ldb
                integer, intent(out) :: ipiv(*), info
                real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            end subroutine dgesv
        end interface

        nv = size(vertices, 2)
        allocate (stiff(nv, nv), rhs(nv), ipiv(nv), boundary(nv))
        stiff = 0.0_dp
        rhs = 0.0_dp
        do t = 1, size(triangles, 2)
            call p1_gradients(t, gradients, det)
            do j = 1, 3
                do i = 1, 3
                    stiff(triangles(i, t), triangles(j, t)) = &
                        stiff(triangles(i, t), triangles(j, t)) + 0.5_dp*det &
                        *dot_product(gradients(:, i), matmul(coef, gradients(:, j)))
                end do
            end do
            do q = 1, size(w)
                call sample(t, q, x)
                hat = [1.0_dp - xi(q) - eta(q), xi(q), eta(q)]
                source = pi*pi*sin(pi*x(1))*sin(pi*x(2))*(coef(1, 1) + coef(2, 2)) &
                    - 2*pi*pi*coef(1, 2)*cos(pi*x(1))*cos(pi*x(2))
                rhs(triangles(:, t)) = rhs(triangles(:, t)) + w(q)*det*source*hat
            end do
        end do
        boundary = vertices(1, :) < 1.0e-12_dp .or. vertices(1, :) > 1 - 1.0e-12_dp
        boundary = boundary .or. vertices(2, :) < 1.0e-12_dp
        boundary = boundary .or. vertices(2, :) > 1 - 1.0e-12_dp
        do i = 1, nv
            if (.not. boundary(i)) cycle
            stiff(i, :) = 0.0_dp
            stiff(:, i) = 0.0_dp
            stiff(i, i) = 1.0_dp
            rhs(i) = 0.0_dp
        end do
        call dgesv(nv, 1, stiff, nv, ipiv, rhs, nv, info)
        call check_condition(info == 0, "P1 Galerkin solve")
        u = rhs
    end subroutine p1_galerkin

    subroutine rejected_lower_bound(stat)
        integer, intent(out) :: stat
        real(dp), allocatable :: tensor(:, :, :, :), f(:, :), flux(:, :, :), lower(:)
        type(equilibrated_flux_estimate_t) :: res
        integer :: nt

        nt = size(triangles, 2)
        allocate (tensor(2, 2, size(w), nt), f(size(w), nt), flux(2, size(w), nt))
        allocate (lower(nt))
        tensor(1, 1, :, :) = 1.0_dp
        tensor(2, 2, :, :) = 1.0e-3_dp
        tensor(1, 2, :, :) = 0.0_dp
        tensor(2, 1, :, :) = 0.0_dp
        f = 1.0_dp
        flux = 0.0_dp
        lower = 1.0e-3_dp
        lower(nt) = 2.0e-3_dp
        call equilibrated_flux_estimate_tensor(vertices, triangles, 1, xi, eta, w, &
            tensor, f, flux, lower, 1.0_dp, res, stat)
    end subroutine rejected_lower_bound

end program test_equilibrated_flux_tensor_2d
