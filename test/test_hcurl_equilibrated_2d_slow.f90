program test_hcurl_equilibrated_2d_slow
    !! A posteriori energy bound for curl(nu curl u) + kappa u = f with u x n = 0
    !! on the order-one Nedelec space. (1) An exact field inside the space
    !! (constant curl, nu and kappa) gives a bound at rounding level. (2) For
    !! variable nu and kappa and a manufactured solution vanishing on the
    !! boundary, the bound dominates the exact energy error (quadrature oracle,
    !! Nedelec edge interpolant), stays under a measured efficiency cap and
    !! decays at the energy rate.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_nedelec_field_2d, only: evaluate_nedelec_field_2d
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_hcurl_equilibrated_2d, only: hcurl_equilibrated_estimate, &
        hcurl_equilibrated_estimate_t, HCURL_EQUILIBRATION_OK
    implicit none
    ! Measured efficiency indices (variable / constant coefficients):
    ! n = 8: 2.26 / 1.97, n = 16: 1.40 / 1.32, n = 32: 1.18 / 1.17,
    ! n = 64: 1.14; the caps leave a few percent of headroom. The bound
    ! decays at the energy rate.
    real(dp), parameter :: efficiency_cap(4) = [2.4_dp, 1.5_dp, 1.25_dp, 1.2_dp]
    real(dp), parameter :: rate_floor = 0.95_dp
    integer, parameter :: sizes(4) = [8, 16, 32, 64]

    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: vertices(:, :)
    type(mesh_2d_t) :: mesh
    type(hcurl_equilibrated_estimate_t) :: est
    real(dp), allocatable :: xi(:), eta(:), w(:), xo(:), eo(:), wo(:)
    real(dp) :: bound_prev, error, eff, rate
    integer :: status, i, n, previous_n, icase
    logical :: first

    call triangle_duffy_quadrature(12, xi, eta, w, status)
    call check_condition(status == 0, "estimator quadrature")
    call triangle_duffy_quadrature(16, xo, eo, wo, status)
    call check_condition(status == 0, "oracle quadrature")

    ! (1) Exact field u = (1, 2) + 0.9 (-y, x) lies in the Nedelec space.
    call square_mesh(4)
    call run_case(1, est, error, status)
    call check_condition(status == HCURL_EQUILIBRATION_OK, "exact case runs")
    print '(a, es10.2, a, es10.2)', "exact case: bound ", est%bound, &
        " error ", error
    call check_condition(est%bound < 1.0e-10_dp, &
        "exact Nedelec field gives a bound at rounding level")

    ! (2) Manufactured solution, variable (case 2) and constant (case 3)
    ! coefficients.
    do icase = 2, 3
    print '(a, i2)', "manufactured case ", icase
    first = .true.
    previous_n = 0
    bound_prev = huge(1.0_dp)
    do i = 1, size(sizes)
        n = sizes(i)
        call square_mesh(n)
        call run_case(icase, est, error, status)
        call check_condition(status == HCURL_EQUILIBRATION_OK, &
            "manufactured case runs")
        eff = est%bound/error
        call check_condition(est%bound >= error, &
            "bound dominates the exact energy error")
        call check_condition(eff <= efficiency_cap(i), "efficiency index within cap")
        if (.not. first) then
            rate = log(bound_prev/est%bound)/log(real(n, dp)/real(previous_n, dp))
            call check_condition(rate >= rate_floor, "bound decays at energy rate")
            print '(a, f6.2)', "  observed bound rate ", rate
        end if
        print '(a, i3, 4(a, es10.2))', "n", n, " error ", error, " bound ", &
            est%bound, " flux ", est%flux_defect, " field ", est%field_defect
        print '(a, f6.3)', "  efficiency index ", eff
        first = .false.
        previous_n = n
        bound_prev = est%bound
    end do
    end do
    call check_summary("hcurl_equilibrated_2d")

contains

    subroutine square_mesh(n)
        integer, intent(in) :: n
        integer :: i, j, v, c

        if (allocated(vertices)) deallocate (vertices, triangles)
        mesh = mesh_2d_t()
        allocate (vertices(2, (n + 1)**2), triangles(3, 2*n*n))
        do j = 0, n
            do i = 0, n
                vertices(:, j*(n + 1) + i + 1) = [real(i, dp), real(j, dp)]/n
            end do
        end do
        c = 0
        do j = 0, n - 1
            do i = 0, n - 1
                v = j*(n + 1) + i + 1
                ! Alternate diagonals so both edge orientations occur.
                if (mod(i + j, 2) == 0) then
                    triangles(:, c + 1) = [v, v + 1, v + n + 2]
                    triangles(:, c + 2) = [v, v + n + 2, v + n + 1]
                else
                    triangles(:, c + 1) = [v, v + 1, v + n + 1]
                    triangles(:, c + 2) = [v + 1, v + n + 2, v + n + 1]
                end if
                c = c + 2
            end do
        end do
        mesh%n_vertices = size(vertices, 2)
        mesh%n_triangles = size(triangles, 2)
        mesh%has_triangles = .true.
        allocate (mesh%vertices(2, mesh%n_vertices))
        allocate (mesh%triangles(3, mesh%n_triangles))
        mesh%vertices = vertices
        mesh%triangles = triangles
        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
    end subroutine square_mesh

    subroutine exact_fields(x, y, case, u, curl, nu, kappa, f)
        !! case 1: u = (1, 2) + 0.9 (-y, x), nu = 2, kappa = 1.5, f = kappa u.
        !! case 2: manufactured solution with variable nu and kappa.
        !! case 3: the same u with nu = 2 and kappa = 3.
        real(dp), intent(in) :: x, y
        integer, intent(in) :: case
        real(dp), intent(out) :: u(2), curl, nu, kappa, f(2)

        if (case == 1) then
            u = [1.0_dp - 0.9_dp*y, 2.0_dp + 0.9_dp*x]
            curl = 1.8_dp
            nu = 2.0_dp
            kappa = 1.5_dp
            f = kappa*u
            return
        end if
        ! Cases 2 and 3, from sympy: beta = x(1-x)y(1-y), u = beta(1 + y, x + y^2),
        ! nu = 1 + x/2, kappa = 1 + y, f = rot(nu curl u) + kappa u with
        ! rot(s) = (d_y s, -d_x s) (independent of the FortSym kernels).
        u(1) = x**2*y**3 - x**2*y - x*y**3 + x*y
        u(2) = x**3*y**2 - x**3*y + x**2*y**4 - x**2*y**3 - x**2*y**2 + x**2*y - &
              x*y**4 + x*y**3
        curl = -3*x**2*y + x**2 + 2*x*y**4 - 2*x*y**3 + x*y**2 + 2*x*y - x - y** &
              4 + y**3
        nu = (1.0d0/2.0d0)*x + 1
        kappa = y + 1
        if (case == 3) then
            ! f = nu rot(curl u) + kappa u with rot(c) = (d_y c, -d_x c).
            nu = 2.0_dp
            kappa = 3.0_dp
            f(1) = nu*(-3*x**2 + 8*x*y**3 - 6*x*y**2 + 2*x*y + 2*x - 4*y**3 &
                + 3*y**2) + kappa*u(1)
            f(2) = -nu*(-6*x*y + 2*x + 2*y**4 - 2*y**3 + y**2 + 2*y - 1) &
                + kappa*u(2)
            return
        end if
        f(1) = -3.0d0/2.0d0*x**3 + x**2*y**4 + 5*x**2*y**3 - 4*x**2*y**2 - 2*x** &
              2 - x*y**4 + 5*x*y**3 - 7.0d0/2.0d0*x*y**2 + 3*x*y + 2*x - 4*y**3 &
              + 3*y**2
        f(2) = x**3*y**3 - x**3*y + x**2*y**5 - 2*x**2*y**3 + (11.0d0/2.0d0)*x** &
              2*y - 3.0d0/2.0d0*x**2 - x*y**5 - 2*x*y**4 + 3*x*y**3 - x*y**2 + &
              4*x*y - x - 3.0d0/2.0d0*y**4 + (3.0d0/2.0d0)*y**3 - y**2 - 2*y + &
              1
    end subroutine exact_fields

    subroutine run_case(case, res, energy_error, stat)
        !! Nedelec edge interpolant: the global edge DOF is the tangential moment
        !! along the global orientation, i.e. orientation * (moment along the
        !! local edge direction from local vertex e to e + 1).
        integer, intent(in) :: case
        type(hcurl_equilibrated_estimate_t), intent(out) :: res
        real(dp), intent(out) :: energy_error
        integer, intent(out) :: stat
        complex(dp), allocatable :: dofs(:)
        complex(dp) :: value(2), curl_value
        real(dp), allocatable :: uh(:, :, :), curlh(:, :), nu(:, :), kap(:, :)
        real(dp), allocatable :: src(:, :, :), flux(:, :)
        real(dp) :: p0(2), d(2), mom, x, y, det, ue(2), ce, ne, ke, fe(2)
        real(dp) :: err_cell
        integer :: t, e, k, q, nt, nqt, edge_dofs(3), edge_orient(3)
        integer :: v1, v2, v3

        nt = size(triangles, 2)
        nqt = size(xi)
        allocate (dofs(mesh%n_edges))
        allocate (uh(2, nqt, nt), curlh(nqt, nt), nu(nqt, nt), kap(nqt, nt))
        allocate (src(2, nqt, nt), flux(nqt, nt))
        dofs = (0.0_dp, 0.0_dp)
        do t = 1, nt
            call mesh%get_triangle_edge_dofs(t, edge_dofs, edge_orient)
            do e = 1, 3
                v1 = triangles(e, t)
                v2 = triangles(mod(e, 3) + 1, t)
                p0 = vertices(:, v1)
                d = vertices(:, v2) - p0
                mom = 0.0_dp
                do k = 1, 5
                    x = p0(1) + 0.5_dp*(1.0_dp + gauss_x(k))*d(1)
                    y = p0(2) + 0.5_dp*(1.0_dp + gauss_x(k))*d(2)
                    call exact_fields(x, y, case, ue, ce, ne, ke, fe)
                    mom = mom + 0.5_dp*gauss_w(k)*dot_product(ue, d)
                end do
                dofs(edge_dofs(e) + 1) = cmplx(real(edge_orient(e), dp)*mom, &
                    0.0_dp, dp)
            end do
        end do

        do t = 1, nt
            do q = 1, nqt
                x = vertex_point(t, xi(q), eta(q), 1)
                y = vertex_point(t, xi(q), eta(q), 2)
                call exact_fields(x, y, case, ue, ce, ne, ke, fe)
                call evaluate_nedelec_field_2d(mesh, t, x, y, dofs, value, curl_value)
                uh(:, q, t) = real(value)
                curlh(q, t) = real(curl_value)
                nu(q, t) = ne
                kap(q, t) = ke
                src(:, q, t) = fe
                flux(q, t) = ne*curlh(q, t)
            end do
        end do

        ! Exact energy error |||u - u_h||| on the oracle rule.
        energy_error = 0.0_dp
        do t = 1, nt
            det = cell_det(t)
            err_cell = 0.0_dp
            do q = 1, size(xo)
                x = vertex_point(t, xo(q), eo(q), 1)
                y = vertex_point(t, xo(q), eo(q), 2)
                call exact_fields(x, y, case, ue, ce, ne, ke, fe)
                call evaluate_nedelec_field_2d(mesh, t, x, y, dofs, value, curl_value)
                err_cell = err_cell + wo(q)*det*(ne*(ce - real(curl_value))**2 &
                    + ke*sum((ue - real(value))**2))
            end do
            energy_error = energy_error + err_cell
        end do
        energy_error = sqrt(energy_error)

        call hcurl_equilibrated_estimate(vertices, triangles, xi, eta, w, nu, kap, &
            src, uh, flux, res, stat)
    end subroutine run_case

    real(dp) function vertex_point(t, xr, er, comp) result(xp)
        integer, intent(in) :: t, comp
        real(dp), intent(in) :: xr, er
        real(dp) :: p(2, 3)

        p = vertices(:, triangles(:, t))
        xp = p(comp, 1) + xr*(p(comp, 2) - p(comp, 1)) + er*(p(comp, 3) - p(comp, 1))
    end function vertex_point

    real(dp) function cell_det(t) result(det)
        integer, intent(in) :: t
        real(dp) :: p(2, 3)

        p = vertices(:, triangles(:, t))
        det = (p(1, 2) - p(1, 1))*(p(2, 3) - p(2, 1)) &
            - (p(1, 3) - p(1, 1))*(p(2, 2) - p(2, 1))
    end function cell_det

    real(dp) function gauss_x(k) result(xk)
        integer, intent(in) :: k
        real(dp), parameter :: nodes(5) = [-0.9061798459386640_dp, &
            -0.5384693101056831_dp, 0.0_dp, 0.5384693101056831_dp, &
            0.9061798459386640_dp]

        xk = nodes(k)
    end function gauss_x

    real(dp) function gauss_w(k) result(wk)
        integer, intent(in) :: k
        real(dp), parameter :: wts(5) = [0.2369268850561891_dp, &
            0.4786286704993665_dp, 0.5688888888888889_dp, &
            0.4786286704993665_dp, 0.2369268850561891_dp]

        wk = wts(k)
    end function gauss_w

end program test_hcurl_equilibrated_2d_slow
