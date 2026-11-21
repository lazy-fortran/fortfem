module spline_thb_problem
    !! Wall layer on the unit disc: u = A(r) G(phi) with A = exp(-(1 - r)/eps)
    !! and G = exp(kappa (cos(phi - phi0) - 1)); u - lap u = f with the exact
    !! flux du/dn on the wall.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: layer_kernel_t, exact, source

    real(dp), parameter, public :: center(2) = [0.3_dp, -0.2_dp]
    real(dp), parameter :: eps = 0.04_dp, kappa = 6.0_dp, phi0 = 0.3_dp

    type, extends(pointwise_energy_t) :: layer_kernel_t
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type layer_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(layer_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            f = self%f(first_point + q - 1)
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2 + jet(1, q)**2) - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [jet(1, q) - f, jet(2, q), jet(3, q)]
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        end do
        status = 0
    end subroutine evaluate

    pure subroutine exact(s, u, grad)
        real(dp), intent(in) :: s(2)
        real(dp), intent(out) :: u, grad(2)

        real(dp) :: d(2), r, phi, a, g, dg

        d = s - center
        r = max(norm2(d), 1.0e-12_dp)
        phi = atan2(d(2), d(1))
        a = exp(-(1.0_dp - r)/eps)
        g = exp(kappa*(cos(phi - phi0) - 1.0_dp))
        dg = -kappa*sin(phi - phi0)*g
        u = a*g
        grad = (a/eps)*g*d/r + a*dg/r*[-d(2), d(1)]/r
    end subroutine exact

    pure real(dp) function source(s) result(f)
        real(dp), intent(in) :: s(2)

        real(dp) :: d(2), r, phi, a, g, d2g, u, grad(2)

        d = s - center
        r = max(norm2(d), 1.0e-12_dp)
        phi = atan2(d(2), d(1))
        a = exp(-(1.0_dp - r)/eps)
        g = exp(kappa*(cos(phi - phi0) - 1.0_dp))
        d2g = (kappa**2*sin(phi - phi0)**2 - kappa*cos(phi - phi0))*g
        call exact(s, u, grad)
        f = u - (a/eps**2*g + a/eps*g/r + a*d2g/r**2)
    end function source

end module spline_thb_problem

program test_spline_thb
    !! THB refinement of the C0 O-grid with periodic annulus (gate 2 of
    !! doc/design/multipatch_spline_spaces.md): partition of unity, exact
    !! projection from coarser spaces, at most two levels per truncated
    !! function, and adaptive radial-then-isotropic refinement of a wall
    !! layer reaching a smaller H1 error than uniform knots with more DOFs.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_side_quadrature_t, spline_side_quadrature, &
        spline_point_matrix, spline_geometry_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use fortfem_spline_thb, only: spline_thb_t, spline_thb_init, spline_thb_refine, &
        spline_thb_mark, spline_project, THB_ISOTROPIC, THB_RADIAL
    use spline_thb_problem, only: layer_kernel_t, exact, source, center
    implicit none

    type(block_layout_t) :: layout
    type(spline_thb_t) :: thb
    type(spline_space_t) :: coarse, uniform, before
    type(spline_quadrature_t) :: quad
    real(dp), allocatable :: c0(:, :), c1(:, :), c2(:, :), u(:, :), eta(:)
    real(dp) :: uniform_error(3), error, worst
    integer, allocatable :: marked(:)
    integer :: status, step, level, uniform_dofs(3), direction

    call layout_ogrid(center, 1.0_dp, 0.3_dp, layout, status, &
        interface_fraction=0.6_dp)
    call spline_thb_init(layout, options(4), thb, status)
    call check_condition(status == 0, "THB hierarchy on the O-grid with annulus")
    coarse = thb%space
    allocate(c0(1, coarse%dof_count))
    call random_number(c0)

    call spline_thb_refine(thb, wall_cells(thb%space, 0.8_dp, 0.6_dp), status)
    call spline_thb_refine(thb, wall_cells(thb%space, 0.9_dp, 0.4_dp), status, &
        THB_RADIAL)
    call spline_thb_refine(thb, wall_cells(thb%space, 0.95_dp, 0.3_dp), status)
    call check_condition(status == 0 .and. size(thb%level) == 4, &
        "three levels of local refinement, the second radial")
    call spline_quadrature_build(thb%space, 4, quad, status)
    allocate(u(1, thb%space%dof_count), c1(1, thb%space%dof_count))
    u = 1.0_dp
    call check_condition(max_deviation(thb%space, quad, u, 1.0_dp) < 1.0e-12_dp, &
        "truncated functions form a partition of unity")
    call spline_project(coarse, c0, thb%space, quad, c1, status)
    worst = projection_defect(coarse, c0, thb%space, quad, c1)
    call check_condition(status == 0 .and. worst < 1.0e-10_dp, &
        "projection of a level-0 field into the hierarchy is exact")
    call check_condition(max_levels(thb%space) <= 2, &
        "admissible refinement: every function spans at most two levels")
    before = thb%space
    call spline_thb_refine(thb, wall_cells(thb%space, 0.97_dp, 0.2_dp), status)
    call spline_quadrature_build(thb%space, 4, quad, status)
    allocate(c2(1, thb%space%dof_count))
    call spline_project(before, c1, thb%space, quad, c2, status)
    call check_condition(status == 0 .and. projection_defect(before, c1, &
        thb%space, quad, c2) < 1.0e-10_dp, &
        "projection from the coarser hierarchy is exact (nested)")

    do level = 1, 3
        call spline_space_build(layout, options(2**(level + 1)), uniform, status)
        call solve(uniform, error, eta)
        uniform_error(level) = error
        uniform_dofs(level) = uniform%dof_count
    end do
    call spline_thb_init(layout, options(4), thb, status)
    do step = 1, 7
        call solve(thb%space, error, eta)
        print '(a,i2,a,i6,a,es10.3)', "THB step", step, " dofs", &
            thb%space%dof_count, " H1 error", error
        if (step == 7) exit
        call spline_thb_mark(eta, 0.5_dp, marked)
        direction = merge(THB_RADIAL, THB_ISOTROPIC, size(thb%level) <= 2)
        call spline_thb_refine(thb, marked, status, direction)
        if (status /= 0) error stop "refine"
    end do
    print '(a,3i7,a,3es10.3)', "uniform dofs", uniform_dofs, " H1 errors", &
        uniform_error
    call check_condition(thb%space%dof_count < uniform_dofs(3) .and. &
        error < uniform_error(3), &
        "adaptive THB beats uniform knots: smaller H1 error with fewer DOFs")
    call check_summary("spline_thb")

contains

    function options(cells) result(opt)
        integer, intent(in) :: cells
        type(spline_space_options_t) :: opt

        opt%degree = 3
        allocate(opt%knots(2, 6))
        opt%knots(:, :)%element_count = cells
    end function options

    function wall_cells(space, rmin, width) result(cells)
        !! Cells whose centre lies in the wall sector r > rmin, |phi| < width.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: rmin, width
        integer, allocatable :: cells(:)

        real(dp) :: jet(2, 6), xi(2), d(2)
        integer :: cell, patch, e(2), status

        allocate(cells(0))
        do cell = 1, space%cell_count
            patch = space%cell_block(cell)
            e = space%cell_element(:, cell)
            xi(1) = 0.5_dp*sum(space%blocks(patch)%breaks1(e(1):e(1) + 1))
            xi(2) = 0.5_dp*sum(space%blocks(patch)%breaks2(e(2):e(2) + 1))
            call spline_geometry_at(space, space%patch_block(patch), xi, jet, status)
            d = jet(:, 1) - center
            if (norm2(d) > rmin .and. abs(atan2(d(2), d(1))) < width) &
                cells = [cells, cell]
        end do
    end function wall_cells

    real(dp) function max_deviation(space, quad, coef, value) result(worst)
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: coef(:, :), value

        real(dp), allocatable :: jet(:, :)
        integer :: status

        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, coef, 0, jet, status)
        worst = maxval(abs(jet(1, :) - value))
    end function max_deviation

    real(dp) function projection_defect(source, cs, target, quad, ct) result(worst)
        use fortfem_spline_eval, only: spline_evaluate_at
        type(spline_space_t), intent(in) :: source, target
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: cs(:, :), ct(:, :)

        real(dp), allocatable :: jet(:, :)
        real(dp) :: v(1)
        integer :: q, status

        allocate(jet(1, quad%point_count))
        call spline_evaluate(target, quad, ct, 0, jet, status)
        worst = 0.0_dp
        do q = 1, quad%point_count
            call spline_evaluate_at(source, cs, quad%block(q), quad%xi(:, q), 0, v, &
                status)
            worst = max(worst, abs(v(1) - jet(1, q)))
        end do
    end function projection_defect

    integer function max_levels(space) result(span)
        type(spline_space_t), intent(in) :: space

        integer, allocatable :: lo(:), hi(:)
        integer :: j, q, patch, level

        allocate(lo(space%dof_count), hi(space%dof_count))
        lo = huge(1)
        hi = -1
        do patch = 1, space%block_count
            level = space%patch_level(patch)
            do j = space%blocks(patch)%offset + 1, space%blocks(patch)%offset + &
                    product(space%blocks(patch)%n)
                do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                    lo(space%broken_dof(q)) = min(lo(space%broken_dof(q)), level)
                    hi(space%broken_dof(q)) = max(hi(space%broken_dof(q)), level)
                end do
            end do
        end do
        span = maxval(hi - lo) + 1
    end function max_levels

    subroutine solve(space, error, eta)
        !! Discrete solution, its H1 error and the cellwise H1 errors.
        type(spline_space_t), intent(in) :: space
        real(dp), intent(out) :: error
        real(dp), allocatable, intent(out) :: eta(:)

        type(spline_quadrature_t) :: quad
        type(spline_side_quadrature_t) :: side
        type(layer_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :), val(:), b(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp) :: energy, normal(2), ue, ge(2), e2
        integer :: status, q, k, cell

        call spline_quadrature_build(space, 5, quad, status)
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = source(quad%s(:, q))
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        b = -g(1, :)
        call spline_side_quadrature(space, LAYOUT_WALL, 6, side, status)
        call spline_point_matrix(space, side%block, side%xi, 0, rowptr, colidx, val, &
            status)
        do q = 1, side%point_count
            normal = [side%tangent(2, q), -side%tangent(1, q)]
            call exact(side%s(:, q), ue, ge)
            do k = rowptr(q), rowptr(q + 1) - 1
                b(colidx(k)) = b(colidx(k)) + side%weight(q)*val(k)* &
                    dot_product(ge, normal)
            end do
        end do
        call pcg(matrix, b, c(1, :))
        allocate(jet(3, quad%point_count), eta(space%cell_count))
        call spline_evaluate(space, quad, c, 1, jet, status)
        eta = 0.0_dp
        do cell = 1, space%cell_count
            do q = quad%point_offset(cell), quad%point_offset(cell) + &
                    quad%points_per_cell - 1
                call exact(quad%s(:, q), ue, ge)
                e2 = (jet(1, q) - ue)**2 + sum((jet(2:3, q) - ge)**2)
                eta(cell) = eta(cell) + quad%weight(q)*e2
            end do
        end do
        error = sqrt(sum(eta))
        eta = sqrt(eta)
    end subroutine solve

    subroutine pcg(matrix, b, x)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)

        real(dp), allocatable :: r(:), z(:), p(:), ap(:), diag(:)
        real(dp) :: rz, rz_new, alpha, norm_b
        integer :: it, i, k

        allocate(r(size(b)), z(size(b)), p(size(b)), ap(size(b)), diag(size(b)))
        do i = 1, matrix%n
            diag(i) = 1.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                if (matrix%colidx(k) == i) diag(i) = matrix%val(k)
            end do
        end do
        x = 0.0_dp
        r = b
        norm_b = norm2(r)
        z = r/diag
        p = z
        rz = dot_product(r, z)
        do it = 1, 20*size(b)
            call matvec(matrix, p, ap)
            alpha = rz/dot_product(p, ap)
            x = x + alpha*p
            r = r - alpha*ap
            if (norm2(r) < 1.0e-12_dp*norm_b) exit
            z = r/diag
            rz_new = dot_product(r, z)
            p = z + (rz_new/rz)*p
            rz = rz_new
        end do
    end subroutine pcg

    subroutine matvec(matrix, v, y)
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer :: row, k

        do row = 1, matrix%n
            y(row) = 0.0_dp
            do k = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                y(row) = y(row) + matrix%val(k)*v(matrix%colidx(k))
            end do
        end do
    end subroutine matvec

end program test_spline_thb
