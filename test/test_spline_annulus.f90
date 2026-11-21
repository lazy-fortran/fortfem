module spline_annulus_kernel
    !! W = |grad u|^2/2 + u^2/2 - f u on the disc |s - c| <= 1, so the
    !! minimiser solves u - lap u = f with zero normal derivative. Exact
    !! solution u = cos(pi r) + (r - r^3/3) cos(theta), du/dr = 0 at r = 1.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: annulus_kernel_t, exact_solution, source

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter, public :: center(2) = [0.3_dp, -0.2_dp]

    type, extends(pointwise_energy_t) :: annulus_kernel_t
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type annulus_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(annulus_kernel_t), intent(in) :: self
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

    pure real(dp) function exact_solution(s) result(u)
        real(dp), intent(in) :: s(2)

        real(dp) :: d(2), r

        d = s - center
        r = norm2(d)
        u = cos(pi*r) + (1.0_dp - r**2/3.0_dp)*d(1)
    end function exact_solution

    pure real(dp) function source(s) result(f)
        !! f = u - lap u with lap cos(pi r) = -pi^2 cos(pi r) - pi sin(pi r)/r
        !! and lap(g(r) cos(theta)) = (g'' + g'/r - g/r^2) cos(theta).
        real(dp), intent(in) :: s(2)

        real(dp) :: d(2), r, lap

        d = s - center
        r = norm2(d)
        if (r < 1.0e-12_dp) then
            lap = -2.0_dp*pi**2
        else
            lap = -pi**2*cos(pi*r) - pi*sin(pi*r)/r
        end if
        ! g = r - r^3/3: g'' + g'/r - g/r^2 = -8 r/3, times cos(theta) = d1/r
        lap = lap - 8.0_dp*d(1)/3.0_dp
        f = exact_solution(s) - lap
    end function source

end module spline_annulus_kernel

program test_spline_annulus
    !! Periodic outer annulus around the five-block O-grid core: the
    !! annulus is C^(p-1) across its periodic seam, the core meets it in C0
    !! with identical maps on the interface, the disc area and the solution
    !! of u - lap u = f converge at order p + 1, and the radial grading
    !! toward the wall is independent of the core.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_extraction_1d, only: SPLINE_GRID_BUMP_EDGE
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, &
        LAYOUT_ANNULUS_INTERFACE, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C1
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_evaluate_at, spline_geometry_at, &
        spline_boundary_nodes, spline_locate
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use spline_annulus_kernel, only: annulus_kernel_t, exact_solution, source, &
        center
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    real(dp) :: errors(3), areas(3), rate
    integer :: status, level

    call layout_ogrid(center, 1.0_dp, 0.3_dp, layout, status, &
        interface_fraction=0.6_dp)
    call check_condition(status == 0 .and. layout%block_count == 6 .and. &
        count(layout%side_tag == LAYOUT_ANNULUS_INTERFACE) == 5, &
        "O-grid core with a periodic annulus")
    call build(3, 4, 4, .true., space, status)
    call check_condition(status == 0, "graded annulus space builds")
    call check_condition(seam_jump(space) < 1.0e-9_dp, &
        "annulus fields are C^(p-1) across the periodic seam")
    call check_condition(interface_jump(space) < 1.0e-12_dp, &
        "core and annulus agree in values and maps on the interface")
    call check_condition(dof_count_ok(space), &
        "core interface functions are replaced by the annulus traces")
    call check_condition(wall_nodes_ok(space), &
        "the wall is one closed loop of panel nodes on the unit circle")
    call check_condition(locate_ok(space), &
        "spline_locate recovers block and parameters in core and annulus")
    do level = 1, 3
        call solve(3, 2**level, errors(level), areas(level))
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    print '(a,3es11.3,a,f6.2,a,3es11.3)', "p = 3 L2 errors", errors, " rate", rate, &
        " area errors", abs(areas - pi)
    call check_condition(rate > 3.7_dp, &
        "annulus solution converges at order p + 1 in L2")
    call check_condition(abs(areas(3) - pi) < abs(areas(2) - pi)/12.0_dp, &
        "disc area converges at order p + 1")
    call build_c1(status)
    call check_condition(status == 6, "C1 spaces reject periodic blocks")
    call check_summary("spline_annulus")

contains

    subroutine build(degree, cells, radial, graded, space, status)
        integer, intent(in) :: degree, cells, radial
        logical, intent(in) :: graded
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = degree
        allocate(options%knots(2, 6))
        options%knots(:, :)%element_count = cells
        options%knots(2, 6)%element_count = radial
        if (graded) then
            options%knots(2, 6)%grading = SPLINE_GRID_BUMP_EDGE
            options%knots(2, 6)%reversed = .true.
        end if
        call spline_space_build(layout, options, space, status)
    end subroutine build

    subroutine build_c1(status)
        integer, intent(out) :: status

        type(spline_space_options_t) :: options
        type(spline_space_t) :: c1_space

        options%degree = 3
        options%regularity = SPLINE_C1
        allocate(options%knots(2, 6))
        options%knots(:, :)%element_count = 4
        call spline_space_build(layout, options, c1_space, status)
    end subroutine build_c1

    real(dp) function seam_jump(space) result(jump)
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: coef(:, :)
        real(dp) :: left(6), right(6), t
        integer :: q, status

        allocate(coef(1, space%dof_count))
        call random_number(coef)
        jump = 0.0_dp
        do q = 0, 8
            t = real(q, dp)/8.0_dp
            call spline_evaluate_at(space, coef, 6, [0.0_dp, t], 2, left, status)
            call spline_evaluate_at(space, coef, 6, [1.0_dp, t], 2, right, status)
            jump = max(jump, maxval(abs(left - right))/maxval(abs(left)))
        end do
    end function seam_jump

    real(dp) function interface_jump(space) result(jump)
        !! Core side 1 of blocks 2-5 runs counter-clockwise over a quarter of
        !! the interface, starting at the annulus seam for block 2.
        type(spline_space_t), intent(in) :: space

        real(dp), allocatable :: coef(:, :)
        real(dp) :: core(1), ring(1), gcore(2, 6), gring(2, 6), t, xi1
        integer :: k, q, status

        allocate(coef(1, space%dof_count))
        call random_number(coef)
        jump = 0.0_dp
        do k = 2, 5
            do q = 0, 10
                t = real(q, dp)/10.0_dp
                xi1 = modulo(0.25_dp*(k - 2 + t), 1.0_dp)
                call spline_evaluate_at(space, coef, k, [t, 0.0_dp], 0, core, status)
                call spline_evaluate_at(space, coef, 6, [xi1, 1.0_dp], 0, ring, status)
                call spline_geometry_at(space, k, [t, 0.0_dp], gcore, status)
                call spline_geometry_at(space, 6, [xi1, 1.0_dp], gring, status)
                jump = max(jump, abs(core(1) - ring(1)), &
                    maxval(abs(gcore(:, 1) - gring(:, 1))))
            end do
        end do
    end function interface_jump

    logical function wall_nodes_ok(space) result(ok)
        type(spline_space_t), intent(in) :: space

        integer, allocatable :: block(:)
        real(dp), allocatable :: xi(:, :), s(:, :)
        integer :: status

        call spline_boundary_nodes(space, LAYOUT_WALL, 2, block, xi, s, status)
        ok = status == 0 .and. size(block) == 2*2*16 .and. all(block == 6)
        if (.not. ok) return
        ok = maxval(abs(norm2(s - spread(center, 2, size(block)), dim=1) - 1.0_dp)) &
            < 1.0e-3_dp .and. (s(1, 1) - center(1))*(s(2, 2) - center(2)) - &
            (s(2, 1) - center(2))*(s(1, 2) - center(1)) > 0.0_dp
    end function wall_nodes_ok

    logical function locate_ok(space) result(ok)
        type(spline_space_t), intent(in) :: space

        type(spline_quadrature_t) :: quad
        real(dp) :: xi(2), worst
        integer :: q, block, status, misses

        call spline_quadrature_build(space, 3, quad, status)
        misses = 0
        worst = 0.0_dp
        do q = 1, quad%point_count, 5
            call spline_locate(space, 1, quad%s(:, q), block, xi, status)
            if (status /= 0 .or. block /= quad%block(q)) then
                misses = misses + 1
            else
                worst = max(worst, maxval(abs(xi - quad%xi(:, q))))
            end if
        end do
        ok = misses == 0 .and. worst < 1.0e-10_dp
        if (.not. ok) print *, "locate misses", misses, " worst", worst
    end function locate_ok

    logical function dof_count_ok(space) result(ok)
        !! C0 core with n + p functions per block direction, minus its
        !! 4 (n + p) - 4 interface functions, plus the n_theta x (m + p)
        !! annulus functions, n_theta = 4 n.
        type(spline_space_t), intent(in) :: space

        integer :: n, core

        n = 4 + 3
        core = 5*(n - 2)**2 + 12*(n - 2) + 8
        ok = space%dof_count == core - (4*n - 4) + 16*(4 + 3)
    end function dof_count_ok

    subroutine solve(degree, cells, error, area)
        integer, intent(in) :: degree, cells
        real(dp), intent(out) :: error, area

        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(annulus_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :)
        real(dp) :: energy
        integer :: status, q

        call build(degree, cells, cells, .false., space, status)
        if (status /= 0) error stop "space"
        call spline_quadrature_build(space, degree + 2, quad, status)
        if (status /= 0) error stop "quadrature"
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = source(quad%s(:, q))
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        call pcg(matrix, -g(1, :), c(1, :))
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = 0.0_dp
        do q = 1, quad%point_count
            error = error + quad%weight(q)*(jet(1, q) - exact_solution(quad%s(:, q)))**2
        end do
        error = sqrt(error)
        area = sum(quad%weight)
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
            if (norm2(r) < 1.0e-13_dp*norm_b) exit
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

end program test_spline_annulus
