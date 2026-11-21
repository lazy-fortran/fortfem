module spline_separatrix_problem
    !! Doublet bounded by the lemniscate F = (x^2 + z^2)^2 - a^2 (z^2 - x^2)
    !! = 0, x = R - r0, with the X-point (r0, 0), lobes F < 0 above and below,
    !! and an elliptic wall. Near X, F = -a^2 (z - x)(z + x), so log|F| has
    !! the logarithmic singularity of a Clebsch potential at a separatrix.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    use fortfem_spline_separatrix, only: doublet_curves_t
    implicit none
    private

    public :: projection_kernel_t, doublet, flux, smooth

    real(dp), parameter, public :: r0 = 1.5_dp
    real(dp), parameter, public :: wall_a = 0.8_dp, wall_b = 1.1_dp

    type, extends(pointwise_energy_t) :: projection_kernel_t
        !! w = u^2/2 - f u: the Hessian is the mass matrix.
        real(dp), allocatable :: f(:)
    contains
        procedure :: evaluate
    end type projection_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(projection_kernel_t), intent(in) :: self
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
            w(q) = 0.5_dp*jet(1, q)**2 - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [jet(1, q) - f, 0.0_dp, 0.0_dp]
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        end do
        status = 0
    end subroutine evaluate

    pure real(dp) function flux(s, a) result(f)
        real(dp), intent(in) :: s(2), a

        f = ((s(1) - r0)**2 + s(2)**2)**2 - a**2*(s(2)**2 - (s(1) - r0)**2)
    end function flux

    pure real(dp) function smooth(s) result(u)
        real(dp), intent(in) :: s(2)

        u = cos(2.0_dp*s(1))*sin(1.7_dp*s(2)) + s(1)*s(2)
    end function smooth

    function doublet(a) result(curves)
        !! Lobe 1 above, lobe 2 below, from the lemniscate parametrisation
        !! x = a sin t cos t/(1 + sin^2 t), z = a cos t/(1 + sin^2 t).
        real(dp), intent(in) :: a
        type(doublet_curves_t) :: curves

        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: t, theta
        integer :: i, k
        integer, parameter :: n = 401

        curves%xpoint = [r0, 0.0_dp]
        curves%opoint(:, 1) = [r0, a/sqrt(2.0_dp)]
        curves%opoint(:, 2) = [r0, -a/sqrt(2.0_dp)]
        do k = 1, 2
            allocate(curves%lobe(k)%points(2, n))
            do i = 1, n
                theta = real(i - 1, dp)/(n - 1)
                if (k == 1) then
                    t = 0.5_dp*pi - pi*theta
                else
                    t = 0.5_dp*pi + pi*theta
                end if
                curves%lobe(k)%points(:, i) = [r0, 0.0_dp] + &
                    a*[sin(t)*cos(t), cos(t)]/(1.0_dp + sin(t)**2)
            end do
        end do
        allocate(curves%wall(2, 400))
        do i = 1, 400
            theta = 2.0_dp*pi*(i - 1)/400
            curves%wall(:, i) = [r0 + wall_a*cos(theta), wall_b*sin(theta)]
        end do
    end function doublet

end module spline_separatrix_problem

program test_spline_separatrix
    !! Doublet layout (doc/design/multipatch_spline_spaces.md, separatrix
    !! layouts): 18 blocks with X of valence 8, order p + 1 for a smooth
    !! L2 projection, a cut space that resolves a jump across lobe 1 at
    !! order p + 1 with one-sided traces at the interface, algebraic grading
    !! for log|F|, and a moved doublet on the same space.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, LAYOUT_SEPARATRIX
    use fortfem_spline_separatrix, only: layout_doublet
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_space_update_geometry, spline_grade_knots
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate, spline_side_quadrature_t, spline_interface_quadrature, &
        spline_point_matrix
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr
    use spline_separatrix_problem, only: projection_kernel_t, doublet, flux, &
        smooth, r0, wall_a, wall_b
    implicit none

    interface
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface

    real(dp), parameter :: pi = acos(-1.0_dp), a = 0.6_dp
    integer, parameter :: smooth_case = 1, jump_case = 2, log_case = 3
    type(block_layout_t) :: layout, moved
    type(spline_space_t) :: space
    real(dp) :: errors(3), cut_errors(3), whole(3), area(3), rate, cut_rate
    real(dp) :: uniform(2), graded(2), trace_gap, jump_error, error
    integer :: status, level, dofs, x_vertex

    call layout_doublet(doublet(a), 0.4_dp, layout, status)
    x_vertex = minloc(norm2(layout%vertices - spread([r0, 0.0_dp], 2, &
        layout%vertex_count), dim=1), 1)
    call check_condition(status == 0 .and. layout%block_count == 18 .and. &
        layout%valence(x_vertex) == 8 .and. &
        count(layout%side_tag == LAYOUT_SEPARATRIX) == 16, &
        "doublet layout: 18 blocks, X of valence 8, 8 separatrix pieces")

    do level = 1, 3
        call build(2**(level + 1), .false., .false., space, status)
        if (status /= 0) error stop "space"
        call project(space, smooth_case, errors(level), area(level))
        call project(space, jump_case, whole(level), area(level))
        call build(2**(level + 1), .true., .false., space, status)
        if (status /= 0) error stop "cut space"
        call project(space, jump_case, cut_errors(level), area(level), trace_gap, &
            jump_error)
    end do
    rate = log(errors(2)/errors(3))/log(2.0_dp)
    cut_rate = log(cut_errors(2)/cut_errors(3))/log(2.0_dp)
    print '(a,3es11.3,a,f6.2)', "smooth L2 errors", errors, " rate", rate
    print '(a,3es11.3,a,f6.2)', "jump, cut space ", cut_errors, " rate", cut_rate
    print '(a,3es11.3)', "jump, C0 space  ", whole
    print '(a,es11.3)', "area error", abs(area(3) - pi*wall_a*wall_b)
    call check_condition(rate > 3.5_dp .and. &
        abs(area(3) - pi*wall_a*wall_b) < 1.0e-5_dp, &
        "smooth projection converges at order p + 1 on the doublet")
    call check_condition(cut_rate > 3.5_dp .and. whole(3) > 1.0e3_dp*cut_errors(3), &
        "a jump across the separatrix: order p + 1 in the cut space only")
    print '(a,2es11.3)', "trace gap of the C0 field, jump error", trace_gap, &
        jump_error
    call check_condition(trace_gap < 1.0e-12_dp .and. jump_error < 1.0e-4_dp, &
        "one-sided traces: equal for a C0 field, the exact jump in the cut space")

    do level = 1, 2
        call build(8*(level + 1), .true., .false., space, status)
        call project(space, log_case, uniform(level), area(1))
        call build(8*(level + 1), .true., .true., space, status)
        if (status /= 0) error stop "graded space"
        call project(space, log_case, graded(level), area(1))
    end do
    rate = log(uniform(1)/uniform(2))/log(1.5_dp)
    cut_rate = log(graded(1)/graded(2))/log(1.5_dp)
    print '(a,2es11.3,a,f6.2,a,2es11.3,a,f6.2)', "log|F| uniform", uniform, &
        " rate", rate, " graded", graded, " rate", cut_rate
    call check_condition(rate < 0.7_dp .and. cut_rate > 1.5_dp .and. &
        graded(2) < 0.1_dp*uniform(2), &
        "power grading towards the separatrix and X resolves log|F|")
    call build(16, .true., .true., space, status, ratio=1.6_dp)
    call project(space, log_case, error, area(1))
    print '(a,es11.3)', "log|F| geometric ratio 1.6, 16 cells", error
    call check_condition(error < 0.25_dp*uniform(1), &
        "geometric grading lowers the log|F| error")

    call build(16, .true., .false., space, status)
    dofs = space%dof_count
    call layout_doublet(doublet(0.55_dp), 0.4_dp, moved, status)
    call spline_space_update_geometry(space, moved, status)
    call check_condition(status == 0 .and. space%dof_count == dofs, &
        "a moved doublet keeps the space")
    call project(space, smooth_case, error, area(1))
    call check_condition(error < 2.0_dp*errors(3), "the moved doublet projects well")
    call check_summary("spline_separatrix")

contains

    subroutine build(cells, cut, grade, space, status, ratio)
        !! grade: power 4, or geometric with ratio.
        integer, intent(in) :: cells
        logical, intent(in) :: cut, grade
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status
        real(dp), intent(in), optional :: ratio

        type(spline_space_options_t) :: options

        options%degree = 3
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        if (cut) options%cut_tags = [LAYOUT_SEPARATRIX]
        if (grade) then
            if (present(ratio)) then
                call spline_grade_knots(layout, LAYOUT_SEPARATRIX, ratio, &
                    options%knots, status)
            else
                call spline_grade_knots(layout, LAYOUT_SEPARATRIX, 1.0_dp, &
                    options%knots, status, power=4.0_dp)
            end if
            if (status /= 0) return
        end if
        call spline_space_build(layout, options, space, status)
    end subroutine build

    real(dp) function exact(case, s, block) result(u)
        integer, intent(in) :: case, block
        real(dp), intent(in) :: s(2)

        select case (case)
        case (smooth_case)
            u = smooth(s)
        case (jump_case)
            u = smooth(s) + merge(1.0_dp, 0.0_dp, block <= 5)
        case default
            u = log(abs(flux(s, a)))
        end select
    end function exact

    subroutine project(space, case, error, area, trace_gap, jump_error)
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: case
        real(dp), intent(out) :: error, area
        real(dp), intent(out), optional :: trace_gap, jump_error

        type(spline_quadrature_t) :: quad
        type(spline_side_quadrature_t) :: side
        type(projection_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :), b(:), minus(:), plus(:)
        real(dp), allocatable :: c0(:)
        real(dp) :: energy, expected
        integer :: status, q, k

        call spline_quadrature_build(space, 6, quad, status)
        if (status /= 0) error stop "quadrature"
        allocate(kernel%f(quad%point_count))
        do q = 1, quad%point_count
            kernel%f(q) = exact(case, quad%s(:, q), quad%block(q))
        end do
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        b = -g(1, :)
        call solve(matrix, b, c(1, :))
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = 0.0_dp
        do q = 1, quad%point_count
            error = error + quad%weight(q)*(jet(1, q) - kernel%f(q))**2
        end do
        error = sqrt(error)
        area = sum(quad%weight)
        if (.not. present(trace_gap)) return
        call spline_interface_quadrature(space, LAYOUT_SEPARATRIX, 4, side, status)
        if (status /= 0 .or. side%point_count /= &
            8*4*space%blocks(1)%element_count(1)) error stop "interface"
        minus = traces(space, side%block, side%xi, c(1, :))
        plus = traces(space, side%neighbour_block, side%neighbour_xi, c(1, :))
        jump_error = 0.0_dp
        do q = 1, side%point_count
            expected = merge(1.0_dp, 0.0_dp, side%block(q) <= 5) - &
                merge(1.0_dp, 0.0_dp, side%neighbour_block(q) <= 5)
            jump_error = max(jump_error, abs(minus(q) - plus(q) - expected))
        end do
        block_c0: block
            type(spline_space_t) :: whole
            call build(16, .false., .false., whole, status)
            c0 = [(sin(0.37_dp*k), k = 1, whole%dof_count)]
            trace_gap = maxval(abs(traces(whole, side%block, side%xi, c0) - &
                traces(whole, side%neighbour_block, side%neighbour_xi, c0)))
        end block block_c0
    end subroutine project

    function traces(space, block, xi, c) result(u)
        type(spline_space_t), intent(in) :: space
        integer, intent(in) :: block(:)
        real(dp), intent(in) :: xi(:, :), c(:)
        real(dp) :: u(size(block))

        real(dp), allocatable :: val(:)
        integer, allocatable :: rowptr(:), colidx(:)
        integer :: q, status

        call spline_point_matrix(space, block, xi, 0, rowptr, colidx, val, status)
        if (status /= 0) error stop "point matrix"
        do q = 1, size(block)
            u(q) = sum(val(rowptr(q):rowptr(q + 1) - 1)* &
                c(colidx(rowptr(q):rowptr(q + 1) - 1)))
        end do
    end function traces

    subroutine solve(matrix, b, x)
        !! Dense Cholesky of the mass matrix.
        type(csr_matrix_t), intent(in) :: matrix
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)

        real(dp), allocatable :: a(:, :), rhs(:, :)
        integer :: row, k, info

        allocate(a(matrix%n, matrix%n), rhs(matrix%n, 1))
        a = 0.0_dp
        do row = 1, matrix%n
            do k = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                a(row, matrix%colidx(k)) = matrix%val(k)
            end do
        end do
        rhs(:, 1) = b
        call dposv("L", matrix%n, 1, a, matrix%n, rhs, matrix%n, info)
        if (info /= 0) error stop "dposv"
        x = rhs(:, 1)
    end subroutine solve

end program test_spline_separatrix
