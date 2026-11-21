module spline_poisson_kernel
    !! W = |grad u|^2/2 - f u on the O-grid disc, with the exact solution
    !! u = (1 - r^2) cos(s1 + s2/2), so the minimiser solves -lap u = f.
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t, SPLINE_HESSIAN_PACKED
    implicit none
    private

    public :: poisson_kernel_t, exact_solution

    type, extends(pointwise_energy_t) :: poisson_kernel_t
        real(dp), allocatable :: s(:, :)
    contains
        procedure :: evaluate
    end type poisson_kernel_t

contains

    subroutine evaluate(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(poisson_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: f
        integer :: q

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            f = source(self%s(:, first_point + q - 1))
            w(q) = 0.5_dp*(jet(2, q)**2 + jet(3, q)**2) - f*jet(1, q)
            if (need >= 1) gradient(:, q) = [-f, jet(2, q), jet(3, q)]
            if (need >= 2) then
                if (self%hessian_storage == SPLINE_HESSIAN_PACKED) then
                    hessian(:, q) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
                else
                    hessian(:, q) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
                end if
            end if
        end do
        status = 0
    end subroutine evaluate

    pure real(dp) function exact_solution(s) result(u)
        real(dp), intent(in) :: s(2)

        u = (1.0_dp - s(1)**2 - s(2)**2)*cos(s(1) + 0.5_dp*s(2))
    end function exact_solution

    pure real(dp) function source(s) result(f)
        !! f = -lap u = 4 b - 4 sin(phi) (s1 + s2/2) + 1.25 a b with
        !! a = 1 - r^2, b = cos(phi), phi = s1 + s2/2.
        real(dp), intent(in) :: s(2)

        real(dp) :: a, b, phi

        phi = s(1) + 0.5_dp*s(2)
        a = 1.0_dp - s(1)**2 - s(2)**2
        b = cos(phi)
        f = 4.0_dp*b - 4.0_dp*sin(phi)*(s(1) + 0.5_dp*s(2)) + 1.25_dp*a*b
    end function source

end module spline_poisson_kernel

program test_spline_energy_poisson
    !! Poisson on the C0 O-grid disc by minimising a pointwise energy:
    !! the assembled and matrix-free Hessians agree, the gradient matches
    !! finite differences of the energy, and the L2 error converges at
    !! order p + 1 for p = 3, 4, 5.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, spline_boundary_dofs
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_evaluate
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr, &
        spline_energy_hessian_apply, SPLINE_HESSIAN_DENSE, SPLINE_HESSIAN_PACKED
    use spline_poisson_kernel, only: poisson_kernel_t, exact_solution
    implicit none

    type(block_layout_t) :: layout
    real(dp) :: errors(3), rate
    integer :: status, degree, level

    call layout_ogrid([0.0_dp, 0.0_dp], 1.0_dp, 0.5_dp, layout, status)
    call check_condition(consistency(), &
        "CSR and matrix-free Hessians agree and the gradient matches differences")
    do degree = 3, 5
        do level = 1, 3
            errors(level) = solve_error(degree, 2**level)
        end do
        rate = log(errors(2)/errors(3))/log(2.0_dp)
        print '(a,i0,a,3es11.3,a,f6.2)', "p = ", degree, " L2 errors", errors, &
            " rate", rate
        call check_condition(rate > real(degree, dp) + 0.7_dp, &
            "O-grid Poisson converges at order p + 1 in L2")
    end do
    call check_summary("spline_energy_poisson")

contains

    subroutine setup(degree, cells, space, quad, kernel)
        integer, intent(in) :: degree, cells
        type(spline_space_t), intent(out) :: space
        type(spline_quadrature_t), intent(out) :: quad
        type(poisson_kernel_t), intent(out) :: kernel

        type(spline_space_options_t) :: options
        integer :: status

        options%degree = degree
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
        call spline_quadrature_build(space, degree + 1, quad, status)
        if (status /= 0) error stop "quadrature"
        kernel%s = quad%s
    end subroutine setup

    logical function consistency() result(ok)
        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(poisson_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), v(:, :), hv(:, :), hv_packed(:, :)
        real(dp), allocatable :: g(:, :), cv(:)
        real(dp) :: e_plus, e_minus, e0, delta, slope
        integer :: status, i, k

        call setup(3, 2, space, quad, kernel)
        allocate(c(1, space%dof_count), v(1, space%dof_count), &
            hv(1, space%dof_count), hv_packed(1, space%dof_count), &
            g(1, space%dof_count), cv(space%dof_count))
        call random_number(c)
        call random_number(v)
        kernel%hessian_storage = SPLINE_HESSIAN_DENSE
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        ok = status == 0
        call spline_energy_hessian_apply(space, quad, kernel, c, v, hv, status)
        ok = ok .and. status == 0
        kernel%hessian_storage = SPLINE_HESSIAN_PACKED
        call spline_energy_hessian_apply(space, quad, kernel, c, v, hv_packed, status)
        ok = ok .and. status == 0
        do i = 1, matrix%n
            cv(i) = 0.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                cv(i) = cv(i) + matrix%val(k)*v(1, matrix%colidx(k))
            end do
        end do
        ok = ok .and. maxval(abs(cv - hv(1, :))) < 1.0e-11_dp*maxval(abs(cv))
        ok = ok .and. maxval(abs(hv_packed - hv)) < 1.0e-12_dp*maxval(abs(cv))
        call spline_energy(space, quad, kernel, c, e0, status, g)
        delta = 1.0e-6_dp
        call spline_energy(space, quad, kernel, c + delta*v, e_plus, status)
        call spline_energy(space, quad, kernel, c - delta*v, e_minus, status)
        slope = (e_plus - e_minus)/(2*delta)
        ok = ok .and. abs(slope - sum(g*v)) < 1.0e-7_dp*abs(slope)
    end function consistency

    real(dp) function solve_error(degree, cells) result(error)
        !! Minimiser with u = 0 on the wall DOFs, by Jacobi-preconditioned CG
        !! on the assembled Hessian.
        integer, intent(in) :: degree, cells

        type(spline_space_t) :: space
        type(spline_quadrature_t) :: quad
        type(poisson_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), g(:, :), jet(:, :)
        logical, allocatable :: free(:)
        integer, allocatable :: wall(:)
        real(dp) :: energy
        integer :: status, q

        call setup(degree, cells, space, quad, kernel)
        allocate(c(1, space%dof_count), g(1, space%dof_count))
        c = 0.0_dp
        call spline_energy(space, quad, kernel, c, energy, status, g)
        call spline_energy_hessian_csr(space, quad, kernel, c, matrix, status)
        call spline_boundary_dofs(space, LAYOUT_WALL, wall)
        allocate(free(space%dof_count))
        free = .true.
        free(wall) = .false.
        call pcg(matrix, free, -g(1, :), c(1, :))
        allocate(jet(1, quad%point_count))
        call spline_evaluate(space, quad, c, 0, jet, status)
        error = 0.0_dp
        do q = 1, quad%point_count
            error = error + quad%weight(q)*(jet(1, q) - exact_solution(quad%s(:, q)))**2
        end do
        error = sqrt(error)
    end function solve_error

    subroutine pcg(matrix, free, b, x)
        type(csr_matrix_t), intent(in) :: matrix
        logical, intent(in) :: free(:)
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
        r = merge(b, 0.0_dp, free)
        norm_b = norm2(r)
        z = r/diag
        p = z
        rz = dot_product(r, z)
        do it = 1, 20*size(b)
            call matvec(matrix, free, p, ap)
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

    subroutine matvec(matrix, free, v, y)
        type(csr_matrix_t), intent(in) :: matrix
        logical, intent(in) :: free(:)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer :: row, k

        do row = 1, matrix%n
            y(row) = 0.0_dp
            if (.not. free(row)) cycle
            do k = matrix%rowptr(row), matrix%rowptr(row + 1) - 1
                if (free(matrix%colidx(k))) &
                    y(row) = y(row) + matrix%val(k)*v(matrix%colidx(k))
            end do
        end do
    end subroutine matvec

end program test_spline_energy_poisson
