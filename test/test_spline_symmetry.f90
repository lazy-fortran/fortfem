module spline_symmetry_fields
    !! Parametric proxy of A = (y + z, x, x^2 + y z), A(Qx) = -Q A(x) for
    !! Q: (x, y, z) -> (x, -y, -z), on the current space and chart.
    use fortfem_kinds, only: dp
    use fortfem_spline_space, only: spline_space_t
    use fortfem_spline_eval, only: spline_geometry_at
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_map
    implicit none
    private

    type(spline_space_t), public :: field_space
    type(curve_chart_t), public :: field_chart

    public :: symmetric_a

contains

    subroutine symmetric_a(block, xi, t, values)
        integer, intent(in) :: block
        real(dp), intent(in) :: xi(2), t
        real(dp), intent(out) :: values(:)
        real(dp) :: jet(2, 6), x(3), jc(3, 3), j3(3, 3)
        integer :: status

        call spline_geometry_at(field_space, block, xi, jet, status)
        call curve_chart_map(field_chart, jet(:, 1), t, x, jc)
        j3(:, 1:2) = matmul(jc(:, 1:2), jet(:, 2:3))
        j3(:, 3) = jc(:, 3)
        values(1:3) = matmul(transpose(j3), [x(2) + x(3), x(1), x(1)**2 + x(2)*x(3)])
    end subroutine symmetric_a

end module spline_symmetry_fields

module spline_symmetry_problem
    !! u - lap u = f on a torus with circular section, in the chart
    !! coordinates: W = (|J^-T grad u|^2/2 + u^2/2 - f u) |det J|, exact
    !! u = x + y^2 + z^2 + y z, symmetric under (x, y, z) -> (x, -y, -z).
    use fortfem_kinds, only: dp
    use fortfem_spline_toroidal, only: pointwise_energy3_t
    implicit none
    private

    public :: torus_kernel_t, exact

    type, extends(pointwise_energy3_t) :: torus_kernel_t
        real(dp), allocatable :: jdet(:, :), ginv(:, :, :, :), f(:, :)
    contains
        procedure :: evaluate
    end type torus_kernel_t

contains

    subroutine evaluate(self, cell, plane, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(torus_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, plane, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status

        real(dp) :: a, gi(3, 3), h(4, 4)
        integer :: q, p

        associate (unused => cell, unused_weight => weight)
        end associate
        do q = 1, npoint
            p = first_point + q - 1
            a = self%jdet(p, plane)
            gi = self%ginv(:, :, p, plane)
            w(q) = (0.5_dp*dot_product(jet(2:4, q), matmul(gi, jet(2:4, q))) + &
                0.5_dp*jet(1, q)**2 - self%f(p, plane)*jet(1, q))*a
            if (need >= 1) gradient(:, q) = [(jet(1, q) - self%f(p, plane))*a, &
                matmul(gi, jet(2:4, q))*a]
            if (need >= 2) then
                h = 0.0_dp
                h(1, 1) = a
                h(2:4, 2:4) = gi*a
                hessian(:, q) = reshape(h, [16])
            end if
        end do
        status = 0
    end subroutine evaluate

    pure subroutine exact(x, u, grad, f)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u, grad(3), f

        u = x(1) + x(2)**2 + x(3)**2 + x(2)*x(3)
        grad = [1.0_dp, 2.0_dp*x(2) + x(3), 2.0_dp*x(3) + x(2)]
        f = u - 4.0_dp
    end subroutine exact

end module spline_symmetry_problem

program test_spline_symmetry
    !! Stellarator symmetry with Fourier modes in t: the mirror map of an
    !! up-down symmetric section, rejection of an asymmetric one, expanded
    !! reduced fields with u(s1, -s2, -t) = parity u(s1, s2, t), and a
    !! symmetric problem whose reduced solve (P^T H P) equals the full solve
    !! with about half the unknowns.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_side_quadrature_t, spline_side_quadrature, spline_point_matrix, &
        spline_locate
    use fortfem_spline_toroidal, only: spline_toroidal_t, spline_toroidal_fourier, &
        spline_energy3, spline_energy3_hessian_csr, spline_evaluate3_at, &
        spline_symmetry_t, spline_toroidal_symmetry, spline_symmetry_expand, &
        spline_symmetry_reduce, spline_symmetry_reduce_matrix
    use fortfem_curve_chart, only: curve_chart_t, curve_chart_cylindrical, &
        curve_chart_map
    use fortsparse, only: pcg_params_t, pcg_handle_t, fortsparse_status_t, &
        pcg_setup_csr, pcg_solve, PCG_IC0
    use spline_symmetry_problem, only: torus_kernel_t, exact
    use spline_symmetry_fields, only: field_space, field_chart, symmetric_a
    use fortfem_spline_derham, only: spline_derham_t, spline_derham_build, &
        spline_derham_symmetry, spline_derham_project
    implicit none

    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    type(spline_toroidal_t) :: tor
    type(spline_symmetry_t) :: sym
    type(curve_chart_t) :: chart
    integer :: status

    call curve_chart_cylindrical(chart, 1)
    call spline_toroidal_fourier(2, 12, tor, status)
    call layout_ogrid([2.0_dp, 0.1_dp], 0.5_dp, 0.5_dp, layout, status)
    call build(4)
    call spline_toroidal_symmetry(space, tor, [1], sym, status)
    call check_condition(status == 2, "a section off the midplane is rejected")
    call layout_ogrid([2.0_dp, 0.0_dp], 0.5_dp, 0.5_dp, layout, status)
    call build(4)
    call spline_toroidal_symmetry(space, tor, [1, -1], sym, status)
    call check_condition(status == 0 .and. 2*sym%reduced_count < &
        sym%full_count + 2*tor%mode_count*2*count(sym%mirror == &
        [(status, status = 1, space%dof_count)]), "mirror map and reduced basis")
    call check_condition(parity_ok(), &
        "expanded fields satisfy u(s1, -s2, -t) = parity u(s1, s2, t)")
    call check_condition(solve_ok(), &
        "the reduced symmetric solve equals the full solve")
    call check_condition(derham_ok(), &
        "de Rham classes: d maps class into class, projected symmetric A is odd")
    call check_summary("spline_symmetry")

contains

    subroutine build(cells)
        integer, intent(in) :: cells
        type(spline_space_options_t) :: options
        integer :: status

        options%degree = 3
        allocate(options%knots(2, layout%block_count))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
    end subroutine build

    logical function parity_ok() result(ok)
        real(dp), allocatable :: red(:), full(:)
        real(dp) :: s(2), xi(2), xm(2), t, a(2), b(2), worst
        integer :: k, block, mblock, status, q

        allocate(red(sym%reduced_count), full(sym%full_count))
        call random_number(red)
        call spline_symmetry_expand(sym, red, full)
        worst = 0.0_dp
        do q = 1, 40
            call random_number(s)
            call random_number(t)
            s = [2.0_dp, 0.0_dp] + 0.3_dp*(s - 0.5_dp)
            call spline_locate(space, 1, s, block, xi, status)
            call spline_locate(space, 1, [s(1), -s(2)], mblock, xm, status)
            call spline_evaluate3_at(space, tor, reshape(full, [2, tor%mode_count, &
                space%dof_count]), block, xi, t, 0, a, status)
            call spline_evaluate3_at(space, tor, reshape(full, [2, tor%mode_count, &
                space%dof_count]), mblock, xm, -t, 0, b, status)
            worst = max(worst, abs(b(1) - a(1)), abs(b(2) + a(2)))
        end do
        ok = worst < 1.0e-12_dp*maxval(abs(full))
        associate (unused => k)
        end associate
    end function parity_ok

    logical function solve_ok() result(ok)
        type(spline_symmetry_t) :: one
        type(spline_quadrature_t) :: quad
        type(spline_side_quadrature_t) :: side
        type(torus_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix, small
        real(dp), allocatable :: c(:, :, :), g(:, :, :), b(:), val(:), x(:), xr(:)
        real(dp), allocatable :: br(:), back(:)
        integer, allocatable :: rowptr(:), colidx(:)
        real(dp) :: energy, xx(3), jac(3, 3), xl(3), xn(3), nrm(3), u, gr(3), f, det
        integer :: status, q, j, k

        call spline_toroidal_symmetry(space, tor, [1], one, status)
        call spline_quadrature_build(space, 4, quad, status)
        allocate(kernel%jdet(quad%point_count, tor%plane_count), &
            kernel%ginv(3, 3, quad%point_count, tor%plane_count), &
            kernel%f(quad%point_count, tor%plane_count))
        do j = 1, tor%plane_count
            do q = 1, quad%point_count
                call curve_chart_map(chart, quad%s(:, q), tor%t(j), xx, jac)
                det = jac(1, 1)*(jac(2, 2)*jac(3, 3) - jac(3, 2)*jac(2, 3)) - &
                    jac(1, 2)*(jac(2, 1)*jac(3, 3) - jac(3, 1)*jac(2, 3)) + &
                    jac(1, 3)*(jac(2, 1)*jac(3, 2) - jac(3, 1)*jac(2, 2))
                kernel%jdet(q, j) = abs(det)
                kernel%ginv(:, :, q, j) = inverse3(matmul(transpose(jac), jac))
                call exact(xx, u, gr, f)
                kernel%f(q, j) = f
            end do
        end do
        allocate(c(1, tor%mode_count, space%dof_count), &
            g(1, tor%mode_count, space%dof_count))
        c = 0.0_dp
        call spline_energy3(space, quad, tor, kernel, c, energy, status, g)
        call spline_energy3_hessian_csr(space, quad, tor, kernel, c, matrix, status)
        b = -reshape(g, [size(g)])
        call spline_side_quadrature(space, LAYOUT_WALL, 5, side, status)
        call spline_point_matrix(space, side%block, side%xi, 0, rowptr, colidx, val, &
            status)
        do q = 1, side%point_count
            do j = 1, tor%plane_count
                call curve_chart_map(chart, side%s(:, q), tor%t(j), xx, jac)
                xl = jac(:, 1)*side%tangent(1, q) + jac(:, 2)*side%tangent(2, q)
                xn = jac(:, 1)*side%tangent(2, q) - jac(:, 2)*side%tangent(1, q)
                nrm = [xl(2)*jac(3, 3) - xl(3)*jac(2, 3), xl(3)*jac(1, 3) - &
                    xl(1)*jac(3, 3), xl(1)*jac(2, 3) - xl(2)*jac(1, 3)]
                if (dot_product(nrm, xn) < 0.0_dp) nrm = -nrm
                call exact(xx, u, gr, f)
                do k = rowptr(q), rowptr(q + 1) - 1
                    b(1 + tor%mode_count*(colidx(k) - 1):tor%mode_count*colidx(k)) = &
                        b(1 + tor%mode_count*(colidx(k) - 1):tor%mode_count*colidx(k)) &
                        + side%weight(q)*tor%weight(j)*dot_product(gr, nrm)* &
                        tor%table(:, j, 0)*val(k)
                end do
            end do
        end do
        x = solve(matrix, b)
        call spline_symmetry_reduce_matrix(one, matrix, small)
        allocate(br(one%reduced_count), back(one%full_count))
        call spline_symmetry_reduce(one, b, br)
        xr = solve(small, br)
        call spline_symmetry_expand(one, xr, back)
        print '(a,i0,a,i0,a,es10.2)', "full ", size(x), ", reduced ", size(xr), &
            ", max difference ", maxval(abs(back - x))
        ok = maxval(abs(back - x)) < 1.0e-8_dp*maxval(abs(x)) .and. &
            2*size(xr) < size(x) + size(x)/5
    end function solve_ok

    logical function derham_ok() result(ok)
        type(spline_derham_t) :: derham
        type(spline_symmetry_t) :: s(0:3, -1:1)
        real(dp), allocatable :: c(:), d(:)
        integer :: form, par, status

        call spline_derham_build(space, derham, status, 2)
        ok = status == 0
        do form = 0, 3
            do par = -1, 1, 2
                call spline_derham_symmetry(derham, form, par, s(form, par), status)
                ok = ok .and. status == 0
            end do
        end do
        if (.not. ok) return
        do par = -1, 1, 2
            c = random_member(s(0, par))
            ok = ok .and. in_class(s(1, par), apply(derham%grad, c))
            c = random_member(s(1, par))
            ok = ok .and. in_class(s(2, par), apply(derham%curl, c))
            c = random_member(s(2, par))
            ok = ok .and. in_class(s(3, par), apply(derham%div, c))
        end do
        field_space = space
        field_chart = chart
        allocate(d(derham%dimension(1)))
        call spline_derham_project(derham, 1, symmetric_a, d, status, 10)
        ok = ok .and. status == 0 .and. in_class(s(1, -1), d) .and. &
            .not. in_class(s(1, 1), d)
    end function derham_ok

    function random_member(sym) result(full)
        type(spline_symmetry_t), intent(in) :: sym
        real(dp), allocatable :: full(:), red(:)

        allocate(red(sym%reduced_count), full(sym%full_count))
        call random_number(red)
        call spline_symmetry_expand(sym, red, full)
    end function random_member

    logical function in_class(sym, full) result(ok)
        !! full equals its orthogonal projection P P^T full onto the class.
        type(spline_symmetry_t), intent(in) :: sym
        real(dp), intent(in) :: full(:)
        real(dp), allocatable :: red(:), back(:)

        allocate(red(sym%reduced_count), back(sym%full_count))
        call spline_symmetry_reduce(sym, full, red)
        call spline_symmetry_expand(sym, red, back)
        ok = maxval(abs(back - full)) < 1.0e-10_dp*max(1.0_dp, maxval(abs(full)))
    end function in_class

    function apply(a, x) result(y)
        type(csr_matrix_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: i, k

        allocate(y(a%n))
        y = 0.0_dp
        do i = 1, a%n
            do k = a%rowptr(i), a%rowptr(i + 1) - 1
                y(i) = y(i) + a%val(k)*x(a%colidx(k))
            end do
        end do
    end function apply

    function solve(a, b) result(x)
        type(csr_matrix_t), intent(in) :: a
        real(dp), intent(in) :: b(:)
        real(dp), allocatable :: x(:)
        type(pcg_params_t) :: params
        type(pcg_handle_t) :: handle
        type(fortsparse_status_t) :: info
        real(dp) :: residual
        integer :: iterations

        params%preconditioner = PCG_IC0
        params%rtol = 1.0e-13_dp
        params%max_iterations = 10*a%n
        call pcg_setup_csr(a%n, a%rowptr, a%colidx, a%val, params, handle, info)
        allocate(x(a%n))
        x = 0.0_dp
        call pcg_solve(handle, params, b, x, info, iterations, residual)
    end function solve

    pure function inverse3(a) result(b)
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: b(3, 3), det

        b(1, 1) = a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)
        b(1, 2) = a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3)
        b(1, 3) = a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2)
        b(2, 1) = a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3)
        b(2, 2) = a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1)
        b(2, 3) = a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3)
        b(3, 1) = a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1)
        b(3, 2) = a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2)
        b(3, 3) = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
        det = a(1, 1)*b(1, 1) + a(1, 2)*b(2, 1) + a(1, 3)*b(3, 1)
        b = b/det
    end function inverse3

end program test_spline_symmetry
