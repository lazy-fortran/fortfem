module spline_singular_kernels
    use fortfem_kinds, only: dp
    use fortfem_spline_energy, only: pointwise_energy_t
    implicit none
    private

    public :: plain_kernel_t, weighted_kernel_t, exact, source

    type, extends(pointwise_energy_t) :: plain_kernel_t
    contains
        procedure :: evaluate => plain
    end type plain_kernel_t

    !> W = ((|grad u|^2 + u^2)/2 - f u)/rho, rho = |s - X|
    type, extends(pointwise_energy_t) :: weighted_kernel_t
        real(dp), allocatable :: inv_rho(:), f(:)
    contains
        procedure :: evaluate => weighted
    end type weighted_kernel_t

contains

    subroutine plain(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(plain_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        integer :: q

        associate (unused => cell, u2 => weight, u3 => first_point, u4 => self)
        end associate
        do q = 1, npoint
            w(q) = 0.5_dp*sum(jet(1:3, q)**2)
            if (need >= 1) gradient(:, q) = jet(1:3, q)
            if (need >= 2) hessian(:, q) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
        end do
        status = 0
    end subroutine plain

    subroutine weighted(self, cell, first_point, npoint, jet, weight, need, w, &
            gradient, hessian, status)
        class(weighted_kernel_t), intent(in) :: self
        integer, intent(in) :: cell, first_point, npoint, need
        real(dp), intent(in) :: jet(:, :), weight(:)
        real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
        integer, intent(out) :: status
        real(dp) :: a, f
        integer :: q

        associate (unused => cell, u2 => weight)
        end associate
        do q = 1, npoint
            a = self%inv_rho(first_point + q - 1)
            f = self%f(first_point + q - 1)
            w(q) = (0.5_dp*sum(jet(1:3, q)**2) - f*jet(1, q))*a
            if (need >= 1) gradient(:, q) = [(jet(1, q) - f)*a, jet(2:3, q)*a]
            if (need >= 2) hessian(:, q) = [a, 0.0_dp, 0.0_dp, 0.0_dp, a, 0.0_dp, &
                0.0_dp, 0.0_dp, a]
        end do
        status = 0
    end subroutine weighted

    pure subroutine exact(s, x0, u, grad)
        real(dp), intent(in) :: s(2), x0(2)
        real(dp), intent(out) :: u, grad(2)
        real(dp) :: r2

        r2 = sum((s - x0)**2)
        u = cos(r2)
        grad = -2.0_dp*sin(r2)*(s - x0)
    end subroutine exact

    pure real(dp) function source(s, x0) result(f)
        real(dp), intent(in) :: s(2), x0(2)
        real(dp) :: r2

        r2 = sum((s - x0)**2)
        f = cos(r2)*(1.0_dp + 4.0_dp*r2) + 2.0_dp*sin(r2)
    end function source

end module spline_singular_kernels

program test_spline_singular
    !! Polar (Duffy) quadrature about a point X in the cells around it: X
    !! inside a block, on a block side (two cells) and on a valence-3
    !! corner (three cells). The adapted rule keeps the area exact and beats
    !! tensor Gauss on the exact integral of 1/|s - X| over the square (by
    !! 1e2 to 1e3 inside a block, 2 to 5 on sides and corners, where the
    !! nearly singular neighbour cells dominate), and the pointwise cells
    !! reproduce tensor results for polynomial integrands.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_quadrature_point_singular, spline_evaluate, spline_evaluate_at, &
        spline_geometry_at
    use fortfem_spline_energy, only: spline_energy, spline_energy_hessian_csr, &
        spline_energy_hessian_apply
    use spline_singular_kernels, only: plain_kernel_t
    implicit none

    real(dp), parameter :: c0(2) = [0.3_dp, -0.2_dp], a = 1.0_dp/sqrt(2.0_dp)
    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    type(spline_quadrature_t) :: quad, adapted
    real(dp) :: points(2, 3), jet(2, 6), gauss(3), polar(3)
    integer :: status, level, k, cells(3)
    character(len=*), parameter :: names(3) = [character(len=12) :: "interior", &
        "block side", "corner"]
    logical :: ok

    call layout_ogrid(c0, 1.0_dp, 0.5_dp, layout, status, polygon=.true.)
    call build(4)
    points(:, 1) = c0 + [0.137_dp, -0.071_dp]
    call spline_geometry_at(space, 1, [0.37_dp, 0.0_dp], jet, status)
    points(:, 2) = jet(:, 1)
    call spline_geometry_at(space, 1, [0.0_dp, 0.0_dp], jet, status)
    points(:, 3) = jet(:, 1)
    do k = 1, 3
        ok = .true.
        do level = 1, 3
            call build(2**(level + 1))
            call spline_quadrature_build(space, 8, quad, status)
            call spline_quadrature_point_singular(space, quad, points(:, k), adapted, &
                status)
            ok = ok .and. status == 0 .and. abs(sum(adapted%weight) - 2.0_dp) < 1.0e-12_dp
            cells(level) = count(adapted%pointwise)
            gauss(level) = abs(inverse_distance(quad, points(:, k)) - &
                exact_integral(points(:, k)))
            polar(level) = abs(inverse_distance(adapted, points(:, k)) - &
                exact_integral(points(:, k)))
        end do
        print '(a12,a,3es10.2,a,3es10.2,a,3i3)', names(k), " Gauss", gauss, &
            "  polar", polar, "  cells", cells
        call check_condition(ok .and. all(cells == k) .and. all(polar < gauss), &
            "polar cells around X ("//trim(names(k))//"): exact area, smaller error")
    end do

    call build(4)
    call spline_quadrature_build(space, 16, quad, status)
    call spline_quadrature_point_singular(space, quad, points(:, 1), adapted, status)
    call check_condition(status == 0 .and. pointwise_ok(), &
        "pointwise cells reproduce energy, Hessians and evaluation")
    call spline_quadrature_build(space, 3, quad, status)
    call spline_quadrature_point_singular(space, quad, points(:, 1), adapted, status)
    call check_condition(status == 1, "odd point counts are rejected")
    call check_summary("spline_singular")

contains

    subroutine build(cells)
        integer, intent(in) :: cells
        type(spline_space_options_t) :: options
        integer :: status

        options%degree = 3
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
        if (status /= 0) error stop "space"
    end subroutine build

    real(dp) function inverse_distance(q, x) result(total)
        type(spline_quadrature_t), intent(in) :: q
        real(dp), intent(in) :: x(2)
        integer :: i

        total = 0.0_dp
        do i = 1, q%point_count
            total = total + q%weight(i)/norm2(q%s(:, i) - x)
        end do
    end function inverse_distance

    real(dp) function exact_integral(x) result(v)
        !! int over [c0 - a, c0 + a]^2 of 1/|s - x| by the antiderivative
        !! G = x ln(y + r) + y ln(x + r).
        real(dp), intent(in) :: x(2)
        real(dp) :: lo(2), hi(2)

        lo = c0 - a - x
        hi = c0 + a - x
        v = g(hi(1), hi(2)) - g(lo(1), hi(2)) - g(hi(1), lo(2)) + g(lo(1), lo(2))
    end function exact_integral

    real(dp) function g(x, y)
        real(dp), intent(in) :: x, y
        real(dp) :: r

        r = sqrt(x*x + y*y)
        g = 0.0_dp
        if (abs(x) > 1.0e-14_dp) g = g + x*log(y + r)
        if (abs(y) > 1.0e-14_dp) g = g + y*log(x + r)
    end function g

    logical function pointwise_ok() result(ok)
        type(plain_kernel_t) :: kernel
        type(csr_matrix_t) :: matrix
        real(dp), allocatable :: c(:, :), v(:, :), hv(:, :), g1(:, :), g2(:, :)
        real(dp), allocatable :: jet(:, :), cv(:)
        real(dp) :: e1, e2, one(1)
        integer :: status, i, k, q

        allocate(c(1, space%dof_count), v(1, space%dof_count), hv(1, space%dof_count), &
            g1(1, space%dof_count), g2(1, space%dof_count), cv(space%dof_count))
        call random_number(c)
        call random_number(v)
        call spline_energy(space, quad, kernel, c, e1, status, g1)
        call spline_energy(space, adapted, kernel, c, e2, status, g2)
        ok = abs(e1 - e2) < 1.0e-11_dp*abs(e1) .and. &
            maxval(abs(g1 - g2)) < 1.0e-11_dp*maxval(abs(g1))
        call spline_energy_hessian_csr(space, adapted, kernel, c, matrix, status)
        call spline_energy_hessian_apply(space, adapted, kernel, c, v, hv, status)
        do i = 1, matrix%n
            cv(i) = 0.0_dp
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                cv(i) = cv(i) + matrix%val(k)*v(1, matrix%colidx(k))
            end do
        end do
        ok = ok .and. maxval(abs(cv - hv(1, :))) < 1.0e-11_dp*maxval(abs(cv))
        allocate(jet(1, adapted%point_count))
        call spline_evaluate(space, adapted, c, 0, jet, status)
        do q = 1, adapted%point_count, 7
            call spline_evaluate_at(space, c, adapted%block(q), adapted%xi(:, q), 0, &
                one, status)
            ok = ok .and. abs(one(1) - jet(1, q)) < 1.0e-12_dp
        end do
    end function pointwise_ok

end program test_spline_singular
