module fortfem_gll_line
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_GLL
    use fortfem_generated_cardinal_product_jet1_order1, only: &
        generated_cardinal_product_jet1_order1
    implicit none
    private
    public :: build_gll_line, gll_rule, gauss_rule, lagrange_table
contains
    subroutine legendre(n, x, pn, previous)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp), intent(out) :: pn, previous
        real(dp) :: next
        integer :: k
        pn = 1.0_dp
        previous = 0.0_dp
        do k = 1, n
            next = (real(2*k - 1, dp)*x*pn - real(k - 1, dp)*previous)/real(k, dp)
            previous = pn
            pn = next
        end do
    end subroutine legendre

    subroutine gll_rule(p, x, w, status)
        integer, intent(in) :: p
        real(dp), allocatable, intent(out) :: x(:), w(:)
        integer, intent(out) :: status
        real(dp) :: z, step, pn, pm
        integer :: i, it
        status = 1
        if (p < 1) return
        allocate (x(p + 1), w(p + 1))
        x(1) = -1.0_dp
        x(p + 1) = 1.0_dp
        do i = 2, (p + 2)/2
            z = -cos(acos(-1.0_dp)*real(i - 1, dp)/real(p, dp))
            do it = 1, 100
                call legendre(p, z, pn, pm)
                step = (z*pn - pm)/(real(p + 1, dp)*pn)
                z = z - step
                if (abs(step) < epsilon(z)) exit
            end do
            if (it > 100) return
            x(i) = z
            x(p + 2 - i) = -z
        end do
        do i = 1, p + 1
            call legendre(p, x(i), pn, pm)
            w(i) = 2.0_dp/(real(p, dp)*real(p + 1, dp)*pn*pn)
        end do
        status = 0
    end subroutine gll_rule

    subroutine gauss_rule(n, x, w, status)
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: x(:), w(:)
        integer, intent(out) :: status
        real(dp) :: z, step, pn, pm, derivative
        integer :: i, it
        status = 1
        if (n < 1) return
        allocate (x(n), w(n))
        do i = 1, (n + 1)/2
            z = cos(acos(-1.0_dp)*(real(i, dp) - 0.25_dp)/(real(n, dp) + 0.5_dp))
            do it = 1, 100
                call legendre(n, z, pn, pm)
                derivative = real(n, dp)*(pm - z*pn)/(1.0_dp - z*z)
                step = pn/derivative
                z = z - step
                if (abs(step) < epsilon(z)) exit
            end do
            if (it > 100) return
            call legendre(n, z, pn, pm)
            derivative = real(n, dp)*(pm - z*pn)/(1.0_dp - z*z)
            x(i) = -z
            x(n + 1 - i) = z
            w(i) = 2.0_dp/((1.0_dp - z*z)*derivative*derivative)
            w(n + 1 - i) = w(i)
        end do
        status = 0
    end subroutine gauss_rule

    subroutine lagrange_table(nodes, points, b, d)
        real(dp), intent(in) :: nodes(:), points(:)
        real(dp), allocatable, intent(out) :: b(:, :), d(:, :)
        real(dp) :: v, dv, nv, ndv
        integer :: i, j, q
        allocate (b(size(points), size(nodes)), d(size(points), size(nodes)))
        do i = 1, size(nodes)
            do q = 1, size(points)
                v = 1.0_dp
                dv = 0.0_dp
                do j = 1, size(nodes)
                    if (j == i) cycle
                    call generated_cardinal_product_jet1_order1(v, dv, &
                           1.0_dp/(nodes(i) - nodes(j)), 0.0_dp, points(q) - nodes(j), &
                                                                nv, ndv)
                    v = nv
                    dv = ndv
                end do
                b(q, i) = v
                d(q, i) = dv
            end do
        end do
    end subroutine lagrange_table

    subroutine build_gll_line(p, nelem, line, status, periodic, nq, gll, x0, x1)
        integer, intent(in) :: p, nelem
        type(tensor_line_tables_t), intent(out) :: line
        integer, intent(out) :: status
        logical, intent(in), optional :: periodic, gll
        integer, intent(in), optional :: nq
        real(dp), intent(in), optional :: x0, x1
        real(dp), allocatable :: nodes(:), weights(:), x(:), w(:), b(:, :), d(:, :)
        real(dp) :: h
        integer :: e, i, count
        logical :: collocated
        status = 1
        if (p < 1 .or. nelem < 1) return
        if (present(x0)) line%x0 = x0
        if (present(x1)) line%x1 = x1
        if (line%x1 <= line%x0) return
        if (present(periodic)) line%periodic = periodic
        collocated = .false.
        if (present(gll)) collocated = gll
        count = p + 1
        if (present(nq)) count = nq
        if (count < 1) return
        if (collocated .and. count /= p + 1) return
        call gll_rule(p, nodes, weights, status)
        if (status /= 0) return
        if (collocated) then
            x = nodes
            w = weights
        else
            call gauss_rule(count, x, w, status)
            if (status /= 0) return
        end if
        call lagrange_table(nodes, x, b, d)
        line%family = LINE_GLL
        line%degree = p
        line%nelem = nelem
        line%nloc = p + 1
        line%nq = count
        line%ndof = p*nelem + 1
        if (line%periodic) line%ndof = p*nelem
        allocate (line%xq(count, nelem), line%wq(count, nelem))
        allocate (line%b(count, p + 1, nelem), line%d(count, p + 1, nelem))
        allocate (line%dof(p + 1, nelem))
        h = (line%x1 - line%x0)/real(nelem, dp)
        do e = 1, nelem
            line%xq(:, e) = line%x0 + (real(e - 1, dp) + (x + 1.0_dp)/2.0_dp)*h
            line%wq(:, e) = w*h/2.0_dp
            line%b(:, :, e) = b
            line%d(:, :, e) = d*2.0_dp/h
            do i = 1, p + 1
                line%dof(i, e) = modulo((e - 1)*p + i - 1, line%ndof) + 1
            end do
        end do
        status = 0
    end subroutine build_gll_line
end module fortfem_gll_line
