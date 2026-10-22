module fortfem_bspline_lines
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_BSPLINE
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: build_bspline_line, build_bspline_line_pair
contains
    subroutine build_bspline_line(degree, line, status, nelem, knots, &
            continuity, periodic, nq)
        integer, intent(in) :: degree
        type(tensor_line_tables_t), intent(out) :: line
        integer, intent(out) :: status
        integer, optional, intent(in) :: nelem, continuity, nq
        real(dp), optional, intent(in) :: knots(:)
        logical, optional, intent(in) :: periodic
        real(dp), allocatable :: t(:)
        integer :: ndof, order
        logical :: wrap

        call prepare_knots(degree, t, ndof, wrap, status, nelem, knots, &
            continuity, periodic)
        if (status /= 0) return
        order = degree + 1
        if (present(nq)) order = nq
        if (order < 1) then
            status = 1
            return
        end if
        call tabulate(t, degree, ndof, wrap, order, line, status)
    end subroutine build_bspline_line

    subroutine build_bspline_line_pair(degree, v0, v1, incidence, status, &
            nelem, knots, continuity, periodic, nq)
        integer, intent(in) :: degree
        type(tensor_line_tables_t), intent(out) :: v0, v1
        real(dp), allocatable, intent(out) :: incidence(:, :)
        integer, intent(out) :: status
        integer, optional, intent(in) :: nelem, continuity, nq
        real(dp), optional, intent(in) :: knots(:)
        logical, optional, intent(in) :: periodic
        real(dp), allocatable :: t(:)
        real(dp) :: scale
        integer :: ndof, nedge, order, j, e, i, row, col
        logical :: wrap

        status = 1
        if (degree < 1) return
        call prepare_knots(degree, t, ndof, wrap, status, nelem, knots, &
            continuity, periodic)
        if (status /= 0) return
        order = degree + 1
        if (present(nq)) order = nq
        status = 1
        if (order < 1) return
        call tabulate(t, degree, ndof, wrap, order, v0, status)
        if (status /= 0) return
        nedge = ndof - 1
        if (wrap) nedge = ndof
        call tabulate(t(2:size(t) - 1), degree - 1, nedge, wrap, order, &
            v1, status)
        if (status /= 0) return
        ! M_j = p B_{j,p-1} / (t_{j+p+1}-t_{j+1}). The coefficients
        ! of d u are differences u_{j+1}-u_j, including the periodic seam.
        do e = 1, v1%nelem
            do i = 1, v1%nloc
                j = v1%dof(i, e)
                scale = real(degree, dp)/(t(j + degree + 1) - t(j + 1))
                v1%b(:, i, e) = scale*v1%b(:, i, e)
                v1%d(:, i, e) = scale*v1%d(:, i, e)
            end do
        end do
        allocate(incidence(nedge, ndof))
        incidence = 0.0_dp
        do row = 1, nedge
            col = row + 1
            if (wrap) col = modulo(col - 1, ndof) + 1
            incidence(row, row) = incidence(row, row) - 1.0_dp
            incidence(row, col) = incidence(row, col) + 1.0_dp
        end do
    end subroutine build_bspline_line_pair

    subroutine prepare_knots(p, t, ndof, wrap, status, nelem, knots, &
            continuity, periodic)
        integer, intent(in) :: p
        real(dp), allocatable, intent(out) :: t(:)
        integer, intent(out) :: ndof, status
        logical, intent(out) :: wrap
        integer, optional, intent(in) :: nelem, continuity
        real(dp), optional, intent(in) :: knots(:)
        logical, optional, intent(in) :: periodic
        real(dp), allocatable :: breaks(:), core(:)
        real(dp) :: length
        integer :: k, mult, n, j, i, offset, count_equal

        status = 1
        wrap = .false.
        ndof = 0
        if (p < 0) return
        if (present(periodic)) wrap = periodic
        k = p - 1
        if (present(continuity)) k = continuity
        if (k < min(0, p - 1) .or. k > p - 1) return
        mult = p - k
        if (present(nelem) .eqv. present(knots)) return
        if (present(knots)) then
            if (size(knots) < 2) return
            if (.not. all(ieee_is_finite(knots))) return
            if (any(knots(2:) < knots(:size(knots) - 1))) return
            if (.not. wrap) then
                n = size(knots)
                if (n < 2*(p + 1)) return
                if (knots(p + 1) >= knots(n - p)) return
                if (any(knots(:p + 1) > knots(1))) return
                if (any(knots(n - p:) < knots(n))) return
                if (knots(p + 2) <= knots(1)) return
                if (knots(n - p - 1) >= knots(n)) return
                j = p + 2
                do while (j < n - p)
                    i = j + 1
                    do while (i < n - p)
                        if (knots(i) > knots(j)) exit
                        i = i + 1
                    end do
                    count_equal = i - j
                    if (count_equal > max(1, p)) return
                    if (present(continuity)) then
                        if (count_equal /= mult) return
                    end if
                    j = i
                end do
                t = knots
                ndof = n - p - 1
                status = 0
                return
            end if
            ! Periodic input is the strictly increasing span break sequence.
            if (any(knots(2:) <= knots(:size(knots) - 1))) return
            breaks = knots
            n = size(breaks) - 1
        else
            if (nelem < 1) return
            n = nelem
            allocate(breaks(n + 1))
            do j = 1, n + 1
                breaks(j) = real(j - 1, dp)/real(n, dp)
            end do
        end if
        if (wrap) then
            ndof = n*mult
            allocate(core(ndof), t(ndof + 2*p + 1))
            do j = 1, n
                core((j - 1)*mult + 1:j*mult) = breaks(j)
            end do
            length = breaks(n + 1) - breaks(1)
            do j = 1, size(t)
                offset = j - p - 1
                t(j) = core(modulo(offset, ndof) + 1) &
                    + real(floor(real(offset, dp)/real(ndof, dp)), dp)*length
            end do
        else
            allocate(t(2*(p + 1) + (n - 1)*mult))
            t(:p + 1) = breaks(1)
            do j = 1, n - 1
                t(p + 2 + (j - 1)*mult:p + 1 + j*mult) = breaks(j + 1)
            end do
            t(size(t) - p:) = breaks(n + 1)
            ndof = size(t) - p - 1
        end if
        status = 0
    end subroutine prepare_knots

    subroutine tabulate(t, p, ndof, wrap, nq, line, status)
        real(dp), intent(in) :: t(:)
        integer, intent(in) :: p, ndof, nq
        logical, intent(in) :: wrap
        type(tensor_line_tables_t), intent(out) :: line
        integer, intent(out) :: status
        real(dp), allocatable :: values(:), derivatives(:)
        integer :: nbase, ncells, span, e, q, i, j

        nbase = size(t) - p - 1
        ncells = count(t(p + 2:nbase + 1) > t(p + 1:nbase))
        line%family = LINE_BSPLINE
        line%degree = p
        line%nelem = ncells
        line%nloc = p + 1
        line%nq = nq
        line%ndof = ndof
        line%periodic = wrap
        line%x0 = t(p + 1)
        line%x1 = t(nbase + 1)
        allocate(line%xq(nq, ncells), line%wq(nq, ncells), &
            line%b(nq, p + 1, ncells), line%d(nq, p + 1, ncells), &
            line%dof(p + 1, ncells))
        e = 0
        do span = p + 1, nbase
            if (t(span + 1) <= t(span)) cycle
            e = e + 1
            call gauss_legendre_ab(nq, t(span), t(span + 1), &
                line%xq(:, e), line%wq(:, e))
            do i = 1, p + 1
                j = span - p + i - 1
                line%dof(i, e) = j
                if (wrap) line%dof(i, e) = modulo(j - 1, ndof) + 1
            end do
            do q = 1, nq
                call evaluate_bspline_basis(t, p, line%xq(q, e), &
                    values, derivatives, status)
                if (status /= 0) return
                line%b(q, :, e) = values(span - p:span)
                line%d(q, :, e) = derivatives(span - p:span)
            end do
        end do
        status = 0
    end subroutine tabulate
end module fortfem_bspline_lines
