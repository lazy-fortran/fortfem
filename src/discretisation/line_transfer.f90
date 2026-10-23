module fortfem_line_transfer
    ! One-dimensional nested transfers between line tables of the same family,
    ! for multigrid on line-table meshes: B-spline knot insertion (uniform span
    ! halving, any degree/continuity, open or periodic), Fourier mode
    ! truncation and GLL degree lowering on the same cells. Matrices map
    ! coarse line dofs to fine line dofs.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_gll_line, only: gll_rule, lagrange_table
    implicit none
    private
    public :: bspline_line_prolongation, fourier_line_prolongation
    public :: gll_line_prolongation
    interface
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface
contains
    subroutine uniform_knots(p, n, k, wrap, t, ndof)
        ! Same uniform knot vectors as build_bspline_line(nelem=n) on [0, 1].
        integer, intent(in) :: p, n, k
        logical, intent(in) :: wrap
        real(dp), allocatable, intent(out) :: t(:)
        integer, intent(out) :: ndof
        real(dp), allocatable :: core(:)
        integer :: mult, j, offset
        mult = p - k
        if (wrap) then
            ndof = n*mult
            allocate (core(ndof), t(ndof + 2*p + 1))
            do j = 1, n
                core((j - 1)*mult + 1:j*mult) = real(j - 1, dp)/real(n, dp)
            end do
            do j = 1, size(t)
                offset = j - p - 1
                t(j) = core(modulo(offset, ndof) + 1) + &
                       real(floor(real(offset, dp)/real(ndof, dp)), dp)
            end do
        else
            allocate (t(2*(p + 1) + (n - 1)*mult))
            t(:p + 1) = 0
            do j = 1, n - 1
                t(p + 2 + (j - 1)*mult:p + 1 + j*mult) = real(j, dp)/real(n, dp)
            end do
            t(size(t) - p:) = 1
            ndof = size(t) - p - 1
        end if
    end subroutine

    subroutine basis_row(t, p, ndof, wrap, x, row, status)
        real(dp), intent(in) :: t(:), x
        integer, intent(in) :: p, ndof
        logical, intent(in) :: wrap
        real(dp), intent(out) :: row(:)
        integer, intent(out) :: status
        real(dp), allocatable :: values(:), derivatives(:)
        integer :: j, col
        call evaluate_bspline_basis(t, p, x, values, derivatives, status)
        if (status /= 0) return
        row = 0
        do j = 1, size(values)
            col = j
            if (wrap) col = modulo(j - 1, ndof) + 1
            if (col > ndof) cycle
            row(col) = row(col) + values(j)
        end do
    end subroutine

    subroutine bspline_line_prolongation(degree, nelem_coarse, p, status, &
                                         continuity, periodic)
        ! Knot insertion from nelem_coarse to 2*nelem_coarse uniform spans, by an
        ! exact least-squares fit of the nested coarse basis (residual checked).
        integer, intent(in) :: degree, nelem_coarse
        real(dp), allocatable, intent(out) :: p(:, :)
        integer, intent(out) :: status
        integer, intent(in), optional :: continuity
        logical, intent(in), optional :: periodic
        real(dp), allocatable :: tf(:), tc(:), bf(:, :), bc(:, :), gram(:, :)
        integer :: k, nf, nc, m, span, i, row, info
        logical :: wrap
        real(dp) :: x
        status = 1
        k = degree - 1
        if (present(continuity)) k = continuity
        wrap = .false.
        if (present(periodic)) wrap = periodic
        if (degree < 1 .or. nelem_coarse < 1 .or. k < 0 .or. k > degree - 1) return
        call uniform_knots(degree, 2*nelem_coarse, k, wrap, tf, nf)
        call uniform_knots(degree, nelem_coarse, k, wrap, tc, nc)
        m = degree + 2
        allocate (bf(2*nelem_coarse*m, nf), bc(2*nelem_coarse*m, nc))
        row = 0
        do span = 1, 2*nelem_coarse
            do i = 1, m
                x = (real(span - 1, dp) + (real(i, dp) - 0.5_dp)/real(m, dp))/ &
                    real(2*nelem_coarse, dp)
                row = row + 1
                call basis_row(tf, degree, nf, wrap, x, bf(row, :), status)
                if (status /= 0) return
                call basis_row(tc, degree, nc, wrap, x, bc(row, :), status)
                if (status /= 0) return
            end do
        end do
        gram = matmul(transpose(bf), bf)
        p = matmul(transpose(bf), bc)
        call dposv('U', nf, nc, gram, nf, p, nf, info)
        status = 2
        if (info /= 0) return
        where (abs(p) < 1.e-14_dp) p = 0
        ! Nestedness check: the fine expansion reproduces every coarse function.
        if (maxval(abs(matmul(bf, p) - bc)) > 1.e-12_dp) return
        status = 0
    end subroutine

    subroutine gll_line_prolongation(p_fine, p_coarse, nelem, p, status, periodic)
        ! Degree lowering on the same cells (p-multigrid): coarse GLL Lagrange
        ! functions interpolated at the fine GLL nodes of every cell.
        integer, intent(in) :: p_fine, p_coarse, nelem
        real(dp), allocatable, intent(out) :: p(:, :)
        integer, intent(out) :: status
        logical, intent(in), optional :: periodic
        real(dp), allocatable :: xf(:), wf(:), xc(:), wc(:), b(:, :), d(:, :)
        integer :: nf, nc, e, i, j
        logical :: wrap
        wrap = .false.
        if (present(periodic)) wrap = periodic
        status = 1
        if (p_coarse < 1 .or. p_coarse > p_fine .or. nelem < 1) return
        call gll_rule(p_fine, xf, wf, status)
        if (status /= 0) return
        call gll_rule(p_coarse, xc, wc, status)
        if (status /= 0) return
        call lagrange_table(xc, xf, b, d)
        nf = p_fine*nelem + merge(0, 1, wrap)
        nc = p_coarse*nelem + merge(0, 1, wrap)
        allocate (p(nf, nc))
        p = 0
        do e = 1, nelem
            do i = 1, p_fine + 1
                do j = 1, p_coarse + 1
                    p(modulo((e - 1)*p_fine + i - 1, nf) + 1, &
                      modulo((e - 1)*p_coarse + j - 1, nc) + 1) = b(i, j)
                end do
            end do
        end do
        where (abs(p) < 1.e-15_dp) p = 0
        status = 0
    end subroutine

    subroutine fourier_line_prolongation(m_fine, m_coarse, p, status)
        ! Mode truncation: coarse modes (1, cos k, sin k), k <= m_coarse, inject.
        integer, intent(in) :: m_fine, m_coarse
        real(dp), allocatable, intent(out) :: p(:, :)
        integer, intent(out) :: status
        integer :: i
        status = 1
        if (m_coarse < 0 .or. m_coarse > m_fine) return
        allocate (p(2*m_fine + 1, 2*m_coarse + 1))
        p = 0
        do i = 1, 2*m_coarse + 1
            p(i, i) = 1
        end do
        status = 0
    end subroutine
end module
