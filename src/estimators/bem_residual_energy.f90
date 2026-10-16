module fortfem_bem_residual_energy
    !! Residual energy in a FINITE single-layer test space. If V is its SPD
    !! Galerkin matrix and r = b - V lambda, eta = ||L^{-1} r||_2, V = L L^T,
    !! is exactly ||lambda_ref - lambda||_V. For an exact continuous RHS it
    !! is a lower bound on the continuous V-energy error, NOT an upper bound.
    !! No saturation assumption, quadrature enclosure or continuous remainder
    !! is supplied here. The same primitive serves GS and toroidal-mode BEM.
    !!
    !! These are dense linear algebra operations, with no pointwise physical
    !! or quadrature leaves. LAPACK owns factorization and triangular solves;
    !! the indexed residual action is threaded. The caller retains the factor
    !! and supplies scratch, so the repeated action allocates no memory.
    use fortfem_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: bem_residual_factor, bem_residual_energy
    integer, parameter, public :: BEM_RESIDUAL_OK = 0, BEM_RESIDUAL_BAD_INPUT = 1
    integer, parameter, public :: BEM_RESIDUAL_NOT_SPD = 2
    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine
        subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
            import :: dp
            character, intent(in) :: uplo, trans, diag
            integer, intent(in) :: n, lda, incx
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: x(*)
        end subroutine
    end interface
contains
    subroutine bem_residual_factor(v, factor, status)
        !! Factor the supplied finite symmetric operator without changing it.
        !! factor must subsequently be used with the SAME immutable v.
        real(dp), intent(in) :: v(:, :)
        real(dp), allocatable, intent(out) :: factor(:, :)
        integer, intent(out) :: status
        integer :: n, info
        status = BEM_RESIDUAL_BAD_INPUT
        n = size(v, 1)
        if (n < 1) return
        if (size(v, 2) /= n) return
        if (.not. all(ieee_is_finite(v))) return
        if (any(abs(v - transpose(v)) > 0.0_dp)) return
        factor = v
        call dpotrf('L', n, factor, n, info)
        status = BEM_RESIDUAL_NOT_SPD
        if (info /= 0) return
        status = BEM_RESIDUAL_OK
    end subroutine

    subroutine bem_residual_energy(v, factor, rhs, candidate, work, indicator, status)
        !! work returns the correction V^{-1}(rhs - V candidate). The factor
        !! must be the lower Cholesky factor of v, e.g. bem_residual_factor's
        !! output. Neither factor nor v is revalidated on each action;
        !! changed matrices require a new factor.
        real(dp), contiguous, intent(in) :: v(:, :), factor(:, :)
        real(dp), intent(in) :: rhs(:), candidate(:)
        real(dp), contiguous, intent(out) :: work(:)
        real(dp), intent(out) :: indicator
        integer, intent(out) :: status
        integer, parameter :: tile_size = 64
        integer :: n, i, j, first, count
        real(dp) :: tile(tile_size)
        indicator = huge(1.0_dp)
        status = BEM_RESIDUAL_BAD_INPUT
        n = size(v, 1)
        if (n < 1) return
        if (size(v, 2) /= n) return
        if (any(shape(factor) /= [n, n])) return
        if (size(rhs) /= n .or. size(candidate) /= n .or. size(work) /= n) return
        if (.not. all(ieee_is_finite(rhs))) return
        if (.not. all(ieee_is_finite(candidate))) return
        ! A bounded row tile keeps the leftmost (contiguous) index innermost.
        !$omp parallel do default(shared) schedule(static) if(n >= 256) &
        !$omp private(i, j, count, tile)
        do first = 1, n, tile_size
            count = min(tile_size, n - first + 1)
            tile = 0
            do j = 1, n
                !$omp simd
                do i = 1, count
                    tile(i) = tile(i) + v(first + i - 1, j)*candidate(j)
                end do
            end do
            work(first:first + count - 1) = rhs(first:first + count - 1) - tile(:count)
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(work))) return
        call dtrsv('L', 'N', 'N', n, factor, n, work, 1)
        if (.not. all(ieee_is_finite(work))) return
        indicator = norm2(work)
        call dtrsv('L', 'T', 'N', n, factor, n, work, 1)
        if (.not. all(ieee_is_finite(work))) return
        if (.not. ieee_is_finite(indicator)) return
        status = BEM_RESIDUAL_OK
    end subroutine
end module fortfem_bem_residual_energy
