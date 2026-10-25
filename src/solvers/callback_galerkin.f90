module fortfem_callback_galerkin
    !! Cached P (P^T A P)^(-1) P^T for a caller-supplied frozen linear action.
    !! P already contains the caller's constraints and coordinate ordering.
    !! Setup owns P and raw Ac, calls A sequentially once per coarse column,
    !! and factors the admitted raw lower triangle without averaging or repair.
    !! Apply never calls A. Its cache-owned workspace requires a single owner:
    !! concurrent applications need separate caches or external serialization.
    !! The explicit byte cap covers all real array payloads, including setup
    !! and apply work. Descriptors, allocator and LAPACK/BLAS internal storage
    !! belong to the separately guarded process footprint, not this payload cap.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use omp_lib, only: omp_get_wtime
    implicit none
    private

    integer, parameter, public :: GALERKIN_OK = 0, GALERKIN_INVALID_SHAPE = 1
    integer, parameter, public :: GALERKIN_BUDGET_EXCEEDED = 2
    integer, parameter, public :: GALERKIN_ALLOCATION_FAILED = 3
    integer, parameter, public :: GALERKIN_NONFINITE_INPUT = 4
    integer, parameter, public :: GALERKIN_CALLBACK_FAILED = 5
    integer, parameter, public :: GALERKIN_NONFINITE_ACTION = 6
    integer, parameter, public :: GALERKIN_NONSYMMETRIC = 7
    integer, parameter, public :: GALERKIN_NOT_POSITIVE = 8
    integer, parameter, public :: GALERKIN_NOT_READY = 9
    integer, parameter, public :: GALERKIN_SOLVE_FAILED = 10
    real(dp), parameter, public :: GALERKIN_SYMMETRY_TOLERANCE = 1.0e-12_dp

    type, public :: galerkin_cache_t
        private
        real(dp), allocatable :: p(:, :), raw(:, :), factor(:, :), work(:, :)
        real(dp), allocatable :: fine_work(:)
        logical :: ready = .false.
    end type galerkin_cache_t

    type, public :: galerkin_report_t
        integer :: fine_size = 0, coarse_size = 0, callback_count = 0
        integer :: failed_column = 0, callback_status = 0, factor_info = 0
        integer(int64) :: peak_payload_bytes = 0_int64
        integer(int64) :: persistent_payload_bytes = 0_int64
        logical :: symmetry_checked = .false., factor_attempted = .false.
        real(dp) :: symmetry_defect = 0.0_dp, setup_seconds = 0.0_dp
    end type galerkin_report_t

    abstract interface
        subroutine galerkin_fine_action(x, y, status)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
            integer, intent(out) :: status
        end subroutine galerkin_fine_action
    end interface

    public :: galerkin_fine_action, build_callback_galerkin
    public :: apply_galerkin_correction, free_callback_galerkin
    public :: callback_galerkin_ready, callback_galerkin_matrix
    public :: callback_galerkin_storage

    interface
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
            real(dp), intent(inout) :: y(*)
        end subroutine dgemv
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
    end interface

contains

    subroutine build_callback_galerkin(p, apply_fine, cache, budget_bytes, &
                                       status, report)
        !! A failed rebuild clears any previous cache. The callback must give
        !! one frozen deterministic linear action; the API does not infer it.
        real(dp), intent(in) :: p(:, :)
        procedure(galerkin_fine_action) :: apply_fine
        type(galerkin_cache_t), intent(inout) :: cache
        integer(int64), intent(in) :: budget_bytes
        integer, intent(out) :: status
        type(galerkin_report_t), intent(out) :: report
        real(dp) :: start

        start = omp_get_wtime()
        call free_callback_galerkin(cache)
        report = galerkin_report_t()
        call build_cache(p, apply_fine, cache, budget_bytes, status, report)
        if (status /= GALERKIN_OK) call free_callback_galerkin(cache)
        report%setup_seconds = omp_get_wtime() - start
    end subroutine build_callback_galerkin

    subroutine build_cache(p, apply_fine, cache, budget_bytes, status, report)
        real(dp), intent(in) :: p(:, :)
        procedure(galerkin_fine_action) :: apply_fine
        type(galerkin_cache_t), intent(inout) :: cache
        integer(int64), intent(in) :: budget_bytes
        integer, intent(out) :: status
        type(galerkin_report_t), intent(inout) :: report
        real(dp) :: scale, defect
        integer :: n, m, i, j, allocation_status

        n = size(p, 1)
        m = size(p, 2)
        report%fine_size = n
        report%coarse_size = m
        call callback_galerkin_storage(n, m, report%peak_payload_bytes, &
                                      report%persistent_payload_bytes, status)
        if (status /= GALERKIN_OK) return
        status = GALERKIN_BUDGET_EXCEEDED
        if (budget_bytes <= 0_int64) return
        if (report%peak_payload_bytes > budget_bytes) return
        status = GALERKIN_NONFINITE_INPUT
        if (.not. all(ieee_is_finite(p))) return
        allocate(cache%p(n, m), cache%raw(m, m), cache%factor(m, m), &
                 cache%work(m, 1), cache%fine_work(n), stat=allocation_status)
        status = GALERKIN_ALLOCATION_FAILED
        if (allocation_status /= 0) return
        cache%p = p
        cache%raw = 0.0_dp
        cache%work = 0.0_dp
        do j = 1, m
            cache%fine_work = 0.0_dp
            report%callback_count = report%callback_count + 1
            call apply_fine(cache%p(:, j), cache%fine_work, report%callback_status)
            report%failed_column = j
            status = GALERKIN_CALLBACK_FAILED
            if (report%callback_status /= 0) return
            status = GALERKIN_NONFINITE_ACTION
            if (.not. all(ieee_is_finite(cache%fine_work))) return
            call dgemv('T', n, m, 1.0_dp, cache%p, n, cache%fine_work, 1, 0.0_dp, &
                       cache%raw(:, j), 1)
            if (.not. all(ieee_is_finite(cache%raw(:, j)))) return
        end do
        report%failed_column = 0
        scale = max(1.0_dp, maxval(abs(cache%raw)))
        defect = 0.0_dp
        do j = 1, m
            do i = 1, j - 1
                defect = max(defect, abs(cache%raw(i, j)/scale - &
                                         cache%raw(j, i)/scale))
            end do
        end do
        report%symmetry_defect = defect
        report%symmetry_checked = .true.
        status = GALERKIN_NONSYMMETRIC
        if (defect > GALERKIN_SYMMETRY_TOLERANCE) return
        cache%factor = cache%raw
        report%factor_attempted = .true.
        call dpotrf('L', m, cache%factor, m, report%factor_info)
        status = GALERKIN_NOT_POSITIVE
        if (report%factor_info /= 0) return
        if (.not. all(ieee_is_finite(cache%factor))) return
        cache%ready = .true.
        status = GALERKIN_OK
    end subroutine build_cache

    subroutine apply_galerkin_correction(cache, r, z, status, wall_seconds)
        type(galerkin_cache_t), intent(inout) :: cache
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: wall_seconds
        real(dp) :: start

        start = omp_get_wtime()
        call apply_cache(cache, r, z, status)
        if (present(wall_seconds)) wall_seconds = omp_get_wtime() - start
    end subroutine apply_galerkin_correction

    subroutine apply_cache(cache, r, z, status)
        type(galerkin_cache_t), intent(inout) :: cache
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        integer :: n, m, info

        z = 0.0_dp
        status = GALERKIN_NOT_READY
        if (.not. cache%ready) return
        n = size(cache%p, 1)
        m = size(cache%p, 2)
        status = GALERKIN_INVALID_SHAPE
        if (size(r) /= n .or. size(z) /= n) return
        status = GALERKIN_NONFINITE_INPUT
        if (.not. all(ieee_is_finite(r))) return
        ! Copy through owned contiguous work so strided caller vectors never
        ! require unbudgeted BLAS argument-copy arrays.
        cache%fine_work = r
        call dgemv('T', n, m, 1.0_dp, cache%p, n, cache%fine_work, 1, 0.0_dp, &
                   cache%work(:, 1), 1)
        status = GALERKIN_SOLVE_FAILED
        if (.not. all(ieee_is_finite(cache%work))) return
        call dpotrs('L', m, 1, cache%factor, m, cache%work, m, info)
        if (info /= 0) return
        if (.not. all(ieee_is_finite(cache%work))) return
        call dgemv('N', n, m, 1.0_dp, cache%p, n, cache%work(:, 1), 1, &
                   0.0_dp, cache%fine_work, 1)
        if (.not. all(ieee_is_finite(cache%fine_work))) return
        z = cache%fine_work
        status = GALERKIN_OK
    end subroutine apply_cache

    subroutine callback_galerkin_matrix(cache, matrix, status)
        !! Copy the unaveraged, unfactored audit matrix into caller storage.
        type(galerkin_cache_t), intent(in) :: cache
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        matrix = 0.0_dp
        status = GALERKIN_NOT_READY
        if (.not. cache%ready) return
        status = GALERKIN_INVALID_SHAPE
        if (any(shape(matrix) /= shape(cache%raw))) return
        matrix = cache%raw
        status = GALERKIN_OK
    end subroutine callback_galerkin_matrix

    pure logical function callback_galerkin_ready(cache) result(ready)
        type(galerkin_cache_t), intent(in) :: cache
        ready = cache%ready
    end function callback_galerkin_ready

    subroutine free_callback_galerkin(cache)
        type(galerkin_cache_t), intent(inout) :: cache
        cache%ready = .false.
        if (allocated(cache%p)) deallocate(cache%p)
        if (allocated(cache%raw)) deallocate(cache%raw)
        if (allocated(cache%factor)) deallocate(cache%factor)
        if (allocated(cache%work)) deallocate(cache%work)
        if (allocated(cache%fine_work)) deallocate(cache%fine_work)
    end subroutine free_callback_galerkin

    pure subroutine callback_galerkin_storage(n, m, peak_bytes, &
                                              persistent_bytes, status)
        integer, intent(in) :: n, m
        integer(int64), intent(out) :: peak_bytes, persistent_bytes
        integer, intent(out) :: status
        integer(int64) :: count, fine, coarse, bytes_per_real
        logical :: ok

        peak_bytes = 0_int64
        persistent_bytes = 0_int64
        status = GALERKIN_INVALID_SHAPE
        if (m < 1 .or. n < m) return
        status = GALERKIN_BUDGET_EXCEEDED
        fine = int(n, int64)
        coarse = int(m, int64)
        bytes_per_real = int(storage_size(0.0_dp)/8, int64)
        count = 0_int64
        call add_payload(count, fine, coarse, ok)
        if (.not. ok) return
        call add_payload(count, coarse, coarse, ok)
        if (.not. ok) return
        call add_payload(count, coarse, coarse, ok)
        if (.not. ok) return
        call add_payload(count, coarse, 1_int64, ok)
        if (.not. ok) return
        if (count > huge(count)/bytes_per_real) return
        call add_payload(count, fine, 1_int64, ok)
        if (.not. ok) return
        if (count > huge(count)/bytes_per_real) return
        peak_bytes = count*bytes_per_real
        persistent_bytes = peak_bytes
        status = GALERKIN_OK
    end subroutine callback_galerkin_storage

    pure subroutine add_payload(count, rows, columns, ok)
        integer(int64), intent(inout) :: count
        integer(int64), intent(in) :: rows, columns
        logical, intent(out) :: ok

        ok = .false.
        if (rows > (huge(count) - count)/columns) return
        count = count + rows*columns
        ok = .true.
    end subroutine add_payload

end module fortfem_callback_galerkin
