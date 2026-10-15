program bem_residual_study
    !! Cost/scaling fixture for finite single-layer residual energy algebra.
    !! This synthetic SPD fixture is not a physical FEM-BEM efficiency study.
    use fortfem_kinds, only: dp
    use, intrinsic :: iso_fortran_env, only: int64
    !$ use omp_lib, only: omp_set_num_threads
    use fortfem_bem_residual_energy, only: bem_residual_factor, &
        bem_residual_energy, BEM_RESIDUAL_OK
    implicit none
    real(dp), allocatable :: v(:, :), factor(:, :), q(:), rhs(:), candidate(:), work(:)
    real(dp) :: setup, solve, elapsed, indicator, oracle
    integer :: n, i, nt, threads, status, repeat, trial
    integer(int64) :: t0, t1, rate
    interface
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
    call system_clock(count_rate=rate)
    print '(a)', 'n trial threads factor_s solve_s residual_s residual/(factor+solve)'
    do n = 512, 2048, 512
        allocate (v(n, n), q(n), rhs(n), candidate(n), work(n))
        q = [(sin(real(i, dp)), i=1, n)]
        q = q/norm2(q)
        do i = 1, n
            v(:, i) = q*q(i)
            v(i, i) = v(i, i) + 1
        end do
        rhs = q
        candidate = 0.25_dp*q
        ! V q = 2 q; hence the correction is q/4 and its energy is sqrt(1/8).
        oracle = sqrt(0.125_dp)
        do trial = 1, 3
            call system_clock(t0)
            call bem_residual_factor(v, factor, status)
            call system_clock(t1)
            if (status /= BEM_RESIDUAL_OK) error stop 'factor'
            setup = real(t1 - t0, dp)/rate
            call system_clock(t0)
            do repeat = 1, 20
                work = rhs
                call dpotrs('L', n, 1, factor, n, work, n, status)
                if (status /= 0) error stop 'solve'
            end do
            call system_clock(t1)
            solve = real(t1 - t0, dp)/(20*rate)
            do nt = 1, 4
                if (nt == 3) cycle
                threads = 1
                !$ threads = nt
                !$ call omp_set_num_threads(nt)
                call system_clock(t0)
                do repeat = 1, 20
                    call bem_residual_energy(v, factor, rhs, candidate, work, &
                        indicator, status)
                    if (status /= BEM_RESIDUAL_OK) error stop 'residual'
                end do
                call system_clock(t1)
                elapsed = real(t1 - t0, dp)/(20*rate)
                if (abs(indicator - oracle) > 1.0e-13_dp) error stop 'oracle'
                print '(3(i6,1x),4es14.5)', n, trial, threads, setup, solve, &
                    elapsed, elapsed/(setup + solve)
            end do
        end do
        deallocate (v, factor, q, rhs, candidate, work)
    end do
end program bem_residual_study
