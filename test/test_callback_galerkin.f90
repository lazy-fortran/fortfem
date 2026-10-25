module callback_galerkin_test_actions
    !! Module callback avoids a compiler-specific executable-stack trampoline.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    real(dp), public :: a(4, 4), bad_value
    integer, public :: calls = 0, mode = 0
    public :: action
contains
    subroutine action(x, y, local_status)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: local_status
        calls = calls + 1
        y = 0.0_dp
        local_status = 17
        if (mode == 1 .and. calls == 2) return
        local_status = 0
        y = matmul(a, x)
        if (mode == 2) y(1) = bad_value
    end subroutine action
end module callback_galerkin_test_actions

program test_callback_galerkin
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                             ieee_positive_inf
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use check, only: check_condition, check_summary
    use callback_galerkin_test_actions, only: a, calls, mode, bad_value, action
    use fortfem_callback_galerkin, only: galerkin_cache_t, galerkin_report_t, &
        build_callback_galerkin, apply_galerkin_correction, &
        callback_galerkin_ready, callback_galerkin_matrix, &
        free_callback_galerkin, callback_galerkin_storage, GALERKIN_OK, &
        GALERKIN_INVALID_SHAPE, GALERKIN_BUDGET_EXCEEDED, &
        GALERKIN_NONFINITE_INPUT, GALERKIN_CALLBACK_FAILED, &
        GALERKIN_NONFINITE_ACTION, GALERKIN_NONSYMMETRIC, &
        GALERKIN_NOT_POSITIVE, GALERKIN_NOT_READY
    implicit none

    integer(int64), parameter :: budget = 64_int64*1024_int64*1024_int64
    real(dp), parameter :: tol = 5.0e-14_dp
    type(galerkin_cache_t) :: cache
    type(galerkin_report_t) :: report
    real(dp) :: p(4, 2), coarse(2, 2), rhs(4), z(4), zs(4)
    real(dp) :: s(4), expected(4), elapsed, nan, infinity, wrong(3), wrong_mat(1, 1)
    real(dp) :: strided_rhs(8), strided_z(8)
    integer(int64) :: peak, persistent
    integer :: status, i

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    bad_value = nan
    call apply_galerkin_correction(cache, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
                                   z, status)
    call check_condition(status == GALERKIN_NOT_READY .and. all(z == 0.0_dp), &
                         'unbuilt cache returns a defined zero correction')
    call callback_galerkin_storage(4, 2, peak, persistent, status)
    call check_condition(status == GALERKIN_OK .and. peak == 176_int64 .and. &
                         persistent == 176_int64, 'complete array payload accounting')
    call callback_galerkin_storage(huge(i), huge(i), peak, persistent, status)
    call check_condition(status == GALERKIN_BUDGET_EXCEEDED, &
                         'oversized storage arithmetic rejects without overflow')

    ! Independent rational oracle: det(A[1:3,1:3])=18 and det(Ac)=26.
    ! P has columns e1+e3, e2+e3, and fixed fourth row zero.
    a = 0.0_dp
    a(1, :) = [4.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    a(2, :) = [1.0_dp, 3.0_dp, 1.0_dp, 0.0_dp]
    a(3, :) = [0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp]
    a(4, 4) = 1.0_dp
    p(:, 1) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    p(:, 2) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call check_condition(status == GALERKIN_OK .and. callback_galerkin_ready(cache), &
                         'independent rational SPD cache builds')
    call check_condition(calls == 2 .and. report%callback_count == 2 .and. &
                         report%failed_column == 0 .and. report%callback_status == 0, &
                         'setup calls the frozen action once per ordered column')
    call check_condition(report%symmetry_checked .and. &
                         report%symmetry_defect == 0.0_dp .and. &
                         report%factor_attempted .and. report%factor_info == 0, &
                         'raw symmetry and native Cholesky statuses are recorded')
    call check_condition(report%setup_seconds >= 0.0_dp .and. &
                         report%peak_payload_bytes == 176_int64 .and. &
                         report%persistent_payload_bytes == 176_int64, &
                         'setup timing and peak/persistent payload are reported')
    call callback_galerkin_matrix(cache, coarse, status)
    call check_condition(status == GALERKIN_OK .and. &
                         maxval(abs(coarse - reshape([6.0_dp, 4.0_dp, &
                                    4.0_dp, 7.0_dp], [2, 2]))) < tol, &
                         'raw Ac agrees with the independent rational matrix')
    p = 0.0_dp
    rhs = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    expected = [4.0_dp, 7.0_dp, 11.0_dp, 0.0_dp]/13.0_dp
    call apply_galerkin_correction(cache, rhs, z, status, elapsed)
    call check_condition(status == GALERKIN_OK .and. elapsed >= 0.0_dp .and. &
                         maxval(abs(z - expected)) < tol, &
                         'owned P survives caller mutation and exact result agrees')
    call check_condition(calls == 2, 'cache application makes no fine-action call')
    strided_rhs = -99.0_dp
    strided_z = -99.0_dp
    strided_rhs(1:8:2) = rhs
    call apply_galerkin_correction(cache, strided_rhs(1:8:2), &
                                   strided_z(1:8:2), status)
    call check_condition(status == GALERKIN_OK .and. &
                         maxval(abs(strided_z(1:8:2) - expected)) < tol .and. &
                         all(strided_z(2:8:2) == -99.0_dp), &
                         'strided vectors preserve untouched entries')
    s = [2.0_dp, -1.0_dp, 0.5_dp, 7.0_dp]
    call apply_galerkin_correction(cache, s, zs, status)
    call check_condition(status == GALERKIN_OK .and. &
                         abs(dot_product(rhs, zs) - dot_product(s, z)) < tol .and. &
                         dot_product(rhs, z) > 0.0_dp, &
                         'rational correction has symmetry and positive energy')
    call apply_galerkin_correction(cache, [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
                                   z, status)
    call check_condition(status == GALERKIN_OK .and. all(z == 0.0_dp), &
                         'explicit fixed-row-only input gives zero correction')
    call apply_galerkin_correction(cache, rhs(:3), z, status)
    call check_condition(status == GALERKIN_INVALID_SHAPE .and. all(z == 0.0_dp), &
                         'wrong input dimension fails with zero output')
    call apply_galerkin_correction(cache, rhs, wrong, status)
    call check_condition(status == GALERKIN_INVALID_SHAPE .and. &
                         all(wrong == 0.0_dp), 'wrong output dimension is rejected')
    rhs(2) = nan
    call apply_galerkin_correction(cache, rhs, z, status)
    call check_condition(status == GALERKIN_NONFINITE_INPUT .and. all(z == 0.0_dp), &
                         'nonfinite residual fails without an identity fallback')
    call check_condition(callback_galerkin_ready(cache) .and. calls == 2, &
                         'bad residual leaves the immutable factor reusable')
    call callback_galerkin_matrix(cache, wrong_mat, status)
    call check_condition(status == GALERKIN_INVALID_SHAPE .and. &
                         all(wrong_mat == 0.0_dp), 'audit matrix requires exact shape')

    call unit_problem()
    call build_callback_galerkin(p, action, cache, 175_int64, status, report)
    call expect_failure(GALERKIN_BUDGET_EXCEEDED, 'budget below full peak is rejected')
    call check_condition(calls == 0, 'budget admission precedes callback work')
    call build_callback_galerkin(p, action, cache, 0_int64, status, report)
    call expect_failure(GALERKIN_BUDGET_EXCEEDED, 'nonpositive budget is rejected')
    call build_callback_galerkin(p, action, cache, 176_int64, status, report)
    call check_condition(status == GALERKIN_OK .and. calls == 2, &
                         'exact full peak budget admits the cache')
    call build_callback_galerkin(p(:1, :), action, cache, budget, status, report)
    call expect_failure(GALERKIN_INVALID_SHAPE, 'failed rebuild clears old ready state')
    p(1, 1) = nan
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NONFINITE_INPUT, 'NaN P is rejected')
    p(1, 1) = infinity
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NONFINITE_INPUT, 'infinite P is rejected')

    call unit_problem()
    p(:, 2) = p(:, 1)
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NOT_POSITIVE, 'duplicate columns reject singular Ac')
    call check_condition(report%factor_attempted .and. report%factor_info == 2, &
                         'singular factor retains its native pivot failure')
    call unit_problem()
    a(1, 1) = -1.0_dp
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NOT_POSITIVE, &
                        'indefinite Ac is not shifted or repaired')
    call check_condition(report%factor_info == 1, &
                         'indefinite factor retains its native pivot failure')
    call unit_problem()
    a(1, 2) = 0.1_dp
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NONSYMMETRIC, 'nonsymmetric action is not averaged')
    call check_condition(abs(report%symmetry_defect - 0.1_dp) < tol .and. &
                         .not. report%factor_attempted, &
                         'scaled raw asymmetry rejects before native factorization')
    call unit_problem()
    a(1, 2) = 0.5e-12_dp
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call check_condition(status == GALERKIN_OK, &
                         'asymmetry below the declared scaled guard is admitted')
    call callback_galerkin_matrix(cache, coarse, status)
    call check_condition(status == GALERKIN_OK .and. coarse(1, 2) == a(1, 2) .and. &
                         coarse(2, 1) == 0.0_dp, 'audit retains unaveraged raw entries')
    rhs = [2.0_dp, 3.0_dp, 0.0_dp, 0.0_dp]
    call apply_galerkin_correction(cache, rhs, z, status)
    call check_condition(status == GALERKIN_OK .and. &
                         maxval(abs(z - rhs)) < tol, &
                         'native solve uses the declared lower triangle')
    call unit_problem()
    a(1, 2) = 2.0e-12_dp
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NONSYMMETRIC, 'asymmetry above 1e-12 is rejected')

    call unit_problem()
    mode = 1
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_CALLBACK_FAILED, 'callback failure is not retried')
    call check_condition(calls == 2 .and. report%callback_count == 2 .and. &
                         report%failed_column == 2 .and. report%callback_status == 17, &
                         'original callback status and failing column survive')
    call unit_problem()
    mode = 2
    call build_callback_galerkin(p, action, cache, budget, status, report)
    call expect_failure(GALERKIN_NONFINITE_ACTION, 'nonfinite action is rejected')
    call check_condition(calls == 1 .and. report%failed_column == 1 .and. &
                         report%callback_status == 0, &
                         'successful callback with bad output remains distinguishable')
    call callback_galerkin_matrix(cache, coarse, status)
    call check_condition(status == GALERKIN_NOT_READY .and. all(coarse == 0.0_dp), &
                         'failed cache exposes no partial audit matrix')
    call free_callback_galerkin(cache)
    call free_callback_galerkin(cache)
    call check_condition(.not. callback_galerkin_ready(cache), &
                         'cache reset is idempotent')
    call check_summary('callback Galerkin coarse cache')

contains

    subroutine unit_problem()
        integer :: k
        a = 0.0_dp
        do k = 1, 4
            a(k, k) = 1.0_dp
        end do
        p = 0.0_dp
        p(1, 1) = 1.0_dp
        p(2, 2) = 1.0_dp
        calls = 0
        mode = 0
    end subroutine unit_problem

    subroutine expect_failure(expected_status, label)
        integer, intent(in) :: expected_status
        character(*), intent(in) :: label
        call check_condition(status == expected_status .and. &
                             .not. callback_galerkin_ready(cache), label)
        call apply_galerkin_correction(cache, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
                                       z, status)
        call check_condition(status == GALERKIN_NOT_READY .and. all(z == 0.0_dp), &
                             'failure leaves no usable partial correction')
    end subroutine expect_failure

end program test_callback_galerkin
