program test_verified_rt0_majorant
    use iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, operator(+), operator(/)
    use fortfem_verified_rt0_majorant, only: verified_rt0_majorant
    implicit none
    real(dp) :: points(2, 4), recovered(3, 2), coarse_upper
    integer :: triangles(3, 2), status, mode
    type(interval_t) :: bound, flux, residual

    points = reshape([0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                      1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4])
    triangles = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    mode = 1
    ! sigma=(x,y), represented by its exact outward edge integrals.
    recovered(:, 1) = [0.0_dp, 1.0_dp, 0.0_dp]
    recovered(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
    call verified_rt0_majorant(points, triangles, recovered, interval(1.0_dp), &
                               problem, bound, flux, residual, status, 0)
    if (status /= 0) error stop 'constant residual audit failed'
    ! With grad u=0, f=-2, ||sigma||^2=2/3 and residual zero.
    if (flux%lo > 2.0_dp/3 .or. flux%hi < 2.0_dp/3) &
        error stop 'polynomial integral not enclosed'
    if (residual%lo > 0 .or. residual%hi > 1.e-25_dp) &
        error stop 'exact divergence not enclosed tightly'
    if (bound%hi < sqrt(2.0_dp/3)) error stop 'bound too small'
    coarse_upper = flux%hi
    call verified_rt0_majorant(points, triangles, recovered, interval(1.0_dp), &
                               problem, bound, flux, residual, status, 3)
    if (status /= 0) error stop 'refined audit failed'
    if (flux%lo > 2.0_dp/3 .or. flux%hi < 2.0_dp/3) &
        error stop 'refined polynomial integral not enclosed'
    if (flux%hi >= coarse_upper) error stop 'subdivision did not tighten enclosure'
    ! Deliberately inconsistent second-cell edge moments: shared values must
    ! be canonicalized, rather than implicitly assumed to match in arithmetic.
    recovered(1, 2) = 1.e-12_dp
    call verified_rt0_majorant(points, triangles, recovered, interval(1.0_dp), &
                               problem, bound, flux, residual, status, 1)
    if (status /= 0) error stop 'shared edge canonicalization'
    if (residual%hi > 1.e-25_dp) error stop 'shared edge not canonicalized'
    mode = 2
    recovered = 0
    call verified_rt0_majorant(points, triangles, recovered, interval(1.0_dp), &
                               problem, bound, flux, residual, status, 3)
    ! Independent antiderivative: int_0^1 (1+x)^-2 dx = 1/2.
    if (status /= 0) error stop 'rational source audit'
    if (residual%lo > 0.5_dp .or. residual%hi < 0.5_dp) &
        error stop 'rational source integral not enclosed'
    if (bound%hi < sqrt(0.5_dp)) error stop 'rational norm below exact'
    mode = 3
    call verified_rt0_majorant(points, triangles, recovered, interval(1.0_dp), &
                               problem, bound, flux, residual, status, 1)
    if (status == 0) error stop 'nonpositive coefficient certified'
    print '(a)', 'PASS: verified RT0 polynomial containment and shared edges'
contains
    subroutine problem(cell, x, y, coefficient, source, gradient, ok)
        integer, intent(in) :: cell
        type(interval_t), intent(in) :: x, y
        type(interval_t), intent(out) :: coefficient, source, gradient(2)
        logical, intent(out) :: ok
        coefficient = interval(1.0_dp)
        source = interval(-2.0_dp)
        if (mode == 2) source = interval(1.0_dp)/(interval(1.0_dp) + x)
        if (mode == 3) coefficient = interval(-1.0_dp, 1.0_dp)
        gradient = interval(0.0_dp)
        ok = cell > 0 .and. x%lo <= x%hi .and. y%lo <= y%hi
    end subroutine problem
end program test_verified_rt0_majorant
