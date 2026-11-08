program test_meridian_bspline_gradient
    !! R and Z derivatives of B = curl_n A for a random V1 field: against
    !! central differences of meridian_bspline_eval and against
    !! div B = dB_R/dR + B_R/R + i n B_phi/R + dB_Z/dZ = 0.
    use fortfem_kinds, only: dp
    use fortfem_meridian_bspline_n, only: meridian_bspline_t, &
        meridian_bspline_build, meridian_bspline_eval, &
        meridian_bspline_eval_gradient
    implicit none

    integer, parameter :: npts = 50, n = 3
    real(dp), parameter :: h = 1.0e-5_dp
    type(meridian_bspline_t) :: sp
    real(dp), allocatable :: re(:), im(:)
    complex(dp), allocatable :: u(:)
    real(dp) :: x(2, npts), xs(2, 1), scale, err_value, err_fd, err_div
    complex(dp) :: b(3, npts), db(3, 2, npts), a1(3, 1), b1(3, 1), bp(3, 1), &
        bm(3, 1), fd(3), div
    integer :: k, d, status

    call meridian_bspline_build(2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp, 8, 8, 3, &
        sp, status)
    if (status /= 0) error stop 'build failed'
    allocate (re(sp%n1), im(sp%n1), u(sp%n1))
    call random_seed(put=[(20261108 + k, k = 1, 64)])
    call random_number(re)
    call random_number(im)
    u = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
    call random_number(x)
    x(1, :) = 2.05_dp + 0.9_dp*x(1, :)
    x(2, :) = -0.45_dp + 0.9_dp*x(2, :)
    call meridian_bspline_eval_gradient(sp, n, u, x, b, db, status)
    if (status /= 0) error stop 'eval_gradient failed'
    err_value = 0
    err_fd = 0
    err_div = 0
    do k = 1, npts
        xs(:, 1) = x(:, k)
        call meridian_bspline_eval(sp, n, u, xs, a1, b1, status)
        scale = maxval(abs(db(:, :, k)))
        err_value = max(err_value, maxval(abs(b(:, k) - b1(:, 1))) &
            /maxval(abs(b1(:, 1))))
        do d = 1, 2
            xs(:, 1) = x(:, k)
            xs(d, 1) = x(d, k) + h
            call meridian_bspline_eval(sp, n, u, xs, a1, bp, status)
            xs(d, 1) = x(d, k) - h
            call meridian_bspline_eval(sp, n, u, xs, a1, bm, status)
            fd = (bp(:, 1) - bm(:, 1))/(2*h)
            err_fd = max(err_fd, maxval(abs(db(:, d, k) - fd))/scale)
        end do
        div = db(1, 1, k) + b(1, k)/x(1, k) &
            + cmplx(0.0_dp, real(n, dp), dp)*b(2, k)/x(1, k) + db(3, 2, k)
        err_div = max(err_div, abs(div)/scale)
    end do
    call check('B equals meridian_bspline_eval', err_value, 1.0e-13_dp)
    call check('dB equals central differences', err_fd, 1.0e-7_dp)
    call check('div B vanishes', err_div, 1.0e-12_dp)

contains

    subroutine check(label, err, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: err, tol

        print '(a, es10.2, a, es9.2)', label//': ', err, ' <= ', tol
        if (.not. err <= tol) error stop 'test_meridian_bspline_gradient failed'
    end subroutine check

end program test_meridian_bspline_gradient
