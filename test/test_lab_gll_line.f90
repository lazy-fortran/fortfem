program test_lab_gll_line
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_GLL
    use fortfem_gll_line, only: build_gll_line, gll_rule, lagrange_table
    implicit none
    type(tensor_line_tables_t) :: line
    real(dp), allocatable :: x(:), w(:), b(:, :), d(:, :)
    real(dp) :: exact, err
    integer :: p, k, i, status
    logical :: passed

    passed = .true.
    do p = 1, 16
        call gll_rule(p, x, w, status)
        call record(status == 0, 'GLL rule builds')
        call record(all(w > 0.0_dp), 'positive weights')
        do k = 0, 2*p - 1
            exact = real(1 + (-1)**k, dp)/real(k + 1, dp)
            call record(abs(sum(w*x**k) - exact) < 1.0e-15_dp, &
                'GLL polynomial moments through degree 2p-1')
        end do
        call lagrange_table(x, x, b, d)
        do i = 1, p + 1
            b(i, i) = b(i, i) - 1.0_dp
        end do
        call record(maxval(abs(b)) < 2.0e-14_dp, 'cardinality')
        call record(maxval(abs(sum(d, dim=2))) < 1.0e-11_dp, &
            'constant derivative vanishes')
        call build_gll_line(p, 3, line, status, x0=-2.0_dp, x1=4.0_dp)
        call record(status == 0 .and. line%family == LINE_GLL, 'line builds')
        call record(line%ndof == 3*p + 1, 'nonperiodic map size')
        call record(line%dof(p + 1, 1) == line%dof(1, 2), 'shared endpoint')
        err = 0.0_dp
        do i = 1, line%nelem
            do k = 0, 2*p + 1
                exact = ((2.0_dp*real(i, dp) - 2.0_dp)**(k + 1) - &
                    (2.0_dp*real(i, dp) - 4.0_dp)**(k + 1))/real(k + 1, dp)
                err = max(err, abs(sum(line%wq(:, i)*line%xq(:, i)**k) - &
                    exact)/max(1.0_dp, abs(exact)))
            end do
        end do
        call record(err < 3.0e-14_dp, 'Gauss moments on physical cells')
        call build_gll_line(p, 3, line, status, periodic=.true., gll=.true.)
        call record(line%ndof == 3*p, 'periodic map size')
        call record(line%dof(p + 1, 3) == 1, 'periodic endpoint identified')
    end do
    call build_gll_line(0, 3, line, status)
    call record(status /= 0, 'reject degree zero')
    call build_gll_line(2, 0, line, status)
    call record(status /= 0, 'reject empty mesh')
    call check_summary()
    if (.not. passed) error stop 'GLL line tests failed'
contains
    subroutine record(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        call check_condition(ok, message)
        passed = passed .and. ok
    end subroutine record
end program test_lab_gll_line
