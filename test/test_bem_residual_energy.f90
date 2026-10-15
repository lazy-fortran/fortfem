program test_bem_residual_energy
    !! Independent Fourier oracle: on the unit circle the Laplace single
    !! layer has eigenvalues 1/(2 k), k >= 1, in an L2-normalized cosine basis.
    !! A finite residual norm is a LOWER bound, even with exact quadrature.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_bem_residual_energy, only: bem_residual_factor, &
        bem_residual_energy, BEM_RESIDUAL_OK, BEM_RESIDUAL_BAD_INPUT, &
        BEM_RESIDUAL_NOT_SPD
    implicit none
    real(dp), allocatable :: v(:, :), factor(:, :), rhs(:), candidate(:), work(:)
    real(dp) :: indicator, exact, previous, rho, total, retained
    integer :: n, i, status, coarse

    rho = 0.8_dp
    previous = huge(1.0_dp)
    do coarse = 2, 8, 2
        n = 4*coarse
        allocate (v(n, n), rhs(n), candidate(n), work(n))
        v = 0
        do i = 1, n
            v(i, i) = 1.0_dp/(2*i)
            rhs(i) = rho**i/(2*i)
            candidate(i) = 0
            if (i <= coarse) candidate(i) = rho**i
        end do
        call bem_residual_factor(v, factor, status)
        call check_condition(status == BEM_RESIDUAL_OK, 'SPD single layer')
        call bem_residual_energy(v, factor, rhs, candidate, work, indicator, status)
        call check_condition(status == BEM_RESIDUAL_OK, 'residual energy')
        ! Infinite tail: sum rho^(2 k)/(2 k) = -log(1-rho^2)/2.
        total = -log(1 - rho**2)/2
        retained = 0
        do i = 1, coarse
            retained = retained + rho**(2*i)/(2*i)
        end do
        exact = sqrt(total - retained)
        call check_condition(indicator < exact, 'enriched norm is not an upper bound')
        call check_condition(indicator > 0.98_dp*exact, 'resolved tail efficiency')
        call check_condition(indicator < previous, 'tail decreases with resolution')
        print '(a,i3,3(a,es12.4))', 'Fourier coarse modes ', coarse, &
            ' error ', exact, ' indicator ', indicator, ' ratio ', indicator/exact
        previous = indicator
        ! Off-diagonal coverage is supplied independently below.
        candidate = [(rho**i, i=1, n)]
        call bem_residual_energy(v, factor, rhs, candidate, work, indicator, status)
        call check_condition(indicator < 1.0e-15_dp, 'exact recovery at rounding')
        deallocate (v, factor, rhs, candidate, work)
    end do

    ! Missing one boundary test direction can make a measured norm exactly
    ! zero despite a nonzero exterior error. A non-diagonal SPD example also
    ! checks that the Cholesky whitening uses the correct triangular factor.
    allocate (v(2, 2), rhs(2), candidate(2), work(2))
    v = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
    rhs = [0.0_dp, 1.0_dp]
    candidate = 0
    call bem_residual_factor(v, factor, status)
    call bem_residual_energy(v, factor, rhs, candidate, work, indicator, status)
    call check_condition(status == BEM_RESIDUAL_OK, 'non-diagonal SPD residual')
    call check_condition(abs(indicator - sqrt(2.0_dp/3)) < 1.0e-14_dp, &
        'inverse known analytically')
    call check_condition(norm2(work - [-1.0_dp/3, 2.0_dp/3]) < 1.0e-14_dp, &
        'correction solves the reference equation')
    call bem_residual_energy(v(:1, :1), factor(:1, :1), rhs(:1), &
        candidate(:1), work(:1), indicator, status)
    call check_condition(indicator < 1.0e-15_dp, 'restricted residual misses the error')
    v(1, 1) = -1
    call bem_residual_factor(v, factor, status)
    call check_condition(status == BEM_RESIDUAL_NOT_SPD, &
        'indefinite operator rejected')
    call bem_residual_factor(v(:, :1), factor, status)
    call check_condition(status == BEM_RESIDUAL_BAD_INPUT, &
        'nonsquare operator rejected')
    call check_summary('bem_residual_energy')
end program test_bem_residual_energy
