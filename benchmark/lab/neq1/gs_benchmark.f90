program fortfem_lab_neq1
    ! Usage: fortfem_lab_neq1 case.nml n p samples.txt out_prefix [rtol] [rays]
    !        [profiles.txt npoints]
    ! Solves the free-axis Grad-Shafranov problem on an n x n per-block O-grid
    ! of degree p and prints one key=value summary row. Timed phases: setup
    ! (mesh, operator, p-MG, locator), initial guess, Newton. Sampling, errors
    ! and safety-factor contours are post-processing and reported separately.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    use fortfem_newton_krylov, only: newton_options_t, newton_info_t, &
                                     newton_krylov_solve
    use fortfem_tensor_operator, only: tensor_l2_error
    use fortfem_lab_gs, only: gs_case_t, gs_system_t, read_gs_case, solovev_psi, &
                              build_gs_system, gs_initial_guess, gs_point, &
                              gs_safety_factor
    !$ use omp_lib, only: omp_get_max_threads
    implicit none
    type(gs_case_t) :: case
    type(gs_system_t) :: system
    type(newton_options_t) :: options
    type(newton_info_t) :: info
    character(512) :: path, samples, prefix, arg, profiles
    real(dp), allocatable :: u(:), zero(:)
    real(dp) :: t0, t_setup, t_guess, t_newton, t_post, error, norm, rtol, &
                x(2), psi, gradient(2), psin, fsq, b(3), pressure, q, &
                exact_axis_error, psia_error
    integer :: n, p, status, guess_its, rays, unit, out, ios, cell, k, threads, &
               npoints
    call get_command_argument(1, path)
    call get_command_argument(2, arg); read (arg, *) n
    call get_command_argument(3, arg); read (arg, *) p
    call get_command_argument(4, samples)
    call get_command_argument(5, prefix)
    rtol = 1.e-10_dp
    rays = 256
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg); read (arg, *) rtol
    end if
    if (command_argument_count() >= 7) then
        call get_command_argument(7, arg); read (arg, *) rays
    end if
    profiles = ''
    npoints = 201
    if (command_argument_count() >= 8) call get_command_argument(8, profiles)
    if (command_argument_count() >= 9) then
        call get_command_argument(9, arg); read (arg, *) npoints
    end if
    threads = 1
    !$ threads = omp_get_max_threads()
    call read_gs_case(trim(path), case, status)
    if (status /= 0) error stop 'cannot read case namelist'
    t0 = wall()
    call build_gs_system(case, n, p, system, status)
    if (status /= 0) error stop 'cannot build Grad-Shafranov system'
    t_setup = wall() - t0
    allocate (u(system%op%mesh%ndof), zero(system%op%mesh%ndof))
    t0 = wall()
    call gs_initial_guess(system, u, 1.e-3_dp, guess_its, status)
    if (status /= 0) error stop 'initial guess failed'
    t_guess = wall() - t0
    options%rtol = rtol
    options%atol = 1.e-14_dp
    t0 = wall()
    call newton_krylov_solve(system, u, options, info)
    t_newton = wall() - t0
    t0 = wall()
    zero = 0
    error = tensor_l2_error(system%op, u, solovev_psi)
    norm = tensor_l2_error(system%op, zero, solovev_psi)
    exact_axis_error = norm2(system%axis_x - [case%r_axis, case%z_axis])
    psia_error = abs(system%psia - case%psi_axis)/abs(case%psi_axis)
    open (newunit=unit, file=trim(samples), status='old', action='read')
    open (newunit=out, file=trim(prefix)//'samples_out.txt', status='replace', &
          action='write')
    write (out, '(a)') '# R Z psi B_R B_Z B_phi absB p iota'
    cell = 0
    do
        read (unit, '(a)', iostat=ios) arg
        if (ios /= 0) exit
        if (len_trim(arg) == 0) cycle
        if (arg(1:1) == '#') cycle
        read (arg, *) x
        call gs_point(system, u, x, cell, psi, gradient, status)
        if (status /= 0) then
            write (out, '(2es25.16e3,a)') x, ' NaN NaN NaN NaN NaN NaN NaN'
            cycle
        end if
        psin = (psi - system%psia)/(system%psib - system%psia)
        fsq = case%f0**2 + case%cf*(1 - psin)**case%alpha_f
        b = [-gradient(2)/x(1), gradient(1)/x(1), sqrt(fsq)/x(1)]
        pressure = case%p_axis*(1 - psin)**case%alpha_p
        call gs_safety_factor(system, u, psi, rays, q, status, &
                              start=norm2(x - system%axis_x))
        if (status /= 0) q = huge(1.0_dp)
        write (out, '(9es25.16e3)') x, psi, b, norm2(b), pressure, 1/q
    end do
    close (unit)
    close (out)
    if (len_trim(profiles) > 0) call write_profiles(trim(profiles), npoints)
    t_post = wall() - t0
    write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'n=', n, ' p=', p, &
        ' dofs=', system%op%mesh%ndof, ' cells=', system%op%mesh%nelem, &
        ' threads=', threads, ' status=', info%status
    write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') 'newton=', info%newton, &
        ' krylov=', info%krylov, ' residuals=', info%residuals, ' backtracks=', &
        info%backtracks, ' guess_pcg=', guess_its, ' axis_searches=', &
        system%axis_searches
    write (output_unit, '(a,es12.5,a,es12.5)') 'residual0=', info%initial, &
        ' residual=', info%residual
    write (output_unit, '(a,*(es12.5,1x))') 'history=', &
        info%history(0:info%newton)
    write (output_unit, '(a,*(i0,1x))') 'krylov_history=', &
        info%krylov_history(1:info%newton)
    write (output_unit, '(a,es12.5,a,es12.5,a,es12.5,a,es12.5)') 't_setup=', &
        t_setup, ' t_guess=', t_guess, ' t_newton=', t_newton, ' t_solve=', &
        t_setup + t_guess + t_newton
    write (output_unit, '(a,5(es12.5,1x),a,es12.5)') &
        't_residual_linearise_jac_prec_krylov=', info%time, ' t_axis=', &
        system%time_axis
    write (output_unit, '(a,es12.5,a,es12.5,a,es12.5,a,es12.5)') 'rel_l2=', &
        error/norm, ' axis_error=', exact_axis_error, ' psia_rel_error=', &
        psia_error, ' t_post=', t_post
    write (output_unit, '(a,es25.16e3,a,2es25.16e3)') 'psia=', system%psia, &
        ' axis=', system%axis_x
contains
    subroutine write_profiles(file, m)
        ! psiN s p iota on a uniform psiN grid; q by contour integration at
        ! interior levels, cubic extrapolation to axis and edge, and
        ! Phi = 2 pi int q dpsi by fourth-order cumulative quadrature.
        character(*), intent(in) :: file
        integer, intent(in) :: m
        real(dp) :: qv(m), phi(m), psin_k(m), h, level, span
        integer :: k, unit, st
        span = system%psib - system%psia
        h = 1.0_dp/real(m - 1, dp)
        do k = 1, m
            psin_k(k) = h*real(k - 1, dp)
        end do
        do k = 2, m - 1
            level = system%psia + psin_k(k)*span
            call gs_safety_factor(system, u, level, rays, qv(k), st)
            if (st /= 0) error stop 'profile contour failed'
        end do
        qv(1) = 4*qv(2) - 6*qv(3) + 4*qv(4) - qv(5)
        qv(m) = 4*qv(m - 1) - 6*qv(m - 2) + 4*qv(m - 3) - qv(m - 4)
        phi(1) = 0
        phi(2) = h/24*(9*qv(1) + 19*qv(2) - 5*qv(3) + qv(4))
        do k = 2, m - 2
            phi(k + 1) = phi(k) + h/24*(-qv(k - 1) + 13*qv(k) + 13*qv(k + 1) - &
                                        qv(k + 2))
        end do
        phi(m) = phi(m - 1) + h/24*(9*qv(m) + 19*qv(m - 1) - 5*qv(m - 2) + &
                                    qv(m - 3))
        phi = 2*acos(-1.0_dp)*span*phi
        open (newunit=unit, file=file, status='replace', action='write')
        write (unit, '(a,es25.16e3)') '# phiedge = ', phi(m)
        write (unit, '(a)') '# lab FortFEM NEQ1 profile table: q by contour '// &
            'integration, Phi = 2*pi*int q dpsi'
        write (unit, '(a)') '# columns: psiN s p[Pa] iota ; rows uniform in psiN'
        do k = 1, m
            write (unit, '(4es25.16e3)') psin_k(k), phi(k)/phi(m), &
                case%p_axis*(1 - psin_k(k))**case%alpha_p, 1/qv(k)
        end do
        close (unit)
    end subroutine

    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function
end program fortfem_lab_neq1
