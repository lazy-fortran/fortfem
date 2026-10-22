program test_lab_fem_bem_transmission
    !! FEM interior + BEM exterior transmission on the through-axis O-grid
    !! torus (5-block O-grid x phi, axis interior). Exact oracle in
    !! fem_bem_oracle; errors must decrease under p- and h-refinement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_get_wtime
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                       tensor_load, tensor_l2_error, &
                                       apply_tensor_operator
    use fortfem_tensor_solver, only: tensor_solve_info_t
    use fortfem_tensor_fem_bem, only: tensor_fem_bem_t, build_tensor_fem_bem, &
                                      tensor_fem_bem_options_t, &
                                      solve_tensor_fem_bem, tensor_fem_bem_rhs, &
                                      apply_tensor_fem_bem, &
                                      tensor_fem_bem_exterior_flux, &
                                      free_tensor_fem_bem, &
                                      tensor_h1_seminorm_error
    use fortfem_parametric_bem_tensor, only: bem_tensor_flux_l2_error
    use fem_bem_oracle
    implicit none
    integer :: failures = 0
    real(dp) :: e(4, 4)
    integer :: its(4)
    type(tensor_fem_bem_options_t) :: options
    character(32) :: arg
    integer :: values(6), k
    real(dp) :: schur
    print '(a)', 'FB1_TRANSMISSION n,nphi,p,pf,dofs,trace,flux,its,relL2,relH1,'// &
        'relq,relq_dtn,t_bem,t_dtn,t_mg,t_chol,t_solve,t_apply_fem,'// &
        't_apply_steklov,peak_rss_kib'
    if (command_argument_count() >= 3) then
        ! Scan mode: n nphi p [flux_degree target_points fan_points schur_scale]
        values = [0, 0, 0, -1, 0, 0]
        do k = 1, min(6, command_argument_count())
            call get_command_argument(k, arg)
            read (arg, *) values(k)
        end do
        options%flux_degree = values(4)
        options%quadrature%target_points = values(5)
        options%quadrature%fan_points = values(6)
        if (command_argument_count() >= 7) then
            call get_command_argument(7, arg)
            read (arg, *) schur
            options%schur_scale = schur
        end if
        call case(values(1), values(2), values(3), e(:, 1), its(1), options)
        stop
    end if
    call case(1, 6, 2, e(:, 1), its(1), options)
    call case(1, 6, 3, e(:, 2), its(2), options)
    call case(1, 6, 4, e(:, 3), its(3), options)
    call case(2, 12, 2, e(:, 4), its(4), options)
    call check(all(e(:, 2) < e(:, 1)) .and. all(e(:, 3) < e(:, 2)), &
               'p-refinement reduces volume L2/H1 and Neumann errors')
    call check(all(e(:, 4) < e(:, 1)), 'h-refinement reduces all errors')
    call check(e(1, 3) < 5.e-3_dp, 'p=4 relative volume L2 error below 5e-3')
    call check(maxval(its) <= 30, 'block-preconditioned PCG iterations bounded')
    if (failures /= 0) error stop 'test_lab_fem_bem_transmission failed'
    print '(a)', 'test_lab_fem_bem_transmission passed'
contains
    integer function peak_rss_kib()
        character(256) :: line
        integer :: unit, ios
        peak_rss_kib = -1
        open (newunit=unit, file='/proc/self/status', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:6) == 'VmHWM:') read (line(7:), *) peak_rss_kib
        end do
        close (unit)
    end function

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) then
            print '(a,a)', '[PASS] ', label
        else
            print '(a,a)', '[FAIL] ', label
            failures = failures + 1
        end if
    end subroutine

    subroutine case(n, nphi, p, errors, iterations, options)
        !! errors: relative volume L2, H1 seminorm, coupled Neumann, and
        !! DtN-only Neumann error from the interpolated exact exterior trace.
        integer, intent(in) :: n, nphi, p
        real(dp), intent(out) :: errors(4)
        integer, intent(out) :: iterations
        type(tensor_fem_bem_options_t), intent(in) :: options
        type(tensor_mesh_t) :: mesh
        type(tensor_operator_t) :: err_op
        type(tensor_fem_bem_t) :: fb
        type(tensor_solve_info_t) :: info
        real(dp), allocatable :: load(:), rhs(:), x(:), trace(:), flux(:), zero(:)
        real(dp) :: unorm, gnorm, qnorm, solve_time, t_fem, t_coupled
        real(dp), allocatable :: y(:)
        integer :: status, i, nboundary
        errors = huge(1.0_dp)
        iterations = huge(1)
        call ogrid_tensor_mesh(3, n, p, mesh, status, nphi=nphi)
        if (status /= 0) error stop 'mesh'
        nboundary = count(mesh%boundary)
        mesh%boundary = .false.
        call build_tensor_operator(mesh, unit_coefficient, err_op, status, nq=p + 3)
        if (status /= 0) error stop 'error operator'
        call build_tensor_fem_bem(mesh, unit_coefficient, fb, status, options)
        if (status /= 0) then
            print '(a,i0)', 'build_tensor_fem_bem status ', status
            error stop
        end if
        ! The O-grid outer ring is the only boundary: genus-one, 4 n nphi
        ! panels, and exactly the volume nodes the mesh marks as boundary.
        if (fb%trace%genus /= 1 .or. fb%trace%npanel /= 4*n*nphi .or. &
            fb%trace%ntrace /= nboundary) then
            print '(a,3i8)', 'trace topology ', fb%trace%genus, fb%trace%npanel, &
                fb%trace%ntrace
            failures = failures + 1
        end if
        call tensor_load(err_op, interior_forcing, load)
        call tensor_fem_bem_rhs(fb, load, rhs, status, jump_neumann, jump_dirichlet)
        if (status /= 0) error stop 'rhs'
        allocate (x(mesh%ndof), zero(mesh%ndof))
        x = 0
        zero = 0
        solve_time = omp_get_wtime()
        call solve_tensor_fem_bem(fb, rhs, 1.e-11_dp, x, info)
        solve_time = omp_get_wtime() - solve_time
        if (info%status /= 0) print '(a,es10.3)', 'PCG residual ', info%residual
        iterations = info%iterations
        unorm = tensor_l2_error(err_op, zero, interior_u)
        gnorm = tensor_h1_seminorm_error(err_op, zero, interior_gradient)
        errors(1) = tensor_l2_error(err_op, x, interior_u)/unorm
        errors(2) = tensor_h1_seminorm_error(err_op, x, interior_gradient)/gnorm
        allocate (trace(fb%trace%ntrace))
        do i = 1, fb%trace%ntrace
            trace(i) = x(fb%trace%trace_dof(i)) - &
                       jump_dirichlet(mesh%coordinates(:, fb%trace%trace_dof(i)))
        end do
        call tensor_fem_bem_exterior_flux(fb, trace, flux, status)
        if (status /= 0) error stop 'flux'
        call bem_tensor_flux_l2_error(fb%trace%surface, fb%space, flux, &
                                      exterior_normal, errors(3), status)
        flux = 0
        call bem_tensor_flux_l2_error(fb%trace%surface, fb%space, flux, &
                                      exterior_normal, qnorm, status)
        errors(3) = errors(3)/qnorm
        do i = 1, fb%trace%ntrace
            trace(i) = exterior_u(mesh%coordinates(:, fb%trace%trace_dof(i)))
        end do
        call tensor_fem_bem_exterior_flux(fb, trace, flux, status)
        if (status /= 0) error stop 'flux'
        call bem_tensor_flux_l2_error(fb%trace%surface, fb%space, flux, &
                                      exterior_normal, errors(4), status)
        errors(4) = errors(4)/qnorm
        ! Average one matrix-free FEM apply and one coupled apply; their
        ! difference is the dense trace-space Steklov action.
        allocate (y(mesh%ndof))
        t_fem = omp_get_wtime()
        do i = 1, 5
            call apply_tensor_operator(fb%mg%level(1)%op, x, y, eliminate=.false.)
        end do
        t_fem = (omp_get_wtime() - t_fem)/5
        t_coupled = omp_get_wtime()
        do i = 1, 5
            call apply_tensor_fem_bem(fb, x, y)
        end do
        t_coupled = (omp_get_wtime() - t_coupled)/5
        print '(a,8(i0,a),4(es10.3,a),7(es10.3,a),i0)', 'FB1_TRANSMISSION ', n, &
            ',', nphi, ',', p, ',', fb%space%flux_degree, ',', mesh%ndof, ',', &
            fb%trace%ntrace, ',', fb%space%nflux, ',', info%iterations, ',', &
            errors(1), ',', errors(2), ',', errors(3), ',', errors(4), ',', &
            fb%time_bem, ',', fb%time_dtn, ',', fb%time_mg, ',', fb%time_schur, &
            ',', solve_time, ',', t_fem, ',', t_coupled - t_fem, ',', peak_rss_kib()
        call free_tensor_fem_bem(fb)
    end subroutine
end program
