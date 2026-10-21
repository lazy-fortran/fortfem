program test_parametric_bem_high_order
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, &
        bem_exterior_dtn_t, build_exterior_dtn
    use fortfem_parametric_bem_high_order, only: build_p2_trace_connectivity, &
        assemble_parametric_bem_p2, bem_pair_options_t, bem_pair_report_t, &
        integrate_parametric_bem_p2_pair
    use parametric_bem_test_geometry, only: test_sphere_t
    implicit none
    type(test_sphere_t) :: sphere
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    integer, allocatable :: trace_dofs(:, :), flux_dofs(:, :)
    real(dp), allocatable :: u(:), q(:), load(:)
    real(dp) :: v(3,3), k(3,6), w(6,6), correction(3,6)
    integer :: status, n
    logical :: passed
    passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, sphere%vertices, sphere%faces)
    call build_p2_trace_connectivity(sphere%faces, trace_dofs, flux_dofs, status)
    call record(status == 0, 'closed oriented P2 connectivity succeeds')
    if (status /= 0) error stop 1
    n = maxval(trace_dofs)
    call record(n == 18, 'octahedron has six vertex and twelve shared edge DOFs')
    call record(maxval(flux_dofs) == 24, 'P1 flux DOFs remain discontinuous')
    options%relative_tolerance = 0.03_dp
    options%absolute_tolerance = 1.e-7_dp
    options%initial_degree = 8
    options%maximum_degree = 32
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, op, &
        trace_dofs, flux_dofs, report, status)
    write(*,'(a,3i8,es16.7)') 'P2 quadrature status failed degree ratio ', &
        status,report%failed_pairs,report%maximum_degree_used,report%error_ratio
    call record(status == 0, 'tolerance-controlled curved P2/P1 assembly succeeds')
    if (status /= 0) error stop 1
    call record(all(shape(op%v) == [24,24]), 'V acts on discontinuous P1 fluxes')
    call record(all(shape(op%k) == [24,18]), 'K couples P1 fluxes to P2 traces')
    allocate(u(n), load(n), q(24))
    u = 1.0_dp
    call record(norm2(matmul(op%k+0.5_dp*op%mass,u)) < 1.e-11_dp, &
        'constant-subtracted K preserves the exact closed-surface identity')
    call record(norm2(matmul(op%w,u)) < 1.e-10_dp, &
        'P2 hypersingular form annihilates constants')
    call build_exterior_dtn(op, exterior, status)
    call record(status == 0, 'P2/P1 exterior DtN factorization succeeds')
    if (status /= 0) error stop 1
    call exterior%apply(u,q,load,status)
    call record(status == 0, 'P2/P1 exterior action succeeds')
    write(*,'(a,2es16.7)') 'BEM3D P2 capacity and flux error ', &
        sum(load), norm2(q+1.0_dp)/sqrt(real(size(q),dp))
    call record(abs(sum(load)/(4*acos(-1.0_dp))-1) < 0.03_dp, &
        'curved sphere monopole capacity has physical accuracy')
    call record(norm2(q+1.0_dp)/sqrt(real(size(q),dp)) < 0.05_dp, &
        'curved sphere monopole flux has physical accuracy')
    ! A cap failure must not masquerade as an accepted quadrature result.
    options%maximum_evaluations = 1
    call integrate_parametric_bem_p2_pair(sphere, 1, 1, options, &
        v,k,w,correction,report,status)
    call record(status == 3, 'kernel-work cap returns explicit nonconvergence')
    options%maximum_evaluations = 100000000
    options%maximum_degree = options%initial_degree
    call integrate_parametric_bem_p2_pair(sphere, 1, 2, options, &
        v,k,w,correction,report,status)
    call record(status == 3, 'one rule cannot claim a convergence estimate')
    options%relative_tolerance = -1.0_dp
    call integrate_parametric_bem_p2_pair(sphere, 1, 2, options, &
        v,k,w,correction,report,status)
    call record(status == 1, 'invalid tolerances are rejected')
    sphere%faces(1:2,1) = sphere%faces(2:1:-1,1)
    call build_p2_trace_connectivity(sphere%faces,trace_dofs,flux_dofs,status)
    call record(status /= 0, 'inconsistently oriented geometry is rejected')
    call check_summary('Curved P2/P1 BEM')
    if (.not. passed) error stop 1
contains
    subroutine record(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok,label)
    end subroutine
end program
