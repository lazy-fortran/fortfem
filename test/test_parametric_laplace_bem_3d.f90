program test_parametric_laplace_bem_3d
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_parametric_laplace_bem_3d, only: &
        assemble_parametric_laplace_bem, build_exterior_dtn, &
        bem_operators_t, bem_exterior_dtn_t, parametric_bem_residual
    use parametric_bem_test_geometry, only: test_sphere_t
    implicit none
    type(test_sphere_t) :: sphere
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    real(dp), allocatable :: u(:), q(:), residual(:), load(:)
    real(dp), allocatable :: matrix(:,:)
    real(dp) :: error, defect
    integer :: status, n
    logical :: passed

    passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 1, sphere%vertices, sphere%faces)
    call assemble_parametric_laplace_bem(sphere, sphere%faces, 8, op, status)
    call record(status == 0, 'general curved sphere assembly succeeds')
    if (status /= 0) error stop 1
    n = size(sphere%vertices, 2)
    allocate(u(n), q(size(sphere%faces, 2)), residual(size(sphere%faces, 2)), load(n))
    u = 1.0_dp
    defect = norm2(matmul(op%k + 0.5_dp*op%mass, u))/norm2(matmul(op%mass, u))
    write (*, '(a,es12.4)') 'constant double-layer defect: ', defect
    call record(defect < 0.03_dp, 'K 1 = -1/2 on a closed sphere')
    call record(maxval(abs(op%v - transpose(op%v))) < 1.e-12_dp, 'V symmetric')
    call record(maxval(abs(op%w - transpose(op%w))) < 1.e-12_dp, 'W symmetric')
    call record(norm2(matmul(op%w, u)) < 1.e-10_dp, 'W annihilates constants')
    call record(maxval(abs(op%adjoint - transpose(op%k))) < 1.e-14_dp, &
        'adjoint is the Galerkin transpose')
    call build_exterior_dtn(op, exterior, status)
    call record(status == 0, 'fixed-shell DtN setup succeeds')
    if (status /= 0) error stop 1
    call exterior%apply(u, q, load, status)
    call record(status == 0, 'cached exterior DtN applies')
    error = norm2(q + 1.0_dp)/sqrt(real(size(q), dp))
    write (*, '(a,es12.4)') 'sphere monopole flux error: ', error
    call record(error < 0.03_dp, 'exterior 1/r has outward derivative -1')
    residual = matmul(op%v, q) - matmul(op%k - 0.5_dp*op%mass, u)
    call record(norm2(residual) < 1.e-11_dp, 'first exterior Calderon equation')
    call record(dot_product(u, load) > 0.0_dp, 'exterior energy is positive')
    call record(abs(dot_product(u, load) - 4.0_dp*acos(-1.0_dp)) < 0.1_dp, &
        'symmetric Steklov map reproduces sphere capacity')
    call exterior%matrix(matrix, status)
    call record(status == 0, 'cached Steklov matrix is available to FEM coupling')
    call record(norm2(matmul(matrix,u)-load) < 1.e-12_dp, &
        'matrix and apply return the same exterior load')
    call parametric_bem_residual(op,u,q,residual,status)
    call record(status == 0, 'equation-defect indicator accepts valid data')
    call record(norm2(residual) < 1.e-11_dp, 'solved finite equation has zero defect')
    q = q + 0.1_dp
    call parametric_bem_residual(op,u,q,residual,status)
    call record(norm2(residual) > 0.01_dp, 'indicator detects a perturbed flux')
    call check_summary('Parametric 3D Galerkin BEM')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
