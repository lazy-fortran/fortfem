program test_parametric_bem_convergence_slow
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, &
        bem_exterior_dtn_t, assemble_parametric_laplace_bem, &
        build_exterior_dtn, parametric_bem_refinement_estimate, &
        parametric_bem_residual
    use parametric_bem_test_geometry, only: test_sphere_t
    implicit none
    type(test_sphere_t) :: sphere
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    real(dp), allocatable :: trace(:), flux(:), load(:), old_flux(:)
    real(dp), allocatable :: indicators(:), residual(:), xi(:), eta(:), qw(:)
    integer, allocatable :: parent(:)
    real(dp) :: l2_error(0:2), marking(2), energy(2), old_weighted, effectivity
    real(dp) :: x(3), a(3), b(3), normal(3), integral, jacobian
    real(dp) :: panel_error
    real(dp), allocatable :: old_error(:), old_area(:)
    integer :: level, status, panel, t, degree, ios
    character(32) :: argument
    logical :: passed

    passed = .true.
    degree = 16
    call get_command_argument(1, argument)
    if (len_trim(argument) > 0) then
        read (argument, *, iostat=ios) degree
        if (ios /= 0) error stop 'invalid audit quadrature degree'
        if (degree < 2) error stop 'audit quadrature degree must be at least two'
    end if
    write (*, '(a,i0)') 'BEM3D estimator quadrature_degree ', degree
    call triangle_duffy_quadrature(24, xi, eta, qw, status)
    if (status /= 0) error stop 'audit quadrature failed'
    do level = 0, 2
        call generate_sphere_surface_mesh(1.0_dp, level, &
            sphere%vertices, sphere%faces)
        call assemble_parametric_laplace_bem(sphere, sphere%faces, degree, op, status)
        if (status /= 0) error stop 'audit assembly failed'
        call build_exterior_dtn(op, exterior, status)
        if (status /= 0) error stop 'audit DtN setup failed'
        allocate(trace(size(sphere%vertices, 2)), load(size(sphere%vertices, 2)))
        allocate(flux(size(op%area)), residual(size(op%area)))
        ! Exact exterior dipole u=z/r**3: u=z and q=-2*z on r=1.
        trace = sphere%vertices(3, :)
        call exterior%apply(trace, flux, load, status)
        if (status /= 0) error stop 'audit DtN action failed'
        call parametric_bem_residual(op, trace, flux, residual, status)
        call record(status == 0 .and. norm2(residual) < 1.e-10_dp, &
            'each solved mesh has negligible equation defect')
        integral = 0.0_dp
        if (level > 0) then
            allocate(parent(size(flux)), indicators(size(old_flux)))
            parent = [( (panel+3)/4, panel=1,size(flux) )]
            call parametric_bem_refinement_estimate(op, old_flux, flux, &
                parent, indicators, energy(level), status)
            call record(status == 0, 'refined estimator accepts solved pair')
            if (status /= 0) error stop 'audit estimator failed'
            marking(level) = norm2(indicators)
            old_weighted = sqrt(sum(sqrt(old_area)*old_error))
            effectivity = marking(level)/old_weighted
            write (*, '(a,i0,4(1x,es14.6))') &
                'BEM3D estimator level energy marking weighted_error effectivity ', &
                level, energy(level), marking(level), old_weighted, effectivity
            call record(energy(level) > 1.e-4_dp, &
                'refinement detects error hidden by the equation defect')
            deallocate(parent, indicators, old_error, old_area)
        end if
        allocate(old_error(size(flux)))
        do panel = 1, size(flux)
            panel_error = 0.0_dp
            do t = 1, size(qw)
                call sphere%evaluate(panel, [xi(t), eta(t)], x, a, b, status)
                if (status /= 0) error stop 'audit geometry failed'
                normal = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
                    a(1)*b(2)-a(2)*b(1)]
                jacobian = norm2(normal)
                panel_error = panel_error + qw(t)*jacobian*(flux(panel)+2*x(3))**2
                integral = integral + qw(t)*jacobian*4*x(3)**2
            end do
            old_error(panel) = panel_error
        end do
        l2_error(level) = sqrt(sum(old_error)/integral)
        write (*, '(a,2(i0,1x),es14.6)') 'BEM3D continuous_flux level panels error ', &
            level, size(flux), l2_error(level)
        old_flux = flux
        old_area = op%area
        deallocate(trace, flux, load, residual)
    end do
    call record(l2_error(1) < l2_error(0) .and. l2_error(2) < l2_error(1), &
        'continuous flux error converges including within-panel variation')
    call record(marking(2) < marking(1), 'local refinement estimator converges')
    call record(energy(2) < energy(1), 'single-layer energy estimator converges')
    call check_summary('Parametric BEM dipole convergence audit')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
