program test_parametric_bem_torus_slow
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_toroidal_mode_kernel, only: tm_green, tm_green_gradient
    use fortfem_parametric_laplace_bem_3d, only: &
        assemble_parametric_laplace_bem, build_exterior_dtn, &
        bem_operators_t, bem_exterior_dtn_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_axisymmetric_gs_boundary_operators_2d, only: &
        gs_bem_curve_t, gs_bem_operators_t, gs_bem_curve_from_nodes, &
        gs_bem_panel_point
    use fortfem_toroidal_mode_boundary_operators_2d, only: tm_bem_assemble
    use parametric_bem_test_geometry, only: test_torus_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(test_torus_t) :: torus
    type(gs_bem_curve_t) :: curve
    type(gs_bem_operators_t) :: modal_op
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    real(dp), allocatable :: vertices(:,:), trace(:), q(:), load(:), exact(:)
    real(dp), allocatable :: xi(:), eta(:), qw(:)
    real(dp) :: x(3), a(3), b(3), normal(3), grad(4), value, radius, jac
    real(dp) :: err(2), start, assembly, solve_time, charges(3,64), field_grad(3)
    real(dp), allocatable :: curve_nodes(:,:), modal_density(:)
    real(dp) :: angle, p2(2), tangent2(2), speed, energy2d, energy_error
    real(dp), allocatable :: ones(:), jump_residual(:), exact_residual(:)
    real(dp) :: jump_error
    integer :: level, nt, np, status, i, p, t, mode, count_override, degree, charge
    character(24) :: arg
    logical :: passed, deform
    passed = .true.
    count_override = 0
    if (command_argument_count() > 0) then
        call get_command_argument(1,arg)
        read(arg,*) count_override
    end if
    deform = .false.
    degree = 16
    if (command_argument_count() > 1) then
        call get_command_argument(2,arg)
        deform = trim(arg) == 'deformed'
    end if
    if (command_argument_count() > 2) then
        call get_command_argument(3,arg)
        read(arg,*) degree
    end if
    write(*,'(a,i0)') 'BEM3D quadrature_degree ', degree
    call triangle_duffy_quadrature(10, xi, eta, qw, status)
    do level = 1, 2
        np = 4 + 4*level
        if (count_override > 0) np = count_override
        nt = 2*np
        call generate_torus_surface_mesh(2.0_dp, 0.6_dp, np, nt, &
            vertices, torus%faces, torus%parameters)
        if (deform) torus%deformation(1,1) = 1.25_dp
        allocate(trace(size(vertices,2)), q(size(torus%faces,2)), &
            load(size(vertices,2)), exact(size(torus%faces,2)))
        mode = 0
        do charge = 1, size(charges,2)
            angle = 2*acos(-1.0_dp)*real(charge-1,dp)/size(charges,2)
            charges(:,charge) = matmul(torus%deformation, &
                [2*cos(angle),2*sin(angle),0.0_dp])
        end do
        do i = 1, size(vertices,2)
            if (deform) then
                call point_source_field(matmul(torus%deformation,vertices(:,i)), &
                    charges, trace(i), field_grad)
            else
                radius = norm2(vertices(1:2,i))
                trace(i) = tm_green(mode, radius, vertices(3,i), 2.0_dp, 0.0_dp)
            end if
        end do
        start = omp_get_wtime()
        call assemble_parametric_laplace_bem(torus, torus%faces, degree, op, status)
        assembly = omp_get_wtime()-start
        if (status /= 0) error stop 'general torus assembly failed'
        call check_condition(op%genus == 1, 'closed torus genus equals one')
        passed = passed .and. op%genus == 1
        if (.not. deform) then
            allocate(curve_nodes(2,2*np))
            do p = 1, 2*np
                angle = acos(-1.0_dp)*real(p-1,dp)/np
                curve_nodes(:,p) = [2.0_dp+0.6_dp*cos(angle), 0.6_dp*sin(angle)]
            end do
            call gs_bem_curve_from_nodes(curve_nodes, curve, status)
            if (status /= 0) error stop 'modal curve failed'
            call tm_bem_assemble(curve, 0, 2, modal_op, status)
            if (status /= 0) error stop 'modal operator failed'
            allocate(modal_density(modal_op%ny))
            do p = 1, np
                do t = 1, 3
                    call gs_bem_panel_point(curve,p,real(t-1,dp)/2, &
                        p2,tangent2,speed)
                    modal_density(3*(p-1)+t) = p2(1)
                end do
            end do
            energy2d = 2*acos(-1.0_dp)* &
                dot_product(modal_density,matmul(modal_op%v,modal_density))
            energy_error = abs(sum(op%v)/energy2d-1.0_dp)
            write(*,'(a,es14.6)') 'BEM3D per-n V energy relative error ', energy_error
            call check_condition(energy_error < 0.03_dp, &
                '3D single-layer form agrees with the per-n axisymmetric operator')
            passed = passed .and. energy_error < 0.03_dp
            deallocate(curve_nodes,modal_density)
        end if
        start = omp_get_wtime()
        call build_exterior_dtn(op, exterior, status)
        if (status /= 0) error stop 'general torus DtN setup failed'
        call exterior%apply(trace, q, load, status)
        solve_time = omp_get_wtime()-start
        if (status /= 0) error stop 'general torus DtN application failed'
        exact = 0.0_dp
        do p = 1, size(q)
            do t = 1, size(qw)
                call torus%evaluate(p, [xi(t),eta(t)], x, a, b, status)
                normal = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3), &
                    a(1)*b(2)-a(2)*b(1)]
                jac = norm2(normal)
                normal = normal/jac
                if (deform) then
                    call point_source_field(x, charges, value, field_grad)
                    value = dot_product(field_grad,normal)
                else
                    radius = norm2(x(1:2))
                    call tm_green_gradient(mode, radius, x(3), &
                        2.0_dp, 0.0_dp, value, grad)
                    value = grad(1)*dot_product(x(1:2)/radius,normal(1:2)) + &
                        grad(2)*normal(3)
                end if
                exact(p) = exact(p) + qw(t)*jac*value/op%area(p)
            end do
        end do
        allocate(ones(size(trace)))
        ones = 1.0_dp
        jump_residual = matmul(op%k+0.5_dp*op%mass,ones)
        jump_error = norm2(jump_residual)/norm2(matmul(op%mass,ones))
        exact_residual = matmul(op%v,exact)-matmul(op%k-0.5_dp*op%mass,trace)
        write(*,'(a,4(es14.6,1x))') &
            'BEM3D jump_defect cauchy_defect mean_q mean_exact ', &
            jump_error, norm2(exact_residual)/norm2(matmul(op%v,exact)), &
            sum(op%area*q)/sum(op%area),sum(op%area*exact)/sum(op%area)
        call check_condition(jump_error < 0.03_dp, &
            'closed torus constant jump relation')
        passed = passed .and. jump_error < 0.03_dp
        deallocate(ones,jump_residual,exact_residual)
        err(level) = sqrt(sum(op%area*(q-exact)**2)/sum(op%area*exact**2))
        write(*,'(a,3(i0,1x),3(es14.6,1x))') &
            'BEM3D threads npanels mode assembly_s setup_s rel_flux_error ', &
            omp_get_max_threads(), size(q), mode, assembly, solve_time, err(level)
        deallocate(vertices,torus%faces,torus%parameters,trace,q,load,exact)
        if (count_override > 0) exit
    end do
    if (count_override > 0) then
        call check_condition(err(1) < 0.03_dp, &
            'benchmark flux error below three percent')
        passed = passed .and. err(1) < 0.03_dp
    else
        call check_condition(err(2) < err(1), 'torus flux converges under refinement')
        passed = passed .and. err(2) < err(1)
        call check_condition(err(2) < 0.03_dp, 'torus flux error below three percent')
        passed = passed .and. err(2) < 0.03_dp
    end if
    call check_summary('General curved torus BEM against exact exterior field')
    if (.not. passed) error stop 1
contains
    subroutine point_source_field(point, sources, value, gradient)
        real(dp), intent(in) :: point(3), sources(:,:)
        real(dp), intent(out) :: value, gradient(3)
        real(dp) :: d(3), r
        integer :: j
        value = 0.0_dp
        gradient = 0.0_dp
        do j = 1, size(sources,2)
            d = point-sources(:,j)
            r = norm2(d)
            value = value + 1/r
            gradient = gradient-d/r**3
        end do
        value = value/size(sources,2)
        gradient = gradient/size(sources,2)
    end subroutine
end program
