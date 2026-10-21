program test_parametric_bem_p2_convergence_slow
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, &
        bem_exterior_dtn_t,build_exterior_dtn
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
        bem_pair_report_t,assemble_parametric_bem_p2
    use parametric_bem_test_geometry, only: test_sphere_t
    use omp_lib, only: omp_get_wtime,omp_set_num_threads
    implicit none
    type(test_sphere_t) :: sphere
    type(bem_operators_t) :: op
    type(bem_exterior_dtn_t) :: exterior
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    integer, allocatable :: trace_dofs(:,:),flux_dofs(:,:)
    real(dp), allocatable :: u(:),q(:),load(:),xi(:),eta(:),qw(:)
    real(dp) :: errors(1:3,2),energy_error,exact_energy,order,start,elapsed
    real(dp) :: x(3),a(3),b(3),normal(3),basis(3),g(2),h(2,2),value,integral
    real(dp), parameter :: nodes(2,6)=reshape( &
        [0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp, &
        0.5_dp,0.0_dp,0.5_dp,0.5_dp,0.0_dp,0.5_dp],[2,6])
    integer :: level,l,p,t,i,status,threads
    logical :: passed
    passed=.true.
    options%relative_tolerance=1.e-3_dp
    options%absolute_tolerance=1.e-9_dp
    options%maximum_degree=48
    call triangle_duffy_quadrature(32,xi,eta,qw,status)
    if (status/=0) error stop 'audit rule'
    call omp_set_num_threads(8)
    do level=1,3
        call generate_sphere_surface_mesh(1.0_dp,level,sphere%vertices,sphere%faces)
        call assemble_parametric_bem_p2(sphere,sphere%faces,options,op, &
            trace_dofs,flux_dofs,report,status)
        write(*,'(a,4i9,3es16.7)') 'BEM3D P2 assembly level status failed degree ', &
            level,status,report%failed_pairs,report%maximum_degree_used, &
            report%error_ratio,report%raw_v_symmetry,report%raw_w_symmetry
        if (status/=0) error stop 'P2 audit assembly'
        call build_exterior_dtn(op,exterior,status)
        if (status/=0) error stop 'P2 audit factorization'
        allocate(u(maxval(trace_dofs)),q(maxval(flux_dofs)),load(maxval(trace_dofs)))
        do l=1,2
            u=0
            do p=1,size(sphere%faces,2)
                do i=1,6
                    call sphere%evaluate(p,nodes(:,i),x,a,b,status)
                    if (status/=0) error stop 'P2 nodal geometry'
                    u(trace_dofs(i,p))=harmonic(l,x(3))
                end do
            end do
            call exterior%apply(u,q,load,status)
            if (status/=0) error stop 'P2 audit action'
            integral=0
            do p=1,size(sphere%faces,2)
                do t=1,size(qw)
                    call sphere%evaluate(p,[xi(t),eta(t)],x,a,b,status)
                    if (status/=0) error stop 'P2 audit geometry'
                    normal=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3), &
                        a(1)*b(2)-a(2)*b(1)]
                    do i=1,3
                        call generated_p1_jet(i,xi(t),eta(t),basis(i),g,h)
                    end do
                    value=dot_product(q(flux_dofs(:,p)),basis)+(l+1)*harmonic(l,x(3))
                    integral=integral+qw(t)*norm2(normal)*value**2
                end do
            end do
            exact_energy=(l+1)*4*acos(-1.0_dp)/real(2*l+1,dp)
            errors(level,l)=sqrt(integral/(exact_energy*(l+1)))
            energy_error=abs(dot_product(u,load)/exact_energy-1)
            write(*,'(a,3i7,3es17.8)') 'BEM3D P2 harmonic level l panels flux energy ', &
                level,l,size(sphere%faces,2),errors(level,l),energy_error, &
                dot_product(u,load)
            call record(energy_error<0.15_dp,'spherical harmonic energy has correct scale')
            if (level==3) then
                call record(errors(level,l)<0.03_dp,'refined P1 flux reaches physical accuracy')
                call record(energy_error<0.03_dp,'refined P2 energy reaches physical accuracy')
            end if
        end do
        deallocate(u,q,load)
    end do
    do l=1,2
        order=log(errors(2,l)/errors(3,l))/log(2.0_dp)
        write(*,'(a,i3,es17.8)') 'BEM3D P2 observed_flux_order ',l,order
        call record(all(errors(2:3,l)<errors(1:2,l)), &
            'continuous spherical harmonic flux error decreases')
        call record(order>1.2_dp,'observed high-order flux convergence exceeds P0')
    end do
    ! Strong scaling uses the same geometry, tolerance and pair work at each T.
    call generate_sphere_surface_mesh(1.0_dp,1,sphere%vertices,sphere%faces)
    do i=0,3
        threads=2**i
        call omp_set_num_threads(threads)
        start=omp_get_wtime()
        call assemble_parametric_bem_p2(sphere,sphere%faces,options,op, &
            trace_dofs,flux_dofs,report,status)
        elapsed=omp_get_wtime()-start
        if (status/=0) error stop 'P2 scaling assembly'
        write(*,'(a,2(1x,i12),4(1x,es17.8))') 'BEM3D P2 scaling threads evaluations seconds V W ', &
            threads,report%kernel_evaluations,elapsed,sum(op%v),sum(abs(op%w)), &
            report%error_ratio
    end do
    call check_summary('Curved P2 BEM harmonics, convergence and scaling')
    if (.not. passed) error stop 1
contains
    pure real(dp) function harmonic(l,z) result(value)
        integer, intent(in) :: l
        real(dp), intent(in) :: z
        value=z
        if (l==2) value=0.5_dp*(3*z*z-1)
    end function
    subroutine record(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed=passed.and.ok
        call check_condition(ok,label)
    end subroutine
end program
