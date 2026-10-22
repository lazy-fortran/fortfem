module bem_scaling_geometry
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    use fortfem_generated_bem_synthetic_surface, only: bem_synthetic_surface
    use fortfem_generated_bem_kernel, only: bem_kernel
    implicit none
    type, extends(bem_surface_t) :: synthetic_t
        real(dp), allocatable :: parameters(:,:)
        integer, allocatable :: faces(:,:)
    contains
        procedure :: evaluate
    end type
contains
    pure subroutine evaluate(self,panel,uv,x,a,b,status)
        class(synthetic_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3),a(3),b(3)
        integer, intent(out) :: status
        real(dp) :: local(2,3),angles(2),dt(3),df(3),pi
        integer :: i
        pi=acos(-1.0_dp)
        local=self%parameters(:,self%faces(:,panel))
        do i=2,3
            local(:,i)=local(:,1)+modulo(local(:,i)-local(:,1)+pi,2*pi)-pi
        end do
        angles=local(:,1)+uv(1)*(local(:,2)-local(:,1))+ &
            uv(2)*(local(:,3)-local(:,1))
        call bem_synthetic_surface(angles(1),angles(2),5.5_dp,.55_dp,.1_dp, &
            5.0_dp,x(1),x(2),x(3),dt(1),dt(2),dt(3),df(1),df(2),df(3))
        a=dt*(local(1,2)-local(1,1))+df*(local(2,2)-local(2,1))
        b=dt*(local(1,3)-local(1,1))+df*(local(2,3)-local(2,1))
        status=0
    end subroutine

    pure subroutine exact(x,normal,u,q)
        real(dp), intent(in) :: x(3),normal(3)
        real(dp), intent(out) :: u,q
        real(dp) :: phi,source(3),d(3),green,double
        integer :: i
        ! Five interior monopoles. Their potential decays in the exterior.
        u=0; q=0
        do i=0,4
            phi=2*acos(-1.0_dp)*i/5
            source=[5.5_dp*cos(phi),5.5_dp*sin(phi),0.0_dp]
            d=x-source
            call bem_kernel(d(1),d(2),d(3),-normal(1),-normal(2),-normal(3), &
                green,double)
            u=u+green; q=q+double
        end do
    end subroutine
end module

program fortfem_parametric_bem_scaling
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use bem_scaling_geometry, only: synthetic_t, exact
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_parametric_bem_high_order, only: assemble_parametric_bem_p2, &
        assemble_parametric_bem_p2_circulant, bem_pair_options_t, bem_pair_report_t
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t
    use fortfem_bem_circulant, only: bem_circulant_t
    use fortnum_linalg, only: dense_solve
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    type(synthetic_t) :: surface
    type(bem_operators_t) :: op
    type(bem_circulant_t) :: circular
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    integer, allocatable :: td(:,:),fd(:,:)
    real(dp), allocatable :: vertices(:,:),u(:),rhs(:),q(:),tmp(:)
    real(dp), allocatable :: xi(:),eta(:),weights(:),projected(:,:)
    real(dp) :: t0,assembly_seconds,solve_seconds,error,best_error,tolerance
    real(dp) :: x(3),a(3),b(3),jac,value,flux,phi(3)
    real(dp) :: mass(3,3),load(3),solution(3),denominator,numerator
    real(dp), parameter :: nodes(2,6)=reshape( &
        [0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp, &
        .5_dp,0.0_dp,.5_dp,.5_dp,0.0_dp,.5_dp],[2,6])
    integer :: ntheta,nphi,status,p,j,t,nf,unit,ios,audit_degree
    integer(int64) :: stored_bytes
    character(64) :: arg,mode
    character(256) :: line
    ntheta=4; nphi=40; mode='circulant'; tolerance=1.e-4_dp; audit_degree=24
    if (command_argument_count()>=1) then
        call get_command_argument(1,arg); read(arg,*) ntheta
    end if
    if (command_argument_count()>=2) then
        call get_command_argument(2,arg); read(arg,*) nphi
    end if
    if (command_argument_count()>=3) call get_command_argument(3,mode)
    if (command_argument_count()>=4) then
        call get_command_argument(4,arg); read(arg,*) tolerance
    end if
    if (command_argument_count()>=5) then
        call get_command_argument(5,arg); read(arg,*) audit_degree
    end if
    if (ntheta<4 .or. nphi<10 .or. mod(nphi,5)/=0) error stop 'mesh counts'
    if (trim(mode)/='dense' .and. trim(mode)/='circulant' .and. &
        trim(mode)/='projection') error stop 'mode'
    call generate_torus_surface_mesh(5.5_dp,.55_dp,ntheta,nphi,vertices, &
        surface%faces,surface%parameters)
    nf=size(surface%faces,2)
    call triangle_duffy_quadrature(audit_degree,xi,eta,weights,status)
    if (status/=0) error stop 'audit rule'
    allocate(projected(3,nf))
    !$omp parallel do default(shared) private(p,t,mass,load,solution,status,phi,jac,flux)
    do p=1,nf
        mass=0; load=0
        do t=1,size(weights)
            call geometry(surface,p,[xi(t),eta(t)],phi,jac,flux)
            mass=mass+weights(t)*jac*spread(phi,2,3)*spread(phi,1,3)
            load=load+weights(t)*jac*phi*flux
        end do
        call dense_solve(mass,load,solution,status)
        if (status/=0) error stop 'projection mass'
        projected(:,p)=solution
    end do
    !$omp end parallel do
    denominator=0; numerator=0
    !$omp parallel do default(shared) private(p,t,phi,jac,flux) &
    !$omp reduction(+:denominator,numerator)
    do p=1,nf
        do t=1,size(weights)
            call geometry(surface,p,[xi(t),eta(t)],phi,jac,flux)
            denominator=denominator+weights(t)*jac*flux**2
            numerator=numerator+weights(t)*jac*(dot_product(phi,projected(:,p))-flux)**2
        end do
    end do
    !$omp end parallel do
    best_error=sqrt(numerator/denominator)
    write(*,'(a,3i10,es18.9)') 'BEM_COST projection ntheta nphi panels rel_L2 ', &
        ntheta,nphi,nf,best_error
    if (trim(mode)=='projection') stop
    options%relative_tolerance=tolerance
    options%absolute_tolerance=1.e-11_dp
    options%maximum_degree=64
    t0=omp_get_wtime()
    if (trim(mode)=='dense') then
        call assemble_parametric_bem_p2(surface,surface%faces,options,op, &
            td,fd,report,status)
    else
        call assemble_parametric_bem_p2_circulant(surface,surface%faces,5,options, &
            circular,td,fd,report,status)
    end if
    assembly_seconds=omp_get_wtime()-t0
    write(*,'(a,3i14,2es18.9)') 'BEM_COST assembly status failed degree ratio seconds ', &
        status,report%failed_pairs,report%maximum_degree_used,report%error_ratio, &
        assembly_seconds
    if (status/=0) error stop 'unresolved assembly'
    allocate(u(maxval(td)),q(3*nf),rhs(3*nf),tmp(3*nf))
    do p=1,nf
        do j=1,6
            call surface%evaluate(p,nodes(:,j),x,a,b,status)
            call exact(x,[0.0_dp,0.0_dp,0.0_dp],value,flux)
            u(td(j,p))=value
        end do
    end do
    t0=omp_get_wtime()
    if (trim(mode)=='dense') then
        rhs=matmul(op%k-.5_dp*op%mass,u)
        call dense_solve(op%v,rhs,q,status)
        stored_bytes=8_int64*(size(op%v,kind=int64)+size(op%k,kind=int64)+ &
            size(op%w,kind=int64)+size(op%mass,kind=int64)+ &
            size(op%adjoint,kind=int64))
    else
        call circular%k%apply(u,rhs,status)
        if (status/=0) error stop 'K apply'
        call circular%mass%apply(u,tmp,status)
        if (status/=0) error stop 'M apply'
        rhs=rhs-.5_dp*tmp
        call circular%v%solve(rhs,q,status)
        stored_bytes=8_int64*(size(circular%v%blocks,kind=int64)+ &
            size(circular%k%blocks,kind=int64)+size(circular%w%blocks,kind=int64)+ &
            size(circular%mass%blocks,kind=int64))
    end if
    solve_seconds=omp_get_wtime()-t0
    if (status/=0) error stop 'flux solve'
    numerator=0
    !$omp parallel do default(shared) private(p,t,phi,jac,flux) reduction(+:numerator)
    do p=1,nf
        do t=1,size(weights)
            call geometry(surface,p,[xi(t),eta(t)],phi,jac,flux)
            numerator=numerator+weights(t)*jac*(dot_product(phi,q(fd(:,p)))-flux)**2
        end do
    end do
    !$omp end parallel do
    error=sqrt(numerator/denominator)
    if (error<best_error-1.e-10_dp) error stop 'inconsistent projection oracle'
    write(*,'(a,a,3i14,4es18.9)') 'BEM_COST result ',trim(mode),nf, &
        omp_get_max_threads(),stored_bytes,assembly_seconds,solve_seconds,error,best_error
    write(*,'(a,2i18,i8)') 'BEM_COST work evaluations accepted degree ', &
        report%kernel_evaluations,report%accepted_pairs,report%maximum_degree_used
    open(newunit=unit,file='/proc/self/status',status='old',iostat=ios)
    if (ios==0) then
        do
            read(unit,'(a)',iostat=ios) line
            if (ios/=0) exit
            if (index(line,'VmHWM:')==1) write(*,'(a)') trim(line)
        end do
        close(unit)
    end if
contains
    subroutine geometry(surface,panel,uv,phi,jac,flux)
        type(synthetic_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: phi(3),jac,flux
        real(dp) :: x(3),a(3),b(3),normal(3),value,g(2),h(2,2)
        integer :: i,status
        call surface%evaluate(panel,uv,x,a,b,status)
        if (status/=0) error stop 'audit geometry'
        call bem_frame(a(1),a(2),a(3),b(1),b(2),b(3),jac, &
            normal(1),normal(2),normal(3))
        call exact(x,normal,value,flux)
        do i=1,3
            call generated_p1_jet(i,uv(1),uv(2),phi(i),g,h)
        end do
    end subroutine
end program
