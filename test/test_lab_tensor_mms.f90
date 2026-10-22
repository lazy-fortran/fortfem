program test_lab_tensor_mms
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t,structured_tensor_mesh, &
        ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t,build_tensor_operator, &
        tensor_load,tensor_l2_error
    use fortfem_tensor_solver, only: tensor_solve_info_t,solve_tensor_diffusion
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op
    type(tensor_solve_info_t) :: info
    real(dp), allocatable :: rhs(:),x(:),bc(:)
    real(dp) :: error,previous,rate,axis_error,radius
    integer :: domain,p,h,level,status,i,dim
    print *, 'domain dim h p dofs iterations L2 rate axis-node-max (criteria 2,4)'
    do domain=1,3
        dim=2
        if(domain==3) dim=3
        do p=2,3
            previous=0
            do level=0,2
                h=2**level
                if(domain==1) then
                    call structured_tensor_mesh(2,[h,8*h,1],p,mesh,status, &
                        map=annulus,periodic=[.false.,.true.,.false.],q=p+2)
                else
                    call ogrid_tensor_mesh(dim,h,p,mesh,status,nphi=8*h,q=p+2)
                end if
                if(status/=0) error stop 'MMS mesh'
                call build_tensor_operator(mesh,coefficient,op,status,nq=p+3)
                if(status/=0) error stop 'MMS geometry'
                call tensor_load(op,forcing,rhs)
                allocate(x(mesh%ndof),bc(mesh%ndof))
                do i=1,mesh%ndof
                    bc(i)=exact(mesh%coordinates(:,i))
                end do
                x=0
                call solve_tensor_diffusion(mesh,coefficient,rhs,1.e-11_dp,x, &
                    info,bc=bc,maxiter=400)
                if(info%status/=0) error stop 'MMS PCG'
                error=tensor_l2_error(op,x,exact)
                rate=0
                if(level>0) rate=log(previous/error)/log(2.0_dp)
                axis_error=0
                do i=1,mesh%ndof
                    if(dim==2) then
                        radius=norm2(mesh%coordinates(:,i))
                    else
                        radius=norm2([norm2(mesh%coordinates(1:2,i))-3, &
                            mesh%coordinates(3,i)])
                    end if
                    if(radius<0.2_dp) axis_error=max(axis_error,abs(x(i)-bc(i)))
                end do
                print '(6i6,es12.3,f7.2,es12.3)',domain,dim,h,p,mesh%ndof, &
                    info%iterations,error,rate,axis_error
                if(level==2 .and. rate<real(p,dp)+0.5_dp) &
                    error stop 'h-convergence rate loss'
                previous=error
                deallocate(x,bc,rhs)
            end do
        end do
        previous=huge(1.0_dp)
        do p=2,6
            if(domain==1) then
                call structured_tensor_mesh(2,[2,16,1],p,mesh,status, &
                    map=annulus,periodic=[.false.,.true.,.false.],q=p+2)
            else
                call ogrid_tensor_mesh(dim,2,p,mesh,status,nphi=16,q=p+2)
            end if
            if(status/=0) error stop 'p-refinement mesh'
            call build_tensor_operator(mesh,coefficient,op,status,nq=p+3)
            if(status/=0) error stop 'p-refinement geometry'
            call tensor_load(op,forcing,rhs)
            allocate(x(mesh%ndof),bc(mesh%ndof))
            do i=1,mesh%ndof
                bc(i)=exact(mesh%coordinates(:,i))
            end do
            x=0
            call solve_tensor_diffusion(mesh,coefficient,rhs,1.e-11_dp,x, &
                info,bc=bc,maxiter=400)
            if(info%status/=0) error stop 'p-refinement solve'
            error=tensor_l2_error(op,x,exact)
            print '(a,2i5,es13.4)', 'p convergence: ',domain,p,error
            if(error>0.7_dp*previous) error stop 'p-convergence loss'
            previous=error
            deallocate(x,bc,rhs)
        end do
    end do
contains
    subroutine annulus(xi,x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: r,theta
        r=0.3_dp+0.7_dp*xi(1)
        theta=2*acos(-1.0_dp)*xi(2)
        x=[r*cos(theta),r*sin(theta)]
    end subroutine
    subroutine coefficient(x,k,c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:,:),c
        integer :: i
        k=0
        do i=1,size(x)
            k(i,i)=1
        end do
        c=0
    end subroutine
    function exact(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value=exp(dot_product([0.1_dp,0.2_dp,0.3_dp], &
            [x(1),x(2),merge(x(size(x)),0.0_dp,size(x)==3)]))
    end function
    function forcing(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value=-0.05_dp*exact(x)
        if(size(x)==3) value=-0.14_dp*exact(x)
    end function
end program
