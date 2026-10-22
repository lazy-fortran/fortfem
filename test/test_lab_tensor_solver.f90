program test_lab_tensor_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t,build_tensor_operator, &
        tensor_load,tensor_l2_error,apply_tensor_operator
    use fortfem_tensor_solver, only: tensor_solve_info_t,solve_tensor_diffusion
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op
    type(tensor_solve_info_t) :: info
    real(dp), allocatable :: rhs(:),x(:),check(:)
    real(dp) :: error,ratio
    integer :: p,h,status,i,anisotropy
    real(dp) :: chi=1
    print *, 'h p dofs iterations residual L2 (criterion 4)'
    do h=2,4,2
        do p=2,8
            call structured_tensor_mesh(2,[h,h,1],p,mesh,status)
            if(status/=0) error stop 'mesh'
            call build_tensor_operator(mesh,coefficient,op,status,nq=p+3)
            if(status/=0) error stop 'operator'
            call tensor_load(op,forcing,rhs)
            allocate(x(mesh%ndof)); x=0
            call solve_tensor_diffusion(mesh,coefficient,rhs,1.e-10_dp,x,info)
            error=tensor_l2_error(op,x,exact)
            print '(4i7,2es13.4)',h,p,mesh%ndof,info%iterations,info%residual,error
            if(info%status/=0) error stop 'p-MG PCG'
            if(info%iterations>40) error stop 'p-MG iteration growth'
            if(p==8 .and. error>1.e-7_dp) error stop 'smooth MMS error'
            deallocate(x,rhs)
        end do
    end do
    do anisotropy=0,6,3
        chi=10.0_dp**anisotropy
        call structured_tensor_mesh(2,[4,4,1],4,mesh,status)
        call build_tensor_operator(mesh,coefficient,op,status)
        call tensor_load(op,forcing,rhs)
        allocate(x(mesh%ndof),check(mesh%ndof)); x=0
        call solve_tensor_diffusion(mesh,coefficient,rhs,1.e-10_dp,x,info, &
            maxiter=800)
        call apply_tensor_operator(op,x,check)
        ratio=norm2(rhs-check)/norm2(rhs)
        print '(a,es10.2,i6,2es13.4)','aligned chi, it, residual, L2: ', &
            chi,info%iterations,ratio,tensor_l2_error(op,x,exact)
        if(info%status/=0 .or. ratio>2.e-10_dp) error stop 'anisotropic solve'
        deallocate(x,rhs,check)
    end do
contains
    subroutine coefficient(x,k,c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:,:),c
        k=0; k(1,1)=chi; k(2,2)=1; c=0
    end subroutine
    function exact(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value=sin(acos(-1.0_dp)*x(1))*sin(acos(-1.0_dp)*x(2))
    end function
    function forcing(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value=(chi+1)*acos(-1.0_dp)**2*exact(x)
    end function
end program
