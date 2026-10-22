program test_lab_tensor_solver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                     tensor_load, tensor_l2_error, apply_tensor_operator
    use fortfem_tensor_solver, only: tensor_solve_info_t, solve_tensor_diffusion, &
                    tensor_multigrid_t, build_tensor_multigrid, free_tensor_multigrid, &
                                     tensor_transfer, tensor_vcycle
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op
    type(tensor_solve_info_t) :: info
    type(tensor_multigrid_t) :: mg
    real(dp), allocatable :: rhs(:), x(:), check(:), coarse(:), fine(:), back(:)
    real(dp) :: error, ratio
    integer :: p, h, status, i, anisotropy
    real(dp) :: chi = 1, angle = 0, k12
    print *, 'h p dofs iterations residual L2 (criterion 4)'
    do h = 2, 4, 2
        do p = 2, 8
            call structured_tensor_mesh(2, [h, h, 1], p, mesh, status)
            if (status /= 0) error stop 'mesh'
            call build_tensor_operator(mesh, coefficient, op, status, nq=p + 3)
            if (status /= 0) error stop 'operator'
            call tensor_load(op, forcing, rhs)
            allocate (x(mesh%ndof)); x = 0
            call solve_tensor_diffusion(mesh, coefficient, rhs, 1.e-10_dp, x, info)
            error = tensor_l2_error(op, x, exact)
           print '(4i7,2es13.4)', h, p, mesh%ndof, info%iterations, info%residual, error
            if (info%status /= 0) then
                print *, 'p-MG setup/solve status:', info%status
                error stop 'p-MG PCG'
            end if
            if (info%iterations > 40) error stop 'p-MG iteration growth'
            if (p == 8 .and. error > 1.e-7_dp) error stop 'smooth MMS error'
            deallocate (x, rhs)
        end do
    end do
    do anisotropy = 0, 6, 3
        chi = 10.0_dp**anisotropy
        call structured_tensor_mesh(2, [4, 4, 1], 4, mesh, status)
        call build_tensor_operator(mesh, coefficient, op, status)
        call tensor_load(op, forcing, rhs)
        allocate (x(mesh%ndof), check(mesh%ndof)); x = 0
        call solve_tensor_diffusion(mesh, coefficient, rhs, 1.e-10_dp, x, info, &
                                    maxiter=800)
        call apply_tensor_operator(op, x, check)
        ratio = norm2(rhs - check)/norm2(rhs)
        print '(a,es10.2,i6,2es13.4)', 'aligned chi, it, residual, L2: ', &
            chi, info%iterations, ratio, tensor_l2_error(op, x, exact)
        if (info%status /= 0 .or. ratio > 2.e-10_dp) error stop 'anisotropic solve'
        deallocate (x, rhs, check)
    end do
    angle = acos(-1.0_dp)/4
    call structured_tensor_mesh(2, [4, 4, 1], 4, mesh, status)
    call build_tensor_operator(mesh, coefficient, op, status)
    call tensor_load(op, forcing, rhs)
    allocate (x(mesh%ndof), check(mesh%ndof)); x = 0
    call solve_tensor_diffusion(mesh, coefficient, rhs, 1.e-10_dp, x, info, maxiter=800)
    call apply_tensor_operator(op, x, check)
    ratio = norm2(rhs - check)/norm2(rhs)
    print '(a,i6,2es13.4)', 'unaligned chi=1e6 it, residual, L2: ', &
        info%iterations, ratio, tensor_l2_error(op, x, exact)
    ! Unaligned degradation is reported, rather than hidden by tuning.
    if (info%status /= 0 .or. ratio > 2.e-10_dp) error stop 'unaligned solve'
    deallocate (x, rhs, check)
    chi = 1; angle = 0
    call build_tensor_multigrid(mesh, coefficient, mg, status)
    if (status /= 0) error stop 'transfer hierarchy'
    allocate (fine(mesh%ndof), x(mesh%ndof), check(mesh%ndof), rhs(mesh%ndof), &
              coarse(mg%level(2)%op%mesh%ndof), back(mg%level(2)%op%mesh%ndof))
    do i = 1, size(fine)
        fine(i) = sin(real(i, dp))
    end do
    do i = 1, size(coarse)
        coarse(i) = cos(real(i, dp))
    end do
    where (mesh%boundary) fine = 0
    where (mg%level(2)%op%mesh%boundary) coarse = 0
    call tensor_transfer(mg, 1, coarse, x, .false.)
    call tensor_transfer(mg, 1, back, fine, .true.)
    ratio = abs(dot_product(fine, x) - dot_product(back, coarse))
    if (ratio > 1.e-12_dp*max(1.0_dp, norm2(fine)*norm2(coarse))) &
        error stop 'global transfer adjoint identity'
    call tensor_vcycle(mg, 1, fine, check, status)
    call tensor_vcycle(mg, 1, x, rhs, status)
    ratio = abs(dot_product(fine, rhs) - dot_product(check, x))
    if (ratio > 1.e-11_dp*max(1.0_dp, norm2(fine)*norm2(x))) &
        error stop 'V-cycle symmetry'
    if (dot_product(fine, check) <= 0) error stop 'V-cycle positive energy'
    call free_tensor_multigrid(mg)

contains
    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        if (size(x) /= 2) error stop 'coefficient dimension'
        k = 0
        k(1, 1) = chi*cos(angle)**2 + sin(angle)**2
        k(2, 2) = chi*sin(angle)**2 + cos(angle)**2
        k(1, 2) = (chi - 1)*cos(angle)*sin(angle)
        k(2, 1) = k(1, 2)
        c = 0
    end subroutine
    function exact(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sin(acos(-1.0_dp)*x(1))*sin(acos(-1.0_dp)*x(2))
    end function
    function forcing(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        k12 = (chi - 1)*cos(angle)*sin(angle)
        value = (chi + 1)*acos(-1.0_dp)**2*exact(x) - &
                2*k12*acos(-1.0_dp)**2*cos(acos(-1.0_dp)*x(1))* &
                cos(acos(-1.0_dp)*x(2))
    end function
end program
