program lab_mf_run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_lab_scalar, only: lab_scalar
    use fortfem_generated_lab_transport, only: lab_transport
    use fortfem_generated_lab_geometry, only: lab_geometry
    use fortfem_generated_lab_square_disk, only: lab_square_disk
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh
use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, tensor_load
    use fortfem_tensor_solver, only: tensor_solve_info_t, solve_tensor_diffusion
    use fortfem_field_aligned_constitutive_tensor, only: &
        evaluate_field_aligned_constitutive_tensor
    use fortsparse, only: fortsparse_status_t
    use fortfem_tensor_contractions, only: tensor_contract
    use fortnum_linalg, only: inv3
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op, error_op
    type(tensor_solve_info_t) :: info
    real(dp), allocatable :: rhs(:), solution(:), bc(:), values(:), local(:)
    real(dp), allocatable :: scratch1(:), scratch2(:), b(:, :), d(:, :), jets(:, :)
    real(dp) :: ratio, err, norm, herr, hnorm, aerr, anorm, projection, value
   real(dp) :: exact, gradient(3), eg(3), source, direction(3), x(3), t(9), inv(3, 3), w
    integer :: p, n, status, i, e, q, ws, dir
    character(32) :: problem, preconditioner
    character(4096) :: arg
    call get_command_argument(1, problem)
    call get_command_argument(2, arg)
    read (arg, *) p
    call get_command_argument(3, arg)
    read (arg, *) n
    call get_command_argument(4, preconditioner)
    call get_command_argument(5, arg)
    read (arg, *) ratio
    if (problem /= 'L1' .and. problem /= 'L3' .and. problem /= 'L5') &
        error stop 'MF adapter supports L1/L3/L5'
    if (preconditioner /= 'pmg') error stop 'MF adapter requires pmg'
    call structured_tensor_mesh(3, [2*n, 2*n, 20*n], p, mesh, status, &
                             map=geometry, periodic=[.false., .false., .true.], q=p + 2)
    if (status /= 0) error stop 'lab MF mesh failed'
    ! Match the load quadrature to the operator used internally by MF1's solver.
    call build_tensor_operator(mesh, coefficient, op, status)
    if (status /= 0) error stop 'lab MF operator failed'
    call tensor_load(op, forcing, rhs)
    allocate (solution(mesh%ndof), bc(mesh%ndof))
    do i = 1, mesh%ndof
        call physical_data(mesh%coordinates(:, i), exact, eg, source, direction)
        bc(i) = exact
    end do
    solution = 0
    call solve_tensor_diffusion(mesh, coefficient, rhs, 1.e-10_dp, solution, info, &
                                bc=bc, maxiter=1500)
    if (info%status /= 0) then
        print *, info%status, info%iterations, info%residual
        error stop 'lab MF solve failed'
    end if
    call build_tensor_operator(mesh, coefficient, error_op, status, nq=p + 3)
    if (status /= 0) error stop 'lab MF error geometry failed'
    b = error_op%line%b(:, :, 1)
    d = error_op%line%d(:, :, 1)
    ws = maxval([error_op%ni, error_op%no])**3
    allocate (values(product(error_op%no)), jets(3, product(error_op%no)))
    allocate (local(product(error_op%ni)), scratch1(ws), scratch2(ws))
    err = 0
    norm = 0
    herr = 0
    hnorm = 0
    aerr = 0
    anorm = 0
    projection = 0
    !$omp parallel do default(shared) &
    !$omp private(e,local,values,jets,scratch1,scratch2,dir,q,x,exact,eg,source, &
    !$omp direction,inv,status,gradient,w,t) &
    !$omp reduction(+:err,norm,herr,hnorm,aerr,anorm,projection) schedule(static)
    do e = 1, mesh%nelem
        local = solution(mesh%dof(:, e))
      call tensor_contract(error_op%ni,error_op%no,b,b,b,local,values,scratch1,scratch2)
        do dir = 1, 3
            select case (dir)
            case (1)
                call tensor_contract(error_op%ni, error_op%no, d, b, b, local, &
                                     jets(dir, :), scratch1, scratch2)
            case (2)
                call tensor_contract(error_op%ni, error_op%no, b, d, b, local, &
                                     jets(dir, :), scratch1, scratch2)
            case (3)
                call tensor_contract(error_op%ni, error_op%no, b, b, d, local, &
                                     jets(dir, :), scratch1, scratch2)
            end select
        end do
        do q = 1, size(values)
            x = error_op%points(:, q, e)
            call physical_data(x, exact, eg, source, direction)
            call inv3(error_op%jac(:, :, q, e), inv, status)
            if (status /= 0) error stop 'lab MF error inverse failed'
            gradient = matmul(transpose(inv), jets(:, q))
            w = error_op%volume(q, e)
            err = err + w*(values(q) - exact)**2
            norm = norm + w*exact**2
            herr = herr + w*sum((gradient - eg)**2)
            hnorm = hnorm + w*sum(eg**2)
            projection = projection + w*values(q)*exact
            call lab_transport(sqrt(sum(x(:2)**2)), atan2(x(2), x(1)), x(3), t)
            if (t(9) < .04_dp) then
                aerr = aerr + w*(values(q) - exact)**2
                anorm = anorm + w*exact**2
            end if
        end do
    end do
    !$omp end parallel do
    value = 0
    if (problem == 'L3') value = projection/norm
    write (*, '(a,1x,i0,1x,6(es24.16,1x),i0)') 'LAB', mesh%ndof, sqrt(err/norm), &
        sqrt(herr/hnorm), sqrt(aerr/max(anorm, tiny(1._dp))), value, info%residual, &
        0._dp, info%iterations
contains
    subroutine geometry(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: disk(2), g(12), phi
        call lab_square_disk(2*xi(1) - 1, 2*xi(2) - 1, disk)
        phi = 2*acos(-1._dp)*xi(3)
        call lab_geometry(disk(1), -disk(2), phi, g)
        x = [g(1)*cos(phi), g(1)*sin(phi), g(2)]
    end subroutine
    subroutine physical_data(x, u, gradient, forcing, direction)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u, gradient(3), forcing, direction(3)
        real(dp) :: data(9), phi
        direction = 0
        if (problem == 'L3') then
            phi = atan2(x(2), x(1))
            call lab_transport(sqrt(sum(x(:2)**2)), phi, x(3), data)
            u = data(1)
            gradient = [data(2)*cos(phi) - data(3)*sin(phi), &
                        data(2)*sin(phi) + data(3)*cos(phi), data(4)]
            forcing = data(5)
            direction = [data(6)*cos(phi) - data(7)*sin(phi), &
                         data(6)*sin(phi) + data(7)*cos(phi), data(8)]
        else
            call lab_scalar(x(1), x(2), x(3), data(:8))
            u = data(1)
            gradient = data(2:4)
            forcing = data(5)
        end if
    end subroutine
    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        real(dp) :: u, grad(3), f, b(3)
        type(fortsparse_status_t) :: status
        call physical_data(x, u, grad, f, b)
        ! L1 has no field term. Supply any unit vector in its isotropic case.
        if (problem /= 'L3') b = [1._dp, 0._dp, 0._dp]
        call evaluate_field_aligned_constitutive_tensor(ratio, 1._dp, b, k, status)
        if (status%code /= 0) error stop 'lab MF tensor failed'
        c = 0
    end subroutine
    function forcing(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value, u, grad(3), b(3)
        call physical_data(x, u, grad, value, b)
    end function
end program
