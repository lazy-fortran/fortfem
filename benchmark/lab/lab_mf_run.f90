program lab_mf_run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_lab_scalar, only: lab_scalar
    use fortfem_generated_lab_transport, only: lab_transport
    use fortfem_generated_lab_geometry, only: lab_geometry
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh, &
                                   ogrid_tensor_mesh, build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, &
        build_tensor_operator, tensor_load, &
                                       assemble_tensor_operator, apply_tensor_operator
    use fortfem_tensor_solver, only: tensor_solve_info_t, solve_tensor_diffusion
    use fortfem_field_aligned_constitutive_tensor, only: &
        evaluate_field_aligned_constitutive_tensor
    use fortsparse, only: fortsparse_status_t, csc_t, csc_from_triplet
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_advanced_solvers, only: solve_sparse, solver_options_t, solver_stats_t
    use fortfem_tensor_contractions, only: tensor_contract
    use fortnum_linalg, only: inv3
    implicit none
    type(tensor_mesh_t) :: mesh, original
    type(tensor_operator_t) :: op, error_op
    type(tensor_solve_info_t) :: info
    type(csr_matrix_t) :: csr
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: sparse_status
    type(solver_options_t) :: options
    type(solver_stats_t) :: stats
    integer, allocatable :: rows(:), cols(:)
    real(dp), allocatable :: lift(:), action(:), mapped(:, :, :)
    real(dp), allocatable :: rhs(:), solution(:), bc(:), values(:), local(:)
    real(dp), allocatable :: scratch1(:), scratch2(:), b(:, :), d(:, :), jets(:, :)
    real(dp) :: ratio, err, norm, herr, hnorm, aerr, anorm, projection, value
    real(dp) :: ballerr, ballnorm, g(12), radius, phi
   real(dp) :: exact, gradient(3), eg(3), source, direction(3), x(3), t(9), inv(3, 3), w
    integer :: p, n, status, i, e, q, ws, dir, node, fine_nnz
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
    if (preconditioner /= 'pmg' .and. preconditioner /= 'ichol_row') &
        error stop 'MF adapter requires pmg or polar ichol_row'
    if (preconditioner == 'pmg') then
        call ogrid_tensor_mesh(3, n, p, original, status, nphi=20*n, q=p + 2)
        if (status /= 0) error stop 'lab O-grid failed'
        mapped = original%geometry
        do e = 1, original%nelem
            do node = 1, size(mapped, 2)
                x = original%geometry(:, node, e)
                radius = norm2(x(:2))
                phi = atan2(x(2), x(1))
                call lab_geometry(radius - 3, x(3), phi, g)
                mapped(:, node, e) = [g(1)*cos(phi), g(1)*sin(phi), g(2)]
            end do
        end do
        call build_tensor_mesh(3, p, p + 2, original%vertices, mapped, mesh, status, &
                               gather=original%dof)
        if (status == 0) mesh%boundary = original%boundary
    else
        if (problem /= 'L5') error stop 'polar comparison requires L5'
        call structured_tensor_mesh(3, [n, 8*n, 20*n], p, mesh, status, &
                        map=polar_geometry, periodic=[.false., .true., .true.], q=p + 2)
        if (status == 0) call collapse_polar_axis()
    end if
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
    fine_nnz = 0
    if (preconditioner == 'pmg') then
        call solve_tensor_diffusion(mesh, coefficient, rhs, 1.e-10_dp, solution, info, &
                                    bc=bc, maxiter=1500)
    else
        ! The polar comparison uses the same tensor operator's assembled
        ! oracle, avoiding a multigrid topology assumption at collapsed cells.
        call assemble_tensor_operator(op, csr)
        allocate (rows(size(csr%val)), cols(size(csr%val)))
        do i = 1, csr%n
            rows(csr%rowptr(i):csr%rowptr(i + 1) - 1) = i
        end do
        cols = csr%colidx
        call csc_from_triplet(csr%n, csr%n, rows, cols, csr%val, matrix, sparse_status)
        if (sparse_status%code /= 0) error stop 'polar CSC failed'
        fine_nnz = matrix%nnz
        allocate (lift(mesh%ndof), action(mesh%ndof))
        lift = 0
        where (mesh%boundary) lift = bc
        call apply_tensor_operator(op, lift, action, eliminate=.false.)
        rhs = rhs - action
        where (mesh%boundary) rhs = bc
        options%method = 'pcg'
        options%preconditioner = 'ichol_row'
        options%fill_level = 64
        options%tolerance = 1.e-10_dp
        options%max_iterations = 20000
        call solve_sparse(matrix, rhs, solution, options, stats)
        info%status = 1
        if (stats%converged) info%status = 0
        info%iterations = stats%iterations
        call apply_tensor_operator(op, solution, action, eliminate=.false.)
        action = action - forcing_vector(op)
        where (mesh%boundary) action = solution - bc
        info%residual = norm2(action)/max(norm2(rhs), tiny(1._dp))
    end if
    if (info%status /= 0 .or. info%residual > 2.e-8_dp) then
        print *, info%status, info%iterations, info%residual
        error stop 'lab MF solve failed'
    end if
    call build_tensor_operator(mesh, coefficient, error_op, status, nq=p + 3)
    if (status /= 0) error stop 'lab MF error geometry failed'
    b = error_op%line(1)%b(:, :, 1)
    d = error_op%line(1)%d(:, :, 1)
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
    ballerr = 0
    ballnorm = 0
    !$omp parallel do default(shared) &
    !$omp private(e,local,values,jets,scratch1,scratch2,dir,q,x,exact,eg,source, &
    !$omp direction,inv,status,gradient,w,t) &
    !$omp reduction(+:err,norm,herr,hnorm,aerr,anorm,projection) &
    !$omp reduction(+:ballerr,ballnorm) schedule(static)
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
            if (sum((x - [3.1_dp, 0._dp, 0._dp])**2) < .04_dp) then
                ballerr = ballerr + w*(values(q) - exact)**2
                ballnorm = ballnorm + w*exact**2
            end if
            if (t(9) < .04_dp) then
                aerr = aerr + w*(values(q) - exact)**2
                anorm = anorm + w*exact**2
            end if
        end do
    end do
    !$omp end parallel do
    value = 0
    if (problem == 'L3') value = projection/norm
    write (*, '(a,1x,i0,1x,6(es24.16,1x),i0,1x,es24.16)') &
        'LAB', mesh%ndof, sqrt(err/norm), &
        sqrt(herr/hnorm), sqrt(aerr/max(anorm, tiny(1._dp))), value, info%residual, &
        real(fine_nnz,dp), info%iterations, sqrt(ballerr/max(ballnorm, tiny(1._dp)))
contains
    function forcing_vector(operator) result(load)
        type(tensor_operator_t), intent(in) :: operator
        real(dp), allocatable :: load(:)
        call tensor_load(operator, forcing, load)
    end function
    subroutine polar_geometry(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: angle, phi, g(12)
        angle = 2*acos(-1._dp)*xi(2)
        phi = 2*acos(-1._dp)*xi(3)
        call lab_geometry(xi(1)*cos(angle), -xi(1)*sin(angle), phi, g)
        x = [g(1)*cos(phi), g(1)*sin(phi), g(2)]
    end subroutine
    subroutine collapse_polar_axis()
        integer, allocatable :: mapping(:), representative(:)
        real(dp), allocatable :: coordinates(:, :)
        logical, allocatable :: boundary(:)
        integer :: nr, nt, np, ir, it, ip, old, nd, axis
        nr = n*p + 1
        nt = 8*n*p
        np = 20*n*p
        allocate (mapping(mesh%ndof), representative(mesh%ndof))
        nd = 0
        do ip = 0, np - 1
            nd = nd + 1
            axis = nd
            representative(nd) = 1 + nr*nt*ip
            do it = 0, nt - 1
                mapping(1 + nr*(it + nt*ip)) = axis
                do ir = 1, nr - 1
                    old = 1 + ir + nr*(it + nt*ip)
                    nd = nd + 1
                    mapping(old) = nd
                    representative(nd) = old
                end do
            end do
        end do
        allocate (coordinates(3, nd), boundary(nd))
        coordinates = mesh%coordinates(:, representative(:nd))
        boundary = mesh%boundary(representative(:nd))
        do ip = 0, np - 1
            boundary(mapping(1 + nr*nt*ip)) = .false.
        end do
       mesh%dof = reshape(mapping(reshape(mesh%dof, [size(mesh%dof)])), shape(mesh%dof))
        call move_alloc(coordinates, mesh%coordinates)
        call move_alloc(boundary, mesh%boundary)
        mesh%ndof = nd
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
