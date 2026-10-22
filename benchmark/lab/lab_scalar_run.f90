program lab_scalar_run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_lab_scalar, only: lab_scalar
    use fortfem_generated_lab_transport, only: lab_transport
    use fortfem_generated_lab_guenter, only: lab_guenter
    use fortfem_tetra_lagrange_arbitrary_order, only: tetra_lagrange_t, &
                                                      initialize_tetra_lagrange, evaluate_tetra_lagrange
    use fortfem_tetra_lagrange_global_dof_map, only: &
        build_tetra_lagrange_dof_map
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_advanced_solvers, only: solve_sparse, solver_options_t, &
                                        solver_stats_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, fortsparse_status_t
    use fortnum_linalg, only: det3, inv3
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_feec, only: assemble_tensor_diffusion_matrix
    implicit none
    type(tetra_lagrange_t) :: basis
    type(csc_t) :: matrix
    type(fortsparse_status_t) :: status
    type(native_ldlt_t) :: factor
    type(solver_options_t) :: opts
    type(solver_stats_t) :: stats
    real(dp), allocatable :: vertices(:, :), rhs(:), solution(:), bcvalue(:)
    real(dp), allocatable :: val(:, :), grad(:, :, :), qw(:), qx(:), qy(:), qz(:)
    real(dp), allocatable :: entries(:), local_rhs(:, :), gx(:), gw(:)
    integer, allocatable :: cells(:, :), faces(:, :), maps(:, :), rows(:), cols(:)
    logical, allocatable :: constrained(:)
    real(dp) :: ratio, err, norm, herr, hnorm, aerr, anorm, value, exact
    real(dp) :: gradient(3), eg(3), b(3), source, w, point(3), t(9), center
    real(dp) :: jac(3, 3), inverse(3, 3), det, residual
    real(dp), allocatable :: local_matrix(:, :), physical_grad(:, :), kv(:, :)
    real(dp), allocatable :: qgrad(:, :, :), kt(:, :, :), qweight(:)
    integer :: dim, nv, nc, nd, nl, nq, p, ios, unit, c, i, j, q, d, entry
    integer :: iterations, f, rr, cc, index, ng, k
    character(4096) :: meshfile, arg
    character(32) :: problem, preconditioner
    call get_command_argument(1, problem)
    call get_command_argument(2, arg)
    read (arg, *) p
    call get_command_argument(3, meshfile)
    call get_command_argument(4, preconditioner)
    call get_command_argument(5, arg)
    read (arg, *) ratio
    open (newunit=unit, file=trim(meshfile), status='old')
    read (unit, *) dim, nv, nc
    allocate (vertices(3, nv), cells(dim + 1, nc), faces(dim + 1, nc))
    do i = 1, nv
        read (unit, *) vertices(:, i)
    end do
    do c = 1, nc
        read (unit, *) cells(:, c), faces(:, c)
    end do
    close (unit)
    call initialize_tetra_lagrange(p, basis, ios)
    if (ios /= 0) error stop 'invalid tetrahedral degree'
    if (dim == 3) then
        call build_tetra_lagrange_dof_map(p, cells, maps, nd, ios)
        if (ios /= 0) error stop 'tetrahedral map failed'
        call tetra_duffy_quadrature(2*p + 6, qx, qy, qz, qw, ios)
        if (ios /= 0) error stop 'tetrahedral quadrature failed'
        nl = size(maps, 1)
    else
        if (p /= 1) error stop '2D adapter currently requires P1'
        maps = cells
        nd = nv
        nl = 3
        ng = 6
        allocate (gx(ng), gw(ng), qx(ng**2), qy(ng**2), qz(ng**2), qw(ng**2))
        call gauss_legendre_ab(ng, 0.0_dp, 1.0_dp, gx, gw)
        k = 0
        do j = 1, ng
            do i = 1, ng
                k = k + 1
                qx(k) = gx(i)
                qy(k) = (1 - gx(i))*gx(j)
                qw(k) = gw(i)*gw(j)*(1 - gx(i))
            end do
        end do
        qz = 0
    end if
    nq = size(qw)
    allocate (val(nl, nq), grad(3, nl, nq))
    do q = 1, nq
        if (dim == 3) then
            call evaluate_tetra_lagrange(basis, [qx(q), qy(q), qz(q)], &
                                         val(:, q), grad(:, :, q), ios)
            if (ios /= 0) error stop 'reference basis failed'
        else
            val(:, q) = [1 - qx(q) - qy(q), qx(q), qy(q)]
            grad(:, :, q) = reshape([-1.0_dp, -1.0_dp, 0.0_dp, &
                                     1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [3, 3])
        end if
    end do
    allocate (constrained(nd), bcvalue(nd), rhs(nd), solution(nd))
    constrained = .false.
    bcvalue = 0
    do c = 1, nc
        do i = 1, nl
            do f = 1, dim + 1
                if (faces(f, c) == 0) cycle
                if (dim == 3) then
                    if (basis%barycentric_indices(f, i) /= 0) cycle
                    point = matmul(vertices(:, cells(:, c)), &
                                   real(basis%barycentric_indices(:, i), dp)/p)
                else
                    if (i == f) cycle
                    point = vertices(:, cells(i, c))
                end if
                constrained(maps(i, c)) = .true.
                call problem_data(point, exact, eg, source, b)
                bcvalue(maps(i, c)) = exact
            end do
        end do
    end do
    allocate (rows(nl**2*nc + nd), cols(nl**2*nc + nd))
    allocate (entries(size(rows)), local_rhs(nl, nc))
    !$omp parallel do default(shared) private(c,d,jac,inverse,det,ios, &
    !$omp physical_grad,kv,local_matrix,point,exact,eg,source,b,q,w,i,j, &
    !$omp entry,rr,cc,qgrad,kt,qweight,status) schedule(static)
    do c = 1, nc
        jac = 0
        jac(3, 3) = 1
        do d = 1, dim
            jac(:, d) = vertices(:, cells(d + 1, c)) - vertices(:, cells(1, c))
        end do
        det = det3(jac)
        call inv3(jac, inverse, ios)
        if (ios /= 0 .or. det <= 0) error stop 'invalid lab cell'
        allocate (physical_grad(3, nl), kv(3, nl), local_matrix(nl, nl))
        if (dim == 2) allocate (qgrad(1, nl, 2), kt(1, 2, 2), qweight(1))
        local_matrix = 0
        local_rhs(:, c) = 0
        do q = 1, nq
            point = vertices(:, cells(1, c)) + &
                    matmul(jac, [qx(q), qy(q), qz(q)])
            call problem_data(point, exact, eg, source, b)
            physical_grad = matmul(transpose(inverse), grad(:, :, q))
            w = det*qw(q)
            kv = physical_grad
            if (problem == 'L4') kv = ratio*physical_grad
            if (problem == 'L3' .or. problem == 'L3a') then
                do i = 1, nl
                    kv(:, i) = kv(:, i) + (ratio - 1)*b* &
                               dot_product(b, physical_grad(:, i))
                end do
            end if
            if (dim == 2) then
                qgrad(1, :, 1) = physical_grad(1, :)
                qgrad(1, :, 2) = physical_grad(2, :)
                kt = 0
                do j = 1, 2
                    do i = 1, 2
                        kt(1, i, j) = (ratio - 1)*b(i)*b(j)
                    end do
                    kt(1, j, j) = kt(1, j, j) + 1
                end do
                qweight = w
                block
                    real(dp) :: block_matrix(3, 3)
                    call assemble_tensor_diffusion_matrix(qgrad, kt, qweight, &
                                                          block_matrix, status)
                    if (status%code /= 0) error stop 'tensor contraction failed'
                    local_matrix = local_matrix + block_matrix
                end block
            else
                local_matrix = local_matrix + w* &
                               matmul(transpose(physical_grad), kv)
            end if
            if (problem == 'L4') then
                do j = 1, nl
                    local_matrix(:, j) = local_matrix(:, j) + &
                                         w*val(:, q)*val(j, q)
                end do
            end if
            local_rhs(:, c) = local_rhs(:, c) + w*source*val(:, q)
        end do
        do j = 1, nl
            cc = maps(j, c)
            do i = 1, nl
                rr = maps(i, c)
                entry = (c - 1)*nl**2 + (j - 1)*nl + i
                rows(entry) = rr
                cols(entry) = cc
                entries(entry) = 0
                if (constrained(rr)) cycle
                if (constrained(cc)) then
                    local_rhs(i, c) = local_rhs(i, c) - &
                                      local_matrix(i, j)*bcvalue(cc)
                else
                    entries(entry) = local_matrix(i, j)
                end if
            end do
        end do
        deallocate (physical_grad, kv, local_matrix)
        if (dim == 2) deallocate (qgrad, kt, qweight)
    end do
    !$omp end parallel do
    rhs = 0
    !$omp parallel do private(c,i,index) schedule(static)
    do c = 1, nc
        do i = 1, nl
            index = maps(i, c)
            !$omp atomic update
            rhs(index) = rhs(index) + local_rhs(i, c)
        end do
    end do
    !$omp end parallel do
    entry = nl**2*nc
    do i = 1, nd
        rows(entry + i) = i
        cols(entry + i) = i
        entries(entry + i) = 0
        if (.not. constrained(i)) cycle
        entries(entry + i) = 1
        rhs(i) = bcvalue(i)
    end do
    call csc_from_triplet(nd, nd, rows, cols, entries, matrix, status)
    if (status%code /= 0) error stop 'lab CSC failed'
    deallocate (rows, cols, entries, local_rhs)
    solution = 0
    iterations = 0
    select case (trim(preconditioner))
    case ('ldlt')
        call factor%factor(matrix, status)
        if (status%code /= 0) then
            print *, status%message
            error stop 'LDLT factor failed'
        end if
        call factor%solve(rhs, solution, status)
        if (status%code /= 0) error stop 'LDLT solve failed'
    case ('ichol', 'ilu', 'jacobi')
        opts%method = 'pcg'
        opts%preconditioner = preconditioner
        opts%tolerance = 1.e-11_dp
        opts%max_iterations = 20000
        call solve_sparse(matrix, rhs, solution, opts, stats)
        if (.not. stats%converged) error stop 'PCG did not converge'
        iterations = stats%iterations
    case default
        error stop 'unknown preconditioner'
    end select
    residual = sqrt(sum((csc_matvec(matrix, solution) - rhs)**2))/ &
               max(sqrt(sum(rhs**2)), tiny(1.0_dp))
    if (residual > 2.e-8_dp) error stop 'linear residual failed'
    err = 0
    norm = 0
    herr = 0
    hnorm = 0
    aerr = 0
    anorm = 0
    center = 0
    !$omp parallel do default(shared) private(c,d,jac,inverse,det,ios,q,point, &
    !$omp exact,eg,source,b,value,gradient,w,t) &
    !$omp reduction(+:err,norm,herr,hnorm,aerr,anorm) schedule(static)
    do c = 1, nc
        jac = 0
        jac(3, 3) = 1
        do d = 1, dim
            jac(:, d) = vertices(:, cells(d + 1, c)) - vertices(:, cells(1, c))
        end do
        det = det3(jac)
        call inv3(jac, inverse, ios)
        do q = 1, nq
            point = vertices(:, cells(1, c)) + &
                    matmul(jac, [qx(q), qy(q), qz(q)])
            call problem_data(point, exact, eg, source, b)
            value = dot_product(solution(maps(:, c)), val(:, q))
            gradient = matmul(transpose(inverse), &
                              matmul(grad(:, :, q), solution(maps(:, c))))
            w = det*qw(q)
            err = err + w*(value - exact)**2
            norm = norm + w*exact**2
            herr = herr + w*sum((gradient - eg)**2)
            hnorm = hnorm + w*sum(eg**2)
            if (dim == 3) then
                call lab_transport(sqrt(sum(point(:2)**2)), &
                                   atan2(point(2), point(1)), point(3), t)
                if (t(9) < 0.04_dp) then
                    aerr = aerr + w*(value - exact)**2
                    anorm = anorm + w*exact**2
                end if
            end if
        end do
    end do
    !$omp end parallel do
    if (problem == 'L3a') then
        do i = 1, nv
            if (sum(vertices(:2, i)**2) < 1.e-20_dp) center = solution(i)
        end do
    end if
    write (*, '(a,1x,i0,1x,6(es24.16,1x),i0)') 'LAB', nd, &
        sqrt(err/norm), sqrt(herr/hnorm), &
        sqrt(aerr/max(anorm, tiny(1.0_dp))), center, residual, &
        real(matrix%nnz, dp), iterations
contains
    subroutine problem_data(x, u, gradient, forcing, direction)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: u, gradient(3), forcing, direction(3)
        real(dp) :: data(9), phi, r, bn
        direction = 0
        if (problem == 'L3a') then
            call lab_guenter(x(1), x(2), data(:8))
            u = data(1)
            gradient = [data(2:3), 0.0_dp]
            forcing = data(4)
            bn = sqrt(sum(data(5:6)**2))
            if (bn > 1.e-14_dp) direction = [data(5:6)/bn, 0.0_dp]
        else if (problem == 'L3') then
            r = sqrt(sum(x(:2)**2))
            phi = atan2(x(2), x(1))
            call lab_transport(r, phi, x(3), data)
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
            if (problem == 'L4') forcing = ratio*forcing + u
        end if
    end subroutine
end program
