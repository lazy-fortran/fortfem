program lab_hcurl_run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_lab_vector, only: lab_vector
    use fortfem_tetra_edge_dof_map, only: build_tetra_edge_dof_map
    use fortfem_tetra_nedelec_first_order, only: evaluate_tetra_nedelec_first_order
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortfem_assembly_tetra_nedelec_3d, only: &
        assemble_tetra_nedelec_curl_mass_csc, assemble_tetra_nedelec_vector_load_order
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_advanced_solvers, only: solve_sparse, solver_options_t, solver_stats_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, fortsparse_status_t
    use fortnum_linalg, only: det3
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    type(csc_t) :: full, matrix
    type(fortsparse_status_t) :: status
    type(native_ldlt_t) :: factor
    type(solver_options_t) :: opts
    type(solver_stats_t) :: stats
    integer, parameter :: local_edges(2, 6) = reshape([ &
                                            1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])
    integer, allocatable :: cells(:, :), faces(:, :), edges(:, :), map(:, :), sign(:, :)
    integer, allocatable :: rows(:), cols(:)
    logical, allocatable :: boundary(:)
    real(dp), allocatable :: vertices(:, :), rhs(:), solution(:), bc(:), entries(:)
    real(dp), allocatable :: qx(:), qy(:), qz(:), qw(:)
    real(dp) :: x(3), tangent(3), data(9), gx(8), gw(8), jac(3, 3), det
    real(dp) :: rv(3, 6), rc(3, 6), pv(3, 6), pc(3, 6), value(3), curl(3)
    real(dp) :: err, norm, cerr, cnorm, w, residual
    integer :: dim, nv, nc, nd, p, unit, i, j, c, e, f, q, ios, entry, iterations
    character(4096) :: arg, meshfile
    character(32) :: problem, preconditioner
    call get_command_argument(1, problem)
    call get_command_argument(2, arg)
    read (arg, *) p
    if (problem /= 'L2' .or. p /= 1) error stop 'Hcurl adapter requires L2 order 1'
    call get_command_argument(3, meshfile)
    call get_command_argument(4, preconditioner)
    open (newunit=unit, file=trim(meshfile), status='old')
    read (unit, *) dim, nv, nc
    if (dim /= 3) error stop 'Hcurl requires tetrahedra'
    allocate (vertices(3, nv), cells(4, nc), faces(4, nc))
    do i = 1, nv
        read (unit, *) vertices(:, i)
    end do
    do c = 1, nc
        read (unit, *) cells(:, c), faces(:, c)
    end do
    close (unit)
    call build_tetra_edge_dof_map(cells, edges, map, sign, ios)
    if (ios /= 0) error stop 'edge map failed'
    nd = size(edges, 2)
    allocate (boundary(nd), bc(nd), solution(nd))
    boundary = .false.
    do c = 1, nc
        do f = 1, 4
            if (faces(f, c) == 0) cycle
            do e = 1, 6
                if (any(local_edges(:, e) == f)) cycle
                boundary(map(e, c)) = .true.
            end do
        end do
    end do
    call gauss_legendre_ab(8, 0.0_dp, 1.0_dp, gx, gw)
    bc = 0
    !$omp parallel do private(e,tangent,q,x,data) schedule(static)
    do e = 1, nd
        if (.not. boundary(e)) cycle
        tangent = vertices(:, edges(2, e)) - vertices(:, edges(1, e))
        do q = 1, 8
            x = vertices(:, edges(1, e)) + gx(q)*tangent
            call lab_vector(x(1), x(2), x(3), data)
            bc(e) = bc(e) + gw(q)*dot_product(data(:3), tangent)
        end do
    end do
    !$omp end parallel do
    call assemble_tetra_nedelec_curl_mass_csc(vertices, cells, full, status, &
                                              1.0_dp, 1.0_dp, 1)
    if (status%code /= 0) error stop 'Hcurl assembly failed'
 call assemble_tetra_nedelec_vector_load_order(vertices, cells, 1, forcing, rhs, status)
    if (status%code /= 0) error stop 'Hcurl load failed'
    allocate (rows(full%nnz + nd), cols(full%nnz + nd), entries(full%nnz + nd))
    entry = 0
    do j = 1, nd
        do e = full%col_ptr(j), full%col_ptr(j + 1) - 1
            i = full%row_idx(e)
            entry = entry + 1
            rows(entry) = i
            cols(entry) = j
            entries(entry) = 0
            if (boundary(i)) cycle
            if (boundary(j)) then
                rhs(i) = rhs(i) - full%val(e)*bc(j)
            else
                entries(entry) = full%val(e)
            end if
        end do
    end do
    do i = 1, nd
        entry = entry + 1
        rows(entry) = i
        cols(entry) = i
        entries(entry) = 0
        if (boundary(i)) then
            entries(entry) = 1
            rhs(i) = bc(i)
        end if
    end do
    call csc_from_triplet(nd, nd, rows, cols, entries, matrix, status)
    if (status%code /= 0) error stop 'Hcurl constraints failed'
    solution = 0
    iterations = 0
    if (preconditioner == 'ldlt') then
        call factor%factor(matrix, status)
        if (status%code /= 0) error stop 'Hcurl LDLT failed'
        call factor%solve(rhs, solution, status)
        if (status%code /= 0) error stop 'Hcurl solve failed'
    else
        opts%method = 'pcg'
        opts%preconditioner = preconditioner
        opts%fill_level = 64
        opts%tolerance = 1.e-11_dp
        opts%max_iterations = 20000
        call solve_sparse(matrix, rhs, solution, opts, stats)
        if (.not. stats%converged) error stop 'Hcurl PCG failed'
        iterations = stats%iterations
    end if
    residual = sqrt(sum((csc_matvec(matrix, solution) - rhs)**2))/ &
               max(sqrt(sum(rhs**2)), tiny(1.0_dp))
    if (residual > 2.e-8_dp) error stop 'Hcurl residual failed'
    call tetra_duffy_quadrature(10, qx, qy, qz, qw, ios)
    if (ios /= 0) error stop 'Hcurl quadrature failed'
    err = 0
    norm = 0
    cerr = 0
    cnorm = 0
    !$omp parallel do private(c,j,jac,det,q,x,data,rv,rc,pv,pc,value,curl,w,ios) &
    !$omp reduction(+:err,norm,cerr,cnorm) schedule(static)
    do c = 1, nc
        do j = 1, 3
            jac(:, j) = vertices(:, cells(j + 1, c)) - vertices(:, cells(1, c))
        end do
        det = det3(jac)
        do q = 1, size(qw)
            x = vertices(:, cells(1, c)) + matmul(jac, [qx(q), qy(q), qz(q)])
            call lab_vector(x(1), x(2), x(3), data)
            call evaluate_tetra_nedelec_first_order([qx(q), qy(q), qz(q)], rv, rc, ios)
            if (ios /= 0) error stop 'Hcurl reference basis failed'
            call map_tetra_nedelec_covariant(jac, rv, rc, pv, pc, ios)
            if (ios /= 0) error stop 'Hcurl Piola map failed'
            value = matmul(pv, solution(map(:, c))*sign(:, c))
            curl = matmul(pc, solution(map(:, c))*sign(:, c))
            w = det*qw(q)
            err = err + w*sum((value - data(:3))**2)
            norm = norm + w*sum(data(:3)**2)
            cerr = cerr + w*sum((curl - data(4:6))**2)
            cnorm = cnorm + w*sum(data(4:6)**2)
        end do
    end do
    !$omp end parallel do
    write (*, '(a,1x,i0,1x,6(es24.16,1x),i0)') 'LAB', nd, sqrt(err/norm), &
     sqrt((err+cerr)/(norm+cnorm)),0.0_dp,0.0_dp,residual,real(matrix%nnz,dp),iterations
contains
    pure subroutine forcing(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)
        real(dp) :: data(9)
        call lab_vector(x, y, z, data)
        value = data(7:9)
    end subroutine
end program
