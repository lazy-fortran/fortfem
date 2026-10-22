program lab_bspline_run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_lab_scalar, only: lab_scalar
    use fortfem_generated_lab_square_disk, only: lab_square_disk
    use fortfem_generated_lab_geometry, only: lab_geometry
    use fortfem_generated_lab_transport, only: lab_transport
    use fortfem_bspline_feec, only: evaluate_bspline_basis, &
                                    evaluate_nurbs_volume_geometry
    use fortfem_assembly_bspline_3d, only: assemble_bspline_h1_operator_3d_csc
    use fortfem_advanced_solvers, only: solve_sparse, solver_options_t, solver_stats_t
    use fortsparse_native_ldlt, only: native_ldlt_t
    use fortsparse, only: csc_t, csc_from_triplet, csc_matvec, fortsparse_status_t
    use fortnum_linalg, only: inv3, det3, dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    type(csc_t) :: full, matrix
    type(fortsparse_status_t) :: status
    type(native_ldlt_t) :: factor
    type(solver_options_t) :: opts
    type(solver_stats_t) :: stats
    real(dp), allocatable :: ku(:), kv(:), kp(:), controls(:, :, :, :), weights(:, :, :)
 real(dp), allocatable :: gu(:), gv(:), gp(:)
    real(dp), allocatable :: coeff(:, :, :), rhs(:), solution(:), bc(:)
    real(dp), allocatable :: bu(:, :), bv(:, :), bp(:, :), line(:), output(:)
    real(dp), allocatable :: entries(:)
    integer, allocatable :: map(:), rows(:), cols(:)
    logical, allocatable :: boundary(:)
    real(dp) :: point(3), jac(3, 3), geom(12), data(8), disk(2)
    real(dp) :: u, v, phi, err, norm, herr, hnorm
    real(dp) :: aerr, anorm, residual, pi
    integer :: p, n, nu, nv, np, nd, nfull, nlocal, nq, ios, i, j, k, c, r, entry
    integer :: l, raw, reduced, iterations
    character(4096) :: arg
    character(32) :: problem, preconditioner
    call get_command_argument(1, problem)
    if (problem /= 'L1' .and. problem /= 'L5') error stop 'spline problem unsupported'
    call get_command_argument(2, arg)
    read (arg, *) p
    call get_command_argument(3, arg)
    read (arg, *) n
    call get_command_argument(4, preconditioner)
    if (p < 2 .or. p > 4 .or. n < 1) error stop 'invalid spline parameters'
    pi = acos(-1.0_dp)
    call knots(2*n, p, ku, gu)
    call knots(2*n, p, kv, gv)
    call knots(20*n, p, kp, gp)
    nu = size(gu)
    nv = size(gv)
    np = size(gp)
    nfull = nu*nv*np
    nd = nu*nv*(np - 1)
    nlocal = (p + 1)**3
    nq = p + 3
    allocate (controls(3, nu, nv, np), weights(nu, nv, np), coeff(nu, nv, np))
    allocate (map(nfull), boundary(nd), bc(nd), rhs(nd), solution(nd))
    weights = 1
    boundary = .false.
    do k = 1, np
        phi = 2*pi*gp(k)
        do j = 1, nv
            do i = 1, nu
                ! Square-to-disk map used only to position geometry controls.
                call lab_square_disk(2*gu(i) - 1, 2*gv(j) - 1, disk)
                u = disk(1)
                v = -disk(2)
                call lab_geometry(u, v, phi, geom)
                controls(:, i, j, k) = [geom(1)*cos(phi), geom(1)*sin(phi), geom(2)]
                raw = idx(i, j, k)
                reduced = idx(i, j, 1 + modulo(k - 1, np - 1))
                map(raw) = reduced
                if (i == 1 .or. i == nu .or. j == 1 .or. j == nv) &
                    boundary(reduced) = .true.
            end do
        end do
    end do
    ! Interpolate the exact trace in the tensor spline space, rather than
    ! treating control coordinates or Greville samples as scalar coefficients.
    call collocation(ku, gu, bu)
    call collocation(kv, gv, bv)
    call collocation(kp, gp, bp)
    do k = 1, np
        do j = 1, nv
            do i = 1, nu
                call evaluate_nurbs_volume_geometry(ku, kv, kp, p, p, p, &
                                controls, weights, gu(i), gv(j), gp(k), point, jac, ios)
                if (ios /= 0) error stop 'spline trace geometry failed'
                call lab_scalar(point(1), point(2), point(3), data)
                coeff(i, j, k) = data(1)
            end do
        end do
    end do
    allocate (line(nu), output(nu))
    do k = 1, np
        do j = 1, nv
            line = coeff(:, j, k)
            call dense_solve(bu, line, output, ios)
            if (ios /= 0) error stop 'trace interpolation u failed'
            coeff(:, j, k) = output
        end do
    end do
    deallocate (line, output)
    allocate (line(nv), output(nv))
    do k = 1, np
        do i = 1, nu
            line = coeff(i, :, k)
            call dense_solve(bv, line, output, ios)
            if (ios /= 0) error stop 'trace interpolation v failed'
            coeff(i, :, k) = output
        end do
    end do
    deallocate (line, output)
    allocate (line(np), output(np))
    do j = 1, nv
        do i = 1, nu
            line = coeff(i, j, :)
            call dense_solve(bp, line, output, ios)
            if (ios /= 0) error stop 'trace interpolation phi failed'
            coeff(i, j, :) = output
        end do
    end do
    bc = 0
    do k = 1, np
        do j = 1, nv
            do i = 1, nu
                bc(map(idx(i, j, k))) = coeff(i, j, k)
            end do
        end do
    end do
    call assemble_bspline_h1_operator_3d_csc(ku, kv, kp, p, p, p, controls, &
                                             weights, nq, full, status)
    if (status%code /= 0) error stop 'existing spline H1 assembly failed'
    ! Collapse only the two toroidal endpoint planes. C0 periodicity, no cut.
    allocate (rows(full%nnz + nd), cols(full%nnz + nd), entries(full%nnz + nd))
    entry = 0
    rhs = 0
    do c = 1, nfull
        do l = full%col_ptr(c), full%col_ptr(c + 1) - 1
            r = map(full%row_idx(l))
            reduced = map(c)
            entry = entry + 1
            rows(entry) = r
            cols(entry) = reduced
            entries(entry) = 0
            if (boundary(r)) cycle
            if (boundary(reduced)) then
                rhs(r) = rhs(r) - full%val(l)*bc(reduced)
            else
                entries(entry) = full%val(l)
            end if
        end do
    end do
    do i = 1, nd
        entry = entry + 1
        rows(entry) = i
        cols(entry) = i
        entries(entry) = 0
        if (boundary(i)) entries(entry) = 1
    end do
    call csc_from_triplet(nd, nd, rows, cols, entries, matrix, status)
    if (status%code /= 0) error stop 'spline periodic reduction failed'
    deallocate (rows, cols, entries)
    call integrate(.false.)
    where (boundary) rhs = bc
    solution = 0
    iterations = 0
    if (preconditioner == 'ldlt') then
        call factor%factor(matrix, status)
        if (status%code /= 0) error stop 'spline LDLT factor failed'
        call factor%solve(rhs, solution, status)
        if (status%code /= 0) error stop 'spline LDLT solve failed'
    else
        opts%method = 'pcg'
        opts%preconditioner = preconditioner
        opts%tolerance = 1.e-11_dp
        opts%max_iterations = 20000
        call solve_sparse(matrix, rhs, solution, opts, stats)
        if (.not. stats%converged) error stop 'spline PCG failed'
        iterations = stats%iterations
    end if
    residual = sqrt(sum((csc_matvec(matrix, solution) - rhs)**2))/ &
               max(sqrt(sum(rhs**2)), tiny(1.0_dp))
    if (residual > 2.e-8_dp) error stop 'spline linear residual failed'
    err = 0
    norm = 0
    herr = 0
    hnorm = 0
    aerr = 0
    anorm = 0
    call integrate(.true.)
    write (*, '(a,1x,i0,1x,6(es24.16,1x),i0)') 'LAB', nd, sqrt(err/norm), &
        sqrt(herr/hnorm), sqrt(aerr/max(anorm, tiny(1.0_dp))), 0.0_dp, &
        residual, real(matrix%nnz, dp), iterations
contains
    integer function idx(i, j, k) result(index)
        integer, intent(in) :: i, j, k
        index = i + nu*((j - 1) + nv*(k - 1))
    end function
    subroutine knots(spans, degree, knot, greville)
        integer, intent(in) :: spans, degree
        real(dp), allocatable, intent(out) :: knot(:), greville(:)
        integer :: i, nb
        nb = spans + degree
        allocate (knot(nb + degree + 1), greville(nb))
        knot(:degree + 1) = 0
        do i = 1, spans - 1
            knot(degree + 1 + i) = real(i, dp)/spans
        end do
        knot(nb + 1:) = 1
        do i = 1, nb
            greville(i) = sum(knot(i + 1:i + degree))/degree
        end do
    end subroutine
    subroutine collocation(knot, greville, matrix)
        real(dp), intent(in) :: knot(:), greville(:)
        real(dp), allocatable, intent(out) :: matrix(:, :)
        real(dp), allocatable :: v(:), d(:)
        integer :: i, status
        allocate (matrix(size(greville), size(greville)))
        do i = 1, size(greville)
            call evaluate_bspline_basis(knot, p, greville(i), v, d, status)
            if (status /= 0) error stop 'collocation basis failed'
            matrix(i, :) = v
        end do
    end subroutine
    subroutine integrate(measure)
        logical, intent(in) :: measure
        integer :: cell, ncell, sx, sy, sz, ix, iy, iz, aa, bb, cc, dof
        integer, allocatable :: local(:)
        real(dp), allocatable :: xq(:), yq(:), zq(:), xw(:), yw(:), zw(:)
        real(dp), allocatable :: vx(:), vy(:), vz(:), dx(:), dy(:), dz(:)
        real(dp), allocatable :: bv(:), bg(:, :), load(:)
        real(dp) :: x(3), jac(3, 3), inv(3, 3), determinant, data(8), t(9)
        real(dp) :: w, value, gradient(3)
        integer :: status
        ncell = (2*n)**2*20*n
        !$omp parallel do default(shared) private(cell,sx,sy,sz,ix,iy,iz,aa,bb,cc, &
        !$omp dof,local,xq,yq,zq,xw,yw,zw,vx,vy,vz,dx,dy,dz,bv,bg,load, &
        !$omp x,jac,inv,determinant,data,t,w,value,gradient,status) &
        !$omp reduction(+:err,norm,herr,hnorm,aerr,anorm) schedule(static)
        do cell = 0, ncell - 1
            sx = p + 1 + modulo(cell, 2*n)
            sy = p + 1 + modulo(cell/(2*n), 2*n)
            sz = p + 1 + cell/(2*n)**2
            allocate (xq(nq), yq(nq), zq(nq), xw(nq), yw(nq), zw(nq))
            allocate (local(nlocal), bv(nlocal), bg(3, nlocal), load(nlocal))
            call gauss_legendre_ab(nq, ku(sx), ku(sx + 1), xq, xw)
            call gauss_legendre_ab(nq, kv(sy), kv(sy + 1), yq, yw)
            call gauss_legendre_ab(nq, kp(sz), kp(sz + 1), zq, zw)
            dof = 0
            do cc = sz - p, sz
                do bb = sy - p, sy
                    do aa = sx - p, sx
                        dof = dof + 1
                        local(dof) = map(idx(aa, bb, cc))
                    end do
                end do
            end do
            load = 0
            do iz = 1, nq
                call evaluate_bspline_basis(kp, p, zq(iz), vz, dz, status)
                if (status /= 0) error stop 'spline phi basis failed'
                do iy = 1, nq
                    call evaluate_bspline_basis(kv, p, yq(iy), vy, dy, status)
                    if (status /= 0) error stop 'spline v basis failed'
                    do ix = 1, nq
                        call evaluate_bspline_basis(ku, p, xq(ix), vx, dx, status)
                        if (status /= 0) error stop 'spline u basis failed'
                        call evaluate_nurbs_volume_geometry(ku, kv, kp, p, p, p, &
                              controls, weights, xq(ix), yq(iy), zq(iz), x, jac, status)
                        if (status /= 0) error stop 'spline volume geometry failed'
                        determinant = det3(jac)
                        call inv3(jac, inv, status)
                        if (status /= 0 .or. determinant <= 0) &
                            error stop 'spline Jacobian invalid'
                        w = determinant*xw(ix)*yw(iy)*zw(iz)
                        dof = 0
                        do cc = sz - p, sz
                            do bb = sy - p, sy
                                do aa = sx - p, sx
                                    dof = dof + 1
                                    bv(dof) = vx(aa)*vy(bb)*vz(cc)
                                    bg(:, dof) = [dx(aa)*vy(bb)*vz(cc), &
                                             vx(aa)*dy(bb)*vz(cc), vx(aa)*vy(bb)*dz(cc)]
                                end do
                            end do
                        end do
                        call lab_scalar(x(1), x(2), x(3), data)
                        if (measure) then
                            value = dot_product(solution(local), bv)
                          gradient = matmul(transpose(inv), matmul(bg, solution(local)))
                            err = err + w*(value - data(1))**2
                            norm = norm + w*data(1)**2
                            herr = herr + w*sum((gradient - data(2:4))**2)
                            hnorm = hnorm + w*sum(data(2:4)**2)
                     call lab_transport(sqrt(sum(x(:2)**2)), atan2(x(2), x(1)), x(3), t)
                            if (t(9) < .04_dp) then
                                aerr = aerr + w*(value - data(1))**2
                                anorm = anorm + w*data(1)**2
                            end if
                        else
                            load = load + w*data(5)*bv
                        end if
                    end do
                end do
            end do
            if (.not. measure) then
                do dof = 1, nlocal
                    aa = local(dof)
                    !$omp atomic update
                    rhs(aa) = rhs(aa) + load(dof)
                end do
            end if
            deallocate (xq, yq, zq, xw, yw, zw, local, bv, bg, load)
        end do
        !$omp end parallel do
    end subroutine
end program
