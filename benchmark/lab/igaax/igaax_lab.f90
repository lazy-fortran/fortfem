program igaax_lab
    ! IGA-AX laboratory driver: GLL, B-spline and spline x Fourier lines through
    ! the shared MF operator on the five-patch O-grid through the axis, with
    ! C1 polar splines as a comparison (criterion 2). One CSV row per run.
    ! Usage: fortfem_lab_igaax problem=L1|disk|aniso family=gll|spline|
    !   spline-fourier|gll-fourier|polar p=3 n=4 nphi=16 m=16 cont=-1 smooth=3
    !   chi=1 coarsen=full|semi tol=1e-10 header=1
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use omp_lib, only: omp_get_max_threads
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_GLL, LINE_BSPLINE
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_bspline_lines, only: build_bspline_line
    use fortfem_line_mesh, only: line_tensor_mesh
    use fortfem_line_hierarchy, only: ogrid_line_hierarchy, box_line_hierarchy, &
                                      PHI_NONE, PHI_LOCAL, PHI_FOURIER
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                            tensor_load, tensor_values, tensor_boundary_projection, &
                                       apply_tensor_operator, assemble_tensor_operator
    use fortfem_tensor_solver, only: tensor_multigrid_t, tensor_solve_info_t, &
                                     build_line_multigrid, solve_tensor_pcg
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_bspline_feec, only: evaluate_bspline_basis
    use fortfem_bspline_polar, only: build_bspline_polar_h1_extraction
    use fortfem_assembly_bspline_polar_2d, only: restrict_bspline_polar_operator_csc
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none
    interface
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(len=1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(len=32) :: problem = 'L1', family = 'spline', coarsen = 'full'
    integer :: p = 3, n = 4, nphi = 16, m = 16, cont = -1, smooth = 3, header = 0
    integer :: levels, iterations, status, dim, dofs
    real(dp) :: chi = 1, tol = 1.e-10_dp, residual
    real(dp) :: t0, t1, t2, rel_l2, abs_l2, axis_max, axis_l2, norm_exact
    type(tensor_mesh_t), allocatable :: meshes(:)
    type(csr_matrix_t), allocatable :: transfers(:)
    real(dp), allocatable :: polar_t(:), polar_rho(:)
    integer :: polar_ntheta
    call parse_arguments()
    dim = merge(3, 2, problem == 'L1')
    t0 = wall()
    if (family == 'polar') then
        call polar_run()
    else
        call line_run()
    end if
    if (header == 1) print '(a)', 'problem,family,p,n,nphi,m,continuity,'// &
        'coarsen,chi,dofs,levels,iterations,residual,setup_s,solve_s,total_s,'// &
        'rel_l2,abs_l2,axis_max,axis_l2,threads,rss_kib'
    print '(a,",",a,",",i0,",",i0,",",i0,",",i0,",",i0,",",a,",",es9.2,",",i0,'// &
        '",",i0,",",i0,",",es10.3,",",f0.4,",",f0.4,",",f0.4,",",4(es11.4,","),'// &
        'i0,",",i0)', trim(problem), trim(family), p, n, nphi, m, cont, &
        trim(coarsen), chi, dofs, levels, iterations, residual, t1 - t0, t2 - t1, &
        t2 - t0, rel_l2, abs_l2, axis_max, axis_l2, omp_get_max_threads(), peak_rss()
contains
    function wall() result(t)
        real(dp) :: t
        integer(int64) :: count, rate
        call system_clock(count, rate)
        t = real(count, dp)/real(rate, dp)
    end function

    integer function peak_rss() result(kib)
        ! VmHWM from /proc (Linux); -1 elsewhere.
        character(len=256) :: line
        integer :: unit, ios
        kib = -1
        open (newunit=unit, file='/proc/self/status', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (line(1:6) == 'VmHWM:') read (line(7:), *) kib
        end do
        close (unit)
    end function

    subroutine parse_arguments()
        character(len=64) :: arg, key, value
        integer :: i, eq
        do i = 1, command_argument_count()
            call get_command_argument(i, arg)
            eq = index(arg, '=')
            if (eq == 0) error stop 'arguments are key=value'
            key = arg(:eq - 1)
            value = arg(eq + 1:)
            select case (trim(key))
            case ('problem'); problem = value
            case ('family'); family = value
            case ('coarsen'); coarsen = value
            case ('p'); read (value, *) p
            case ('n'); read (value, *) n
            case ('nphi'); read (value, *) nphi
            case ('m'); read (value, *) m
            case ('cont'); read (value, *) cont
            case ('smooth'); read (value, *) smooth
            case ('header'); read (value, *) header
            case ('chi'); read (value, *) chi
            case ('tol'); read (value, *) tol
            case default; error stop 'unknown argument'
            end select
        end do
        if (cont < 0) cont = p - 1
    end subroutine

    ! ---- L1 torus (benchmark/lab/PROBLEMS.md), disk MMS, aligned anisotropy ----
    subroutine torus(xi, x)
        ! (rx, ry) is the O-grid disk point; (u, v) = (ry, rx) is right-handed.
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: a, rc, zc, r
        a = 1 + 0.1_dp*cos(5*xi(3))
        rc = 3 + 0.1_dp*cos(5*xi(3))
        zc = 0.1_dp*sin(5*xi(3))
        r = rc + a*xi(2)
        x = [r*cos(xi(3)), r*sin(xi(3)), zc + 1.3_dp*xi(1)]
    end subroutine

    real(dp) function axis_distance2(x) result(s)
        ! Squared distance to the axis in the unit-disk chart (s = u^2 + v^2).
        real(dp), intent(in) :: x(:)
        real(dp) :: phi, a, rc, zc
        if (size(x) == 2) then
            s = sum(x**2)
            return
        end if
        phi = atan2(x(2), x(1))
        a = 1 + 0.1_dp*cos(5*phi)
        rc = 3 + 0.1_dp*cos(5*phi)
        zc = 0.1_dp*sin(5*phi)
        s = ((norm2(x(1:2)) - rc)/a)**2 + ((x(3) - zc)/1.3_dp)**2
    end function

    function exact(x) result(u)
        real(dp), intent(in) :: x(:)
        real(dp) :: u
        if (size(x) == 3) then
            u = sin(x(1)) + cos(x(2)) + sin(x(3)) + sum(x)/10
        else if (problem == 'aniso') then
            u = sin(pi*x(1))*sin(pi*x(2))
        else
            u = exp(0.1_dp*x(1) + 0.2_dp*x(2))
        end if
    end function

    function forcing(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        if (size(x) == 3) then
            f = sin(x(1)) + cos(x(2)) + sin(x(3))
        else if (problem == 'aniso') then
            f = (chi + 1)*pi**2*exact(x)
        else
            f = -0.05_dp*exact(x)
        end if
    end function

    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
        if (problem == 'aniso') k(1, 1) = chi
        c = 0
    end subroutine

    subroutine errors(op, x)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: values(:, :)
        real(dp) :: e, u
        integer :: c, q
        call tensor_values(op, x, values)
        abs_l2 = 0; norm_exact = 0; axis_max = 0; axis_l2 = 0
        do c = 1, op%mesh%nelem
            do q = 1, size(values, 1)
                u = exact(op%points(:, q, c))
                e = values(q, c) - u
                abs_l2 = abs_l2 + op%volume(q, c)*e**2
                norm_exact = norm_exact + op%volume(q, c)*u**2
                if (axis_distance2(op%points(:, q, c)) < 0.04_dp) then
                    axis_max = max(axis_max, abs(e))
                    axis_l2 = axis_l2 + op%volume(q, c)*e**2
                end if
            end do
        end do
        abs_l2 = sqrt(abs_l2)
        axis_l2 = sqrt(axis_l2)
        rel_l2 = abs_l2/sqrt(norm_exact)
    end subroutine

    subroutine line_run()
        type(tensor_multigrid_t) :: mg
        type(tensor_solve_info_t) :: info
        real(dp), allocatable :: rhs(:), x(:), bc(:), check(:)
        integer :: phi, fam
        fam = LINE_BSPLINE
        if (family(1:3) == 'gll') fam = LINE_GLL
        phi = PHI_NONE
        if (dim == 3) phi = PHI_LOCAL
        if (index(family, 'fourier') > 0) phi = PHI_FOURIER
        if (problem == 'aniso') then
            call box_line_hierarchy(2, p, [n, n, 1], &
                                    [.true., coarsen /= 'semi', .false.], meshes, &
                                    transfers, status, continuity=cont)
        else if (dim == 3) then
            call ogrid_line_hierarchy(3, fam, p, n, phi, merge(m, nphi, &
                                      phi == PHI_FOURIER), meshes, transfers, status, &
                                      map=torus, continuity=cont)
        else
            call ogrid_line_hierarchy(2, fam, p, n, PHI_NONE, 0, meshes, transfers, &
                                      status, continuity=cont)
        end if
        if (status /= 0) error stop 'hierarchy'
        levels = size(meshes)
        dofs = meshes(1)%ndof
        call build_line_multigrid(meshes, transfers, coefficient, mg, status, &
                                  smoothing=smooth)
        if (status /= 0) error stop 'multigrid'
        call tensor_load(mg%level(1)%op, forcing, rhs)
        call tensor_boundary_projection(mg%level(1)%op, exact, bc, status)
        if (status /= 0) error stop 'boundary projection'
        allocate (x(dofs), check(dofs))
        x = 0
        t1 = wall()
        call solve_tensor_pcg(mg, rhs, tol, x, info, bc=bc, maxiter=2000)
        t2 = wall()
        if (info%status /= 0) error stop 'PCG'
        iterations = info%iterations
        residual = info%residual
        call errors(mg%level(1)%op, x)
    end subroutine

    ! ---- C1 polar splines (Toshniwal-Hughes Type-1 extraction), comparison ----
    subroutine polar_map(xi, x)
        ! Spline polar disk: ring i at radius rho_i, control angles of the
        ! extraction. xi = (theta, r) with theta first (extraction ordering).
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp), allocatable :: bt(:), dt(:), br(:), dr(:)
        real(dp) :: angle, tr(size(polar_rho) + p + 1)
        integer :: i, j, jj, s
        call evaluate_bspline_basis(polar_t, p, modulo(xi(1), 2*pi), bt, dt, s)
        if (s /= 0) error stop 'polar angular basis'
        tr(:p + 1) = 0
        do i = 1, size(polar_rho) - p - 1
            tr(p + 1 + i) = real(i, dp)/(size(polar_rho) - p)
        end do
        tr(size(polar_rho) + 1:) = 1
        call evaluate_bspline_basis(tr, p, xi(2), br, dr, s)
        if (s /= 0) error stop 'polar radial basis'
        x = 0
        do j = 1, size(bt)
            jj = modulo(j - 1, polar_ntheta) + 1
            angle = 2*pi - real(2*jj - 1, dp)*pi/polar_ntheta
            do i = 1, size(br)
                x = x + bt(j)*br(i)*polar_rho(i)*[cos(angle), sin(angle)]
            end do
        end do
    end subroutine

    subroutine polar_run()
        ! Tensor space (periodic angle x open radius) on the polar spline map,
        ! the shared MF operator assembled to CSR, restricted by E (C1 at the
        ! pole) and solved by dense Cholesky (comparison sizes only; the
        ! polar timing is not comparable with the MF rows).
        type(tensor_line_tables_t) :: lines(3)
        type(tensor_mesh_t) :: mesh
        type(tensor_operator_t) :: op
        type(csr_matrix_t) :: a
        type(csc_t) :: tensor, polar
        type(fortsparse_status_t) :: info
        real(dp), allocatable :: extraction(:, :), rhs(:), bc(:), ax(:), xp(:, :), &
                                 bp(:), x(:), dense(:, :)
        integer :: i, j, k, nr, ntheta, ndof
        nr = n
        ntheta = nphi
        polar_ntheta = ntheta
        call build_bspline_line(p, lines(1), status, periodic=.true., &
                                knots=[(2*pi*i/ntheta, i=0, ntheta)])
        if (status /= 0) error stop 'polar angular line'
        call build_bspline_line(p, lines(2), status, nelem=nr)
        if (status /= 0) error stop 'polar radial line'
        ! SPL1's periodic knot vector on [0, 2 pi), for geometry evaluation.
        allocate (polar_t(ntheta + 2*p + 1))
        do i = 1, size(polar_t)
            j = i - p - 1
            polar_t(i) = 2*pi*real(j, dp)/ntheta
        end do
        allocate (polar_rho(lines(2)%ndof))
        do i = 1, lines(2)%ndof
            polar_rho(i) = 0
            do j = 1, p
                polar_rho(i) = polar_rho(i) + &
                               min(1.0_dp, max(0.0_dp, real(i + j - p - 1, dp)/nr))
            end do
            polar_rho(i) = polar_rho(i)/p
        end do
        call line_tensor_mesh(2, lines, mesh, status, map=polar_map)
        if (status /= 0) error stop 'polar mesh'
        ! Only the outer ring is Dirichlet; the pole ring is constrained by E.
        mesh%boundary(1:ntheta) = .false.
        call build_tensor_operator(mesh, coefficient, op, status)
        if (status /= 0) error stop 'polar operator'
        ndof = mesh%ndof
        dofs = ndof - 2*ntheta + 3
        levels = 1
        call tensor_load(op, forcing, rhs)
        call tensor_boundary_projection(op, exact, bc, status)
        if (status /= 0) error stop 'polar boundary projection'
        allocate (ax(ndof))
        call apply_tensor_operator(op, bc, ax, eliminate=.false.)
        rhs = rhs - ax
        where (mesh%boundary) rhs = 0
        call assemble_tensor_operator(op, a)
        tensor%nrow = a%n; tensor%ncol = a%n; tensor%nnz = size(a%val)
        tensor%col_ptr = a%rowptr; tensor%row_idx = a%colidx; tensor%val = a%val
        call build_bspline_polar_h1_extraction(ntheta, lines(2)%ndof, extraction, &
                                               status)
        if (status /= 0) error stop 'polar extraction'
        call restrict_bspline_polar_operator_csc(extraction, tensor, polar, info)
        if (info%code /= 0) error stop 'polar restriction'
        bp = matmul(extraction, rhs)
        allocate (dense(polar%nrow, polar%nrow), xp(polar%nrow, 1))
        dense = 0
        do j = 1, polar%ncol
            do k = polar%col_ptr(j), polar%col_ptr(j + 1) - 1
                dense(polar%row_idx(k), j) = polar%val(k)
            end do
        end do
        dense = (dense + transpose(dense))/2
        xp(:, 1) = bp
        t1 = wall()
        call dposv('U', polar%nrow, 1, dense, polar%nrow, xp, polar%nrow, status)
        if (status /= 0) error stop 'polar Cholesky'
        t2 = wall()
        x = matmul(transpose(extraction), xp(:, 1))
        where (mesh%boundary) x = bc
        iterations = 0
        call apply_tensor_operator(op, x - bc, ax)
        residual = norm2(matmul(extraction, ax - rhs))/norm2(bp)
        call errors(op, x)
    end subroutine
end program
