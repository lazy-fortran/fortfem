program test_lab_line_operator
    ! Spline and Fourier lines through the shared matrix-free tensor operator:
    ! agreement with the existing assembled B-spline H1 operators, with the
    ! unfactored CSR oracle for mixed products, exactness on spline x trig
    ! polynomials and the Galerkin identity P^T A_fine P = A_coarse.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    use fortfem_bspline_lines, only: build_bspline_line
    use fortfem_fourier_lines, only: build_fourier_line
    use fortfem_gll_line, only: build_gll_line
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_line_mesh, only: line_tensor_mesh, line_ogrid_mesh, &
                                 line_mesh_prolongation, line_matrix_t
    use fortfem_line_transfer, only: bspline_line_prolongation, &
                                     fourier_line_prolongation
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                     apply_tensor_operator, assemble_tensor_operator, tensor_diagonal, &
                                       tensor_load, tensor_l2_error
    use fortfem_tensor_solver, only: tensor_multigrid_t, tensor_solve_info_t, &
                                     build_line_multigrid, solve_tensor_pcg
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_bspline_feec, only: evaluate_nurbs_surface_geometry, &
                                    evaluate_nurbs_volume_geometry
    use fortfem_assembly_bspline_2d, only: assemble_bspline_h1_operator_csc
    use fortfem_assembly_bspline_3d, only: assemble_bspline_h1_operator_3d_csc
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: knots(:), net2(:, :, :), net3(:, :, :, :)
    integer :: degree, nspan
    call assembled_spline_agreement()
    call mixed_oracle_agreement()
    call trig_polynomial_exactness()
    call galerkin_nestedness()
    print '(a)', 'PASS: spline/Fourier lines through the shared MF operator'
contains
    subroutine uniform_open_knots(p, n)
        integer, intent(in) :: p, n
        integer :: j
        allocate (knots(2*(p + 1) + n - 1))
        knots(:p + 1) = 0
        do j = 1, n - 1
            knots(p + 1 + j) = real(j, dp)/n
        end do
        knots(p + n + 1:) = 1
    end subroutine

    real(dp) function greville(i, p)
        integer, intent(in) :: i, p
        greville = sum(knots(i + 1:i + p))/p
    end function

    subroutine spline_map2(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: jac(2, 2), w(size(net2, 2), size(net2, 3))
        integer :: status
        w = 1
        call evaluate_nurbs_surface_geometry(knots, knots, degree, degree, net2, w, &
                                             xi(1), xi(2), x, jac, status)
        if (status /= 0) error stop 'surface geometry'
    end subroutine

    subroutine spline_map3(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: jac(3, 3), w(size(net3, 2), size(net3, 3), size(net3, 4))
        integer :: status
        w = 1
        call evaluate_nurbs_volume_geometry(knots, knots, knots, degree, degree, &
                                   degree, net3, w, xi(1), xi(2), xi(3), x, jac, status)
        if (status /= 0) error stop 'volume geometry'
    end subroutine

    pure subroutine anisotropic(point, value)
        real(dp), intent(in) :: point(2)
        real(dp), intent(out) :: value(2, 2)
        value = reshape([1.3_dp + point(1)**2, 0.2_dp, 0.2_dp, 0.8_dp + point(2)], &
                        [2, 2])
    end subroutine

    subroutine anisotropic_coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        call anisotropic(x(1:2), k)
        c = 0.7_dp
    end subroutine

    subroutine isotropic_coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1.6_dp
        end do
        c = 0.45_dp + 0*x(1)
    end subroutine

    subroutine laplace_coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
        c = 0*x(1)
    end subroutine

    function csc_times(a, x) result(y)
        type(csc_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp) :: y(a%nrow)
        integer :: j, k
        y = 0
        do j = 1, a%ncol
            do k = a%col_ptr(j), a%col_ptr(j + 1) - 1
                y(a%row_idx(k)) = y(a%row_idx(k)) + a%val(k)*x(j)
            end do
        end do
    end function

    function csr_times(a, x) result(y)
        type(csr_matrix_t), intent(in) :: a
        real(dp), intent(in) :: x(:)
        real(dp) :: y(a%n)
        integer :: i, k
        y = 0
        do i = 1, a%n
            do k = a%rowptr(i), a%rowptr(i + 1) - 1
                y(i) = y(i) + a%val(k)*x(a%colidx(k))
            end do
        end do
    end function

    subroutine assembled_spline_agreement()
        ! MF spline apply equals FortFEM's assembled isogeometric H1 operators on
        ! a curved polynomial spline geometry (sampled exactly at degree q = p).
        type(tensor_line_tables_t) :: lines(3)
        type(tensor_mesh_t) :: mesh
        type(tensor_operator_t) :: op
        type(csc_t) :: a
        type(fortsparse_status_t) :: info
        real(dp), allocatable :: x(:), y(:), z(:), w2(:, :), w3(:, :, :)
        integer :: status, i, j, k, n, dim
        real(dp) :: err
        do dim = 2, 3
            do degree = 1, 4
                nspan = 3
                if (allocated(knots)) deallocate (knots)
                call uniform_open_knots(degree, nspan)
                n = size(knots) - degree - 1
                do i = 1, dim
                    call build_bspline_line(degree, lines(i), status, nelem=nspan)
                    if (status /= 0) error stop 'spline line'
                end do
                if (dim == 2) then
                    if (allocated(net2)) deallocate (net2)
                    allocate (net2(2, n, n), w2(n, n))
                    w2 = 1
                    do j = 1, n
                        do i = 1, n
                            net2(:, i, j) = [greville(i, degree) + 0.08_dp* &
                                      sin(real(i*j, dp)), 0.7_dp*greville(j, degree) + &
                                             0.2_dp*greville(i, degree)**2]
                        end do
                    end do
                    call line_tensor_mesh(2, lines, mesh, status, map=spline_map2)
                    if (status /= 0) error stop 'spline mesh 2D'
                   call build_tensor_operator(mesh, anisotropic_coefficient, op, status)
                    if (status /= 0) error stop 'spline operator 2D'
                   call assemble_bspline_h1_operator_csc(knots, knots, degree, degree, &
                               net2, w2, degree + 1, a, info, mass_coefficient=0.7_dp, &
                                                  stiffness_tensor_function=anisotropic)
                    deallocate (w2)
                else
                    if (allocated(net3)) deallocate (net3)
                    allocate (net3(3, n, n, n), w3(n, n, n))
                    w3 = 1
                    do k = 1, n
                        do j = 1, n
                            do i = 1, n
                                net3(:, i, j, k) = [greville(i, degree) + 0.05_dp* &
                                          cos(real(i + 2*k, dp)), greville(j, degree), &
                                                  1.2_dp*greville(k, degree) + 0.1_dp* &
                                                greville(i, degree)*greville(j, degree)]
                            end do
                        end do
                    end do
                    call line_tensor_mesh(3, lines, mesh, status, map=spline_map3)
                    if (status /= 0) error stop 'spline mesh 3D'
                    call build_tensor_operator(mesh, isotropic_coefficient, op, status)
                    if (status /= 0) error stop 'spline operator 3D'
                 call assemble_bspline_h1_operator_3d_csc(knots, knots, knots, degree, &
                                        degree, degree, net3, w3, degree + 1, a, info, &
                                 stiffness_coefficient=1.6_dp, mass_coefficient=0.45_dp)
                    deallocate (w3)
                end if
                if (info%code /= 0) error stop 'assembled spline operator'
                if (a%nrow /= mesh%ndof) error stop 'spline dof count'
                allocate (x(mesh%ndof), y(mesh%ndof))
                do i = 1, mesh%ndof
                    x(i) = cos(real(7*i, dp))
                end do
                call apply_tensor_operator(op, x, y, eliminate=.false.)
                z = csc_times(a, x)
                err = norm2(y - z)/norm2(z)
                print '(a,i1,a,i1,a,es10.2)', 'assembled spline dim=', dim, ' p=', &
                    degree, ' relative=', err
                if (err > 1.e-12_dp) error stop 'MF spline versus assembled B-spline'
                deallocate (x, y, z)
            end do
        end do
    end subroutine

    subroutine annulus(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: r
        r = 0.5_dp + 0.6_dp*xi(1)
        x(1) = r*cos(xi(2))
        x(2) = r*sin(xi(2))
        if (size(x) == 3) then
            x(1) = (2.5_dp + r*cos(xi(2)))*cos(xi(3))
            x(2) = -(2.5_dp + r*cos(xi(2)))*sin(xi(3))
            x(3) = r*sin(xi(2)) + 0.1_dp*sin(xi(3))
        end if
    end subroutine

    subroutine mixed_oracle_agreement()
        ! Mixed products against the independent unfactored CSR oracle, with
        ! the exact diagonal and precomputed versus recomputed metric terms.
        type(tensor_line_tables_t) :: lines(3)
        type(tensor_mesh_t) :: mesh
        type(tensor_operator_t) :: op, stream
        type(csr_matrix_t) :: a
        real(dp), allocatable :: x(:), y(:), v(:), z(:), diagonal(:)
        integer :: case, status, i, j, dim
        real(dp) :: err
        do case = 1, 4
            select case (case)
            case (1) ! spline r x Fourier theta (annulus)
                dim = 2
                call build_bspline_line(3, lines(1), status, nelem=3)
                call build_fourier_line(3, lines(2), status)
            case (2) ! GLL r x periodic C1 quadratic spline theta
                dim = 2
                call build_gll_line(4, 2, lines(1), status)
                call build_bspline_line(2, lines(2), status, periodic=.true., &
                                        knots=[(2*pi*i/5, i=0, 5)])
            case (3) ! spline r x periodic spline theta x Fourier phi (torus)
                dim = 3
                call build_bspline_line(2, lines(1), status, nelem=2)
                call build_bspline_line(3, lines(2), status, periodic=.true., &
                                        knots=[(2*pi*i/4, i=0, 4)])
                call build_fourier_line(2, lines(3), status)
            case default ! C0 cubic spline r x Fourier theta x GLL phi
                dim = 3
                call build_bspline_line(3, lines(1), status, nelem=2, continuity=0)
                call build_fourier_line(2, lines(2), status)
                call build_gll_line(2, 3, lines(3), status, periodic=.true., &
                                    x0=0.0_dp, x1=2*pi)
            end select
            if (status /= 0) error stop 'mixed lines'
            call line_tensor_mesh(dim, lines, mesh, status, map=annulus)
            if (status /= 0) error stop 'mixed mesh'
            call build_tensor_operator(mesh, anisotropic_tensor, op, status)
            if (status /= 0) error stop 'mixed operator'
            call build_tensor_operator(mesh, anisotropic_tensor, stream, status, &
                                       precompute=.false.)
            if (status /= 0) error stop 'mixed stream operator'
            allocate (x(mesh%ndof), y(mesh%ndof), v(mesh%ndof))
            do i = 1, mesh%ndof
                x(i) = sin(real(11*i, dp))
            end do
            call apply_tensor_operator(op, x, y)
            call apply_tensor_operator(stream, x, v)
            if (norm2(y - v)/norm2(y) > 1.e-13_dp) error stop 'mixed recompute'
            call assemble_tensor_operator(op, a)
            z = csr_times(a, x)
            err = norm2(y - z)/norm2(y)
            print '(a,i1,a,i7,a,es10.2)', 'mixed case ', case, ' dofs=', mesh%ndof, &
                ' MF/CSR=', err
            if (err > 1.e-12_dp) error stop 'mixed MF versus CSR oracle'
            call tensor_diagonal(op, diagonal)
            do i = 1, a%n
                do j = a%rowptr(i), a%rowptr(i + 1) - 1
                    if (a%colidx(j) /= i) cycle
               if (abs(diagonal(i) - a%val(j)) > 1.e-12_dp*max(1.0_dp, abs(a%val(j)))) &
                        error stop 'mixed exact diagonal'
                end do
            end do
            deallocate (x, y, v, z, diagonal)
        end do
    end subroutine

    subroutine anisotropic_tensor(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0.1_dp
        do i = 1, size(x)
            k(i, i) = 1 + 0.2_dp*x(i)**2
        end do
        c = 0.3_dp
    end subroutine

    subroutine shell(xi, x)
        ! Annulus (2D) or cylindrical shell (3D); phi is the last parameter.
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: r, phi
        r = 0.5_dp + 0.6_dp*xi(1)
        phi = xi(size(xi))
        x(1) = r*cos(phi)
        x(2) = r*sin(phi)
        if (size(x) == 3) x(3) = 1 - xi(2) ! right-handed (r, z, phi)
    end subroutine

    subroutine trig_parts(x, s, ds, d2s, z, d2z, t, d2t, r)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: s, ds, d2s, z, d2z, t, d2t, r
        real(dp) :: phi, xi
        r = norm2(x(1:2))
        phi = atan2(x(2), x(1))
        xi = (r - 0.5_dp)/0.6_dp
        s = xi*(1 - xi); ds = (1 - 2*xi)/0.6_dp; d2s = -2/0.36_dp
        t = 1 + cos(phi) + 0.5_dp*sin(2*phi); d2t = -cos(phi) - 2*sin(2*phi)
        z = 1; d2z = 0
        if (size(x) == 3) then
            z = x(3)*(1 - x(3)); d2z = -2
        end if
    end subroutine

    function trig_exact(x) result(u)
        real(dp), intent(in) :: x(:)
        real(dp) :: u, s, ds, d2s, z, d2z, t, d2t, r
        call trig_parts(x, s, ds, d2s, z, d2z, t, d2t, r)
        u = s*z*t
    end function

    function trig_forcing(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f, s, ds, d2s, z, d2z, t, d2t, r
        call trig_parts(x, s, ds, d2s, z, d2z, t, d2t, r)
        f = -((d2s + ds/r)*z*t + s*d2z*t + s*z*d2t/r**2)
    end function

    subroutine trig_polynomial_exactness()
        ! u = (radial and axial spline-space polynomials) x (trig polynomial) is
        ! reproduced to quadrature roundoff on an annulus / cylindrical shell,
        ! solved with the nested h- and mode-coarsening multigrid.
        type(tensor_line_tables_t) :: lines(3)
        type(tensor_mesh_t) :: meshes(2)
        type(tensor_operator_t) :: op
        type(csr_matrix_t) :: transfers(1)
        type(line_matrix_t) :: p1d(3)
        type(tensor_multigrid_t) :: mg
        type(tensor_solve_info_t) :: info
        real(dp), allocatable :: rhs(:), x(:)
        integer :: dim, level, status, m, n, d
        real(dp) :: error
        do dim = 2, 3
            do level = 1, 2
                n = 4/level
                m = 6/level
                call build_bspline_line(3, lines(1), status, nelem=n, nq=20)
                if (dim == 3) call build_bspline_line(3, lines(2), status, nelem=n)
                call build_fourier_line(m, lines(dim), status, nq=2*m + 3)
                if (status /= 0) error stop 'trig lines'
                call line_tensor_mesh(dim, lines, meshes(level), status, map=shell)
                if (status /= 0) error stop 'trig mesh'
            end do
            call bspline_line_prolongation(3, 2, p1d(1)%a, status)
            if (status /= 0) error stop 'spline prolongation'
            if (dim == 3) p1d(2)%a = p1d(1)%a
            call fourier_line_prolongation(6, 3, p1d(dim)%a, status)
            if (status /= 0) error stop 'Fourier prolongation'
            if (dim == 2) then
                allocate (p1d(3)%a(1, 1))
                p1d(3)%a = 1
            end if
            call line_mesh_prolongation(meshes(1), meshes(2), p1d, &
                                    [2, merge(1, 2, dim == 2), 1], transfers(1), status)
            if (status /= 0) error stop 'mesh prolongation'
            call build_tensor_operator(meshes(1), laplace_coefficient, op, status)
            if (status /= 0) error stop 'trig operator'
            call tensor_load(op, trig_forcing, rhs)
           call build_line_multigrid(meshes, transfers, laplace_coefficient, mg, status)
            if (status /= 0) error stop 'line multigrid'
            allocate (x(meshes(1)%ndof))
            x = 0
            call solve_tensor_pcg(mg, rhs, 1.e-13_dp, x, info, maxiter=100)
            if (info%status /= 0) error stop 'trig PCG'
            error = tensor_l2_error(op, x, trig_exact)
            print '(a,i1,a,i7,a,i4,a,es10.2)', 'spline x Fourier exactness dim=', dim, &
                ' dofs=', meshes(1)%ndof, ' its=', info%iterations, ' L2=', error
            if (error > 1.e-11_dp) error stop 'spline x Fourier exactness'
            deallocate (x, rhs)
            do d = 1, 3
                if (allocated(p1d(d)%a)) deallocate (p1d(d)%a)
            end do
        end do
    end subroutine

    subroutine dense(a, m)
        type(csr_matrix_t), intent(in) :: a
        real(dp), allocatable, intent(out) :: m(:, :)
        integer :: i, k
        allocate (m(a%n, a%n))
        m = 0
        do i = 1, a%n
            do k = a%rowptr(i), a%rowptr(i + 1) - 1
                m(i, a%colidx(k)) = a%val(k)
            end do
        end do
    end subroutine

    subroutine torus_map(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x(1) = (3 + xi(1))*cos(xi(3))
        x(2) = -(3 + xi(1))*sin(xi(3))
        x(3) = 1.2_dp*xi(2) + 0.1_dp*sin(xi(3))*xi(1)
    end subroutine

    subroutine galerkin_nestedness()
        ! Nested hierarchies: interior block of P^T A_f P equals A_c, through the
        ! conforming multipatch O-grid interfaces and Fourier mode truncation.
        type(tensor_line_tables_t) :: lines(3)
        type(tensor_mesh_t) :: fine, coarse
        type(tensor_operator_t) :: opf, opc
        type(csr_matrix_t) :: af, ac, p
        type(line_matrix_t) :: p1d(3)
        real(dp), allocatable :: mf(:, :), mc(:, :), pd(:, :), galerkin(:, :)
        integer :: case, status, i, j, k, dim, deg
        real(dp) :: err
        do case = 1, 3
            dim = merge(2, 3, case < 3)
            deg = merge(2, 3, case == 1)
            ! Degree-14 geometry and 16/48-point rules make both levels' geometry
            ! and quadrature agree to roundoff, so the identity isolates P.
            do k = 1, 2
                call build_bspline_line(deg, lines(1), status, nelem=4/k, nq=16)
                if (dim == 3) call build_fourier_line(4/k, lines(3), status, nq=48)
                if (status /= 0) error stop 'nested lines'
                if (k == 1) then
                    if (dim == 2) then
                        call line_ogrid_mesh(2, lines, fine, status, q=14)
                    else
                       call line_ogrid_mesh(3, lines, fine, status, map=torus_map, q=14)
                    end if
                else
                    if (dim == 2) then
                        call line_ogrid_mesh(2, lines, coarse, status, q=14)
                    else
                     call line_ogrid_mesh(3, lines, coarse, status, map=torus_map, q=14)
                    end if
                end if
                if (status /= 0) error stop 'nested O-grid mesh'
            end do
            call bspline_line_prolongation(deg, 2, p1d(1)%a, status)
            if (status /= 0) error stop 'nested spline prolongation'
            p1d(2)%a = p1d(1)%a
            if (dim == 3) then
                call fourier_line_prolongation(4, 2, p1d(3)%a, status)
            else
                allocate (p1d(3)%a(1, 1))
                p1d(3)%a = 1
            end if
            call line_mesh_prolongation(fine, coarse, p1d, [2, 2, 1], p, status)
            if (status /= 0) error stop 'nested mesh prolongation'
            call build_tensor_operator(fine, anisotropic_tensor, opf, status)
            if (status /= 0) error stop 'nested fine operator'
            call build_tensor_operator(coarse, anisotropic_tensor, opc, status)
            if (status /= 0) error stop 'nested coarse operator'
            call assemble_tensor_operator(opf, af)
            call assemble_tensor_operator(opc, ac)
            call dense(af, mf)
            call dense(ac, mc)
            allocate (pd(p%n, coarse%ndof))
            pd = 0
            do i = 1, p%n
                do j = p%rowptr(i), p%rowptr(i + 1) - 1
                    pd(i, p%colidx(j)) = p%val(j)
                end do
            end do
            galerkin = matmul(transpose(pd), matmul(mf, pd))
            err = 0
            do j = 1, coarse%ndof
                if (coarse%boundary(j)) cycle
                do i = 1, coarse%ndof
                    if (coarse%boundary(i)) cycle
                    err = max(err, abs(galerkin(i, j) - mc(i, j)))
                end do
            end do
            err = err/maxval(abs(mc))
            print '(a,i1,a,i6,a,i6,a,es10.2)', 'nested case ', case, ' fine=', &
                fine%ndof, ' coarse=', coarse%ndof, ' |P^T A P - A_c|=', err
            if (err > 1.e-12_dp) error stop 'Galerkin nestedness'
            ! Spline constants are reproduced exactly by the prolongation.
            if (dim == 2) then
                if (maxval(abs(matmul(pd, [(1.0_dp, i=1, coarse%ndof)]) - 1)) > &
                    1.e-13_dp) error stop 'prolongation partition of unity'
            end if
            deallocate (pd, mf, mc, galerkin)
            do i = 1, 3
                deallocate (p1d(i)%a)
            end do
        end do
    end subroutine
end program
