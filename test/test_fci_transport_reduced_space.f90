program test_fci_transport_reduced_space
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_fci_field_line_tracer, only: trace_fci_field_line_rk4
    use fortfem_fci_interpolation_map, only: &
        build_fci_bilinear_interpolation_map_2d
    use fortfem_fci_parallel_operator, only: assemble_fci_parallel_gradient_csc
    use fortfem_fci_transport_reference, only: fci_transport_reference_t, &
        setup_fci_transport_reference, &
        solve_fci_transport_reference
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t
    use omp_lib, only: omp_get_wtime
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp), dz = 0.5_dp
    character(16) :: mode
    integer :: level, cells

    interface
        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
                work, lwork, info)
            import :: dp
            character, intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
            integer, intent(out) :: info
        end subroutine dgesvd
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            import :: dp
            character, intent(in) :: side, uplo, transa, diag
            integer, intent(in) :: m, n, lda, ldb
            real(dp), intent(in) :: alpha, a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
        end subroutine dtrsm
    end interface

    call get_command_argument(1, mode)
    call norm_oracle()
    write (*, '(a)') 'space,cells,free_dofs,multiplier_dofs,surface_bins,'// &
        'kernel_dim,rank,beta_l2,beta_h1,flux_kernel_error,'// &
        'checkerboard_retained,projected_defect,full_gradient_defect,seconds'
    do level = 1, 3
        cells = 20*2**(level - 1) + 1
        if (trim(mode) /= 'ladder' .and. level > 1) exit
        call probe(cells, trim(mode) == 'gate')
    end do
    call check_summary('FCI reduced constraint spaces')
contains
    subroutine norm_oracle()
        real(dp) :: identity(16, 16), beta, expected
        integer :: i

        identity = 0.0_dp
        do i = 1, 16
            identity(i, i) = 1.0_dp
        end do
        call energy_singular(identity, 5, .true., beta)
        expected = 1.0_dp/sqrt(1.0_dp + 4.0_dp/dz**2 + &
            8.0_dp*5**2*cos(pi/10.0_dp)**2)
        call check_condition(abs(beta - expected) < 1.0e-14_dp, &
            'H1 scaling agrees with analytic highest sine eigenvalue')
    end subroutine norm_oracle

    pure subroutine field_line_rhs(z, point, derivative)
        real(dp), intent(in) :: z, point(:)
        real(dp), intent(out) :: derivative(:)
        real(dp) :: speed, bx, by

        ! Analytic fixture and oracle, independent of transport kernels.
        associate (unused => z)
            speed = 0.15_dp*(1.0_dp + sin(pi*point(1))*sin(pi*point(2)))
            bx = speed*pi*sin(pi*point(1))*cos(pi*point(2))
            by = -speed*pi*cos(pi*point(1))*sin(pi*point(2))
            derivative = [bx, by, sqrt(1.0_dp + bx*bx + by*by)]
        end associate
    end subroutine field_line_rhs

    subroutine traced_blocks(cells, plus, minus, grid, x, y)
        integer, intent(in) :: cells
        real(dp), allocatable, intent(out) :: plus(:, :), minus(:, :)
        real(dp), allocatable, intent(out) :: grid(:), x(:), y(:)
        type(csc_t) :: gradient
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: forward(:, :, :), backward(:, :, :), ell(:, :)
        real(dp), allocatable :: tx(:), ty(:)
        real(dp) :: point(3), endpoint(3), refined(3), trace_defect, step_defect
        integer :: p, side, n, i, j, row, col, entry, free_col, sign
        logical :: trace_ok

        p = cells - 1
        side = cells + 1
        n = p*p
        allocate (grid(side), x(n), y(n), tx(n), ty(n))
        grid = [(real(i - 1, dp)/cells, i=1, side)]
        allocate (forward(n, side*side, 1), backward(n, side*side, 1), ell(n, 1))
        trace_defect = 0.0_dp
        step_defect = 0.0_dp
        trace_ok = .true.
        do j = 1, p
            do i = 1, p
                row = (j - 1)*p + i
                x(row) = grid(i + 1)
                y(row) = grid(j + 1)
                point = [x(row), y(row), 0.0_dp]
                call trace_fci_field_line_rk4(point, dz, 64, &
                    field_line_rhs, endpoint, status)
                trace_ok = trace_ok .and. status%code == 0
                call trace_fci_field_line_rk4(point, dz, 128, &
                    field_line_rhs, refined, status)
                trace_ok = trace_ok .and. status%code == 0
                trace_defect = max(trace_defect, &
                    abs(sin(pi*endpoint(1))*sin(pi*endpoint(2)) - &
                    sin(pi*x(row))*sin(pi*y(row))))
                step_defect = max(step_defect, maxval(abs(endpoint - refined)))
                tx(row) = endpoint(1)
                ty(row) = endpoint(2)
                ell(row, 1) = endpoint(3)
            end do
        end do
        call check_condition(trace_ok, 'native tracing succeeds')
        call check_condition(max(trace_defect, step_defect) < 1.0e-9_dp, &
            'flux invariant and step doubling control tracing')
        call build_fci_bilinear_interpolation_map_2d(grid, grid, tx, ty, &
            forward(:, :, 1), status)
        call check_condition(status%code == 0, 'native traced endpoint map')
        call build_fci_bilinear_interpolation_map_2d(grid, grid, x, y, &
            backward(:, :, 1), status)
        call check_condition(status%code == 0, 'native departure map')
        call assemble_fci_parallel_gradient_csc(forward, backward, ell, &
            gradient, status)
        call check_condition(status%code == 0, 'native interval gradient')
        deallocate (forward, backward, ell)
        allocate (plus(n, n), minus(n, n))
        plus = 0.0_dp
        minus = 0.0_dp
        ! Identical intervals on two periodic planes diagonalize into +/- blocks.
        do col = 1, gradient%ncol
            i = modulo(col - 1, side) + 1
            j = modulo((col - 1)/side, side) + 1
            if (i == 1 .or. i == side .or. j == 1 .or. j == side) cycle
            free_col = (j - 2)*p + i - 1
            sign = 1
            if (col > side*side) sign = -1
            do entry = gradient%col_ptr(col), gradient%col_ptr(col + 1) - 1
                row = gradient%row_idx(entry)
                plus(row, free_col) = plus(row, free_col) + gradient%val(entry)
                minus(row, free_col) = minus(row, free_col) + &
                    sign*gradient%val(entry)
            end do
        end do
    end subroutine traced_blocks

    subroutine test_basis(candidate, cells, grid, x, y, basis)
        integer, intent(in) :: candidate, cells
        real(dp), intent(in) :: grid(:), x(:), y(:)
        real(dp), allocatable, intent(out) :: basis(:, :)
        real(dp), allocatable :: coarse(:)
        type(fortsparse_status_t) :: status
        integer :: p, nq, i, j, row, basis_column

        p = cells - 1
        if (candidate == 1) then
            nq = p/2
            allocate (basis(p*p, nq*nq), source=0.0_dp)
            do j = 1, p
                do i = 1, p
                    row = (j - 1)*p + i
                    basis_column = ((j - 1)/2)*nq + (i - 1)/2 + 1
                    basis(row, basis_column) = 1.0_dp
                end do
            end do
        else
            nq = p/2 + 1
            allocate (coarse(nq), basis(p*p, nq*nq))
            coarse(:nq - 1) = grid(2:cells - 1:2)
            coarse(nq) = grid(cells)
            call build_fci_bilinear_interpolation_map_2d(coarse, coarse, x, y, &
                basis, status)
            call check_condition(status%code == 0, 'native coarse Q1 test basis')
        end if
    end subroutine test_basis

    subroutine row_space(matrix, vectors, singular, rank)
        real(dp), intent(in) :: matrix(:, :)
        real(dp), allocatable, intent(out) :: vectors(:, :), singular(:)
        integer, intent(out) :: rank
        real(dp), allocatable :: copy(:, :), work(:)
        real(dp) :: unused(1, 1), query(1), tolerance
        integer :: m, n, info

        m = size(matrix, 1)
        n = size(matrix, 2)
        copy = matrix
        allocate (vectors(m, n), singular(m))
        call dgesvd('N', 'S', m, n, copy, m, singular, unused, 1, &
            vectors, m, query, -1, info)
        call check_condition(info == 0, 'LAPACK SVD workspace query')
        allocate (work(ceiling(query(1))))
        call dgesvd('N', 'S', m, n, copy, m, singular, unused, 1, &
            vectors, m, work, size(work), info)
        call check_condition(info == 0, 'dense row-space SVD')
        tolerance = 8.0_dp*epsilon(1.0_dp)*max(m, n)*maxval(singular)
        rank = count(singular > tolerance)
    end subroutine row_space

    subroutine energy_singular(matrix, cells, odd_plane, beta)
        real(dp), intent(in) :: matrix(:, :)
        integer, intent(in) :: cells
        logical, intent(in) :: odd_plane
        real(dp), intent(out) :: beta
        real(dp), allocatable :: sine(:, :), transformed(:, :)
        real(dp), allocatable :: vectors(:, :), singular(:)
        integer :: p, i, k, row, rank

        p = cells - 1
        allocate (sine(p, p), transformed(size(matrix, 1), p*p))
        do k = 1, p
            do i = 1, p
                sine(i, k) = sqrt(2.0_dp/cells)*sin(pi*i*k/cells)
            end do
        end do
        !$omp parallel do schedule(static)
        do row = 1, size(matrix, 1)
            call energy_row(matrix(row, :), sine, cells, odd_plane, &
                transformed(row, :))
        end do
        !$omp end parallel do
        call row_space(transformed, vectors, singular, rank)
        beta = 0.0_dp
        if (rank == size(matrix, 1)) beta = minval(singular)
    end subroutine energy_singular

    subroutine energy_row(row, sine, cells, odd_plane, output)
        real(dp), intent(in) :: row(:), sine(:, :)
        integer, intent(in) :: cells
        logical, intent(in) :: odd_plane
        real(dp), intent(out) :: output(:)
        real(dp) :: plane(cells - 1, cells - 1), value, shift
        integer :: i, j, p

        p = cells - 1
        plane = matmul(transpose(sine), matmul(reshape(row, [p, p]), sine))
        shift = 1.0_dp
        if (odd_plane) shift = shift + 4.0_dp/dz**2
        do j = 1, p
            do i = 1, p
                value = shift + 4.0_dp*cells**2*( &
                    sin(pi*i/(2*cells))**2 + sin(pi*j/(2*cells))**2)
                output((j - 1)*p + i) = plane(i, j)/sqrt(value)
            end do
        end do
    end subroutine energy_row

    subroutine probe(cells, enforce_gate)
        integer, intent(in) :: cells
        logical, intent(in) :: enforce_gate
        real(dp), allocatable :: plus(:, :), minus(:, :), basis(:, :)
        real(dp), allocatable :: gram(:, :), cp(:, :), cm(:, :), vt(:, :)
        real(dp), allocatable :: sigma(:), exact(:), checker(:), kp(:), km(:)
        real(dp), allocatable :: grid(:), x(:), y(:)
        real(dp) :: beta_l2, bp, bm, retained, flux_error, defect, full_defect, t
        integer :: candidate, p, n, nq, i, j, info, rp, rm
        character(8) :: name

        call traced_blocks(cells, plus, minus, grid, x, y)
        p = cells - 1
        n = p*p
        exact = sin(pi*x)*sin(pi*y)
        allocate (checker(n))
        do j = 1, p
            do i = 1, p
                checker((j - 1)*p + i) = (-1.0_dp)**(i + j)* &
                    exact((j - 1)*p + i)
            end do
        end do
        exact = exact/norm2(exact)
        checker = checker/norm2(checker)
        do candidate = 1, 2
            t = omp_get_wtime()
            name = 'P0-2h'
            if (candidate == 2) name = 'Q1-2h'
            call test_basis(candidate, cells, grid, x, y, basis)
            nq = size(basis, 2)
            gram = matmul(transpose(basis), basis)
            call dpotrf('L', nq, gram, nq, info)
            call check_condition(info == 0, 'positive multiplier mass matrix')
            cp = matmul(transpose(basis), plus)
            cm = matmul(transpose(basis), minus)
            ! Uniform support and nodal masses cancel in these normalized blocks.
            call dtrsm('L', 'L', 'N', 'N', nq, n, 1.0_dp, gram, nq, cp, nq)
            call dtrsm('L', 'L', 'N', 'N', nq, n, 1.0_dp, gram, nq, cm, nq)
            call row_space(cp, vt, sigma, rp)
            kp = exact - matmul(transpose(vt(:rp, :)), &
                matmul(vt(:rp, :), exact))
            beta_l2 = 0.0_dp
            if (rp == nq) beta_l2 = minval(sigma)
            call row_space(cm, vt, sigma, rm)
            km = checker - matmul(transpose(vt(:rm, :)), &
                matmul(vt(:rm, :), checker))
            if (rm == nq) then
                beta_l2 = min(beta_l2, minval(sigma))
            else
                beta_l2 = 0.0_dp
            end if
            call energy_singular(cp, cells, .false., bp)
            call energy_singular(cm, cells, .true., bm)
            retained = norm2(km)
            flux_error = norm2(kp - exact)
            defect = max(norm2(matmul(cp, kp)), norm2(matmul(cm, km)))
            full_defect = norm2(matmul(minus, km))/max(norm2(km), tiny(1.0_dp))
            if (cells == 21) call write_kernel(name, x, y, km)
            call check_condition(defect < 1.0e-10_dp, &
                'SVD projections satisfy reduced constraints')
            if (cells == 21) call original_oracle(cp, cm, rp + rm, &
                checker, km)
            write (*, '(a,6(",",i0),7(",",es16.8))') trim(name), cells, &
                2*n, 2*nq, p/2, 2*n - rp - rm, rp + rm, beta_l2, &
                min(bp, bm), flux_error, retained, defect, full_defect, &
                omp_get_wtime() - t
            if (enforce_gate) then
                call check_condition(retained < 1.0e-10_dp, &
                    'gate: no odd-plane checkerboard in flux kernel')
            end if
        end do
    end subroutine probe

    subroutine write_kernel(name, x, y, kernel)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(:), y(:), kernel(:)
        character(1024) :: directory
        integer :: unit, i, io_status

        call get_command_argument(2, directory)
        if (len_trim(directory) == 0) return
        open (newunit=unit, file=trim(directory)//'/'//trim(name)//'-kernel.csv', &
            status='replace', action='write', iostat=io_status)
        call check_condition(io_status == 0, 'open optional kernel artifact')
        if (io_status /= 0) return
        write (unit, '(a)') 'x,y,odd_plane_kernel'
        do i = 1, size(x)
            write (unit, '(es24.16,2(",",es24.16))') x(i), y(i), kernel(i)
        end do
        close (unit)
    end subroutine write_kernel

    subroutine original_oracle(cp, cm, rank, checker, projected)
        real(dp), intent(in) :: cp(:, :), cm(:, :), checker(:), projected(:)
        integer, intent(in) :: rank
        type(fci_transport_reference_t) :: solver
        type(fortsparse_status_t) :: status
        type(csc_t) :: gradient, perpendicular
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: values(:), mass(:), weights(:), load(:), u(:)
        real(dp) :: residual
        integer :: m, n, i, j, k

        m = size(cp, 1)
        n = size(cp, 2)
        allocate (rows(2*m*n), cols(2*m*n), values(2*m*n))
        k = 0
        do j = 1, n
            do i = 1, m
                k = k + 1
                rows(k) = i
                cols(k) = j
                values(k) = cp(i, j)
                k = k + 1
                rows(k) = m + i
                cols(k) = n + j
                values(k) = cm(i, j)
            end do
        end do
        call csc_from_triplet(2*m, 2*n, rows, cols, values, gradient, status)
        call check_condition(status%code == 0, 'projected oracle gradient')
        deallocate (rows, cols, values)
        rows = [(i, i=1, 2*n)]
        mass = spread(1.0_dp, 1, 2*n)
        weights = spread(1.0_dp, 1, 2*m)
        call csc_from_triplet(2*n, 2*n, rows, rows, mass, perpendicular, status)
        call check_condition(status%code == 0, 'oracle identity slow operator')
        call setup_fci_transport_reference(gradient, weights, mass, &
            perpendicular, [integer ::], 2*n - rank, &
            solver, status)
        call check_condition(status%code == 0, 'u3a oracle agrees on nullity')
        if (status%code /= 0) return
        allocate (load(2*n), u(2*n))
        load(:n) = 0.0_dp
        load(n + 1:) = checker
        call solve_fci_transport_reference(solver, 0.0_dp, 1.0e13_dp, 1.0_dp, &
            load, [real(dp) ::], u, residual, status)
        call check_condition(status%code == 0, 'u3a penalized mass projection')
        call check_condition(norm2(u(n + 1:) - projected) < 1.0e-10_dp, &
            'u3a oracle retains checkerboard at epsilon=1e-13')
    end subroutine original_oracle
end program test_fci_transport_reduced_space
