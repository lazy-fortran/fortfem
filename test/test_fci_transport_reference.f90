program test_fci_transport_reference
    use, intrinsic :: iso_fortran_env, only: int64
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    use fortfem_sparse_direct, only: sparse_direct_factor_t, &
        sparse_direct_factor_csc, sparse_direct_solve_factored, sparse_direct_free
    use fortsparse, only: sparse_factor_stats, sparse_factor_stats_t
    use fortfem_fci_parallel_operator, only: assemble_fci_parallel_gradient_csc
    use fortfem_fci_transport_reference, only: fci_transport_reference_t, &
        setup_fci_transport_reference, solve_fci_transport_reference
    use fortsparse, only: csc_t, csc_from_triplet, fortsparse_status_t
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none

    character(32) :: mode

    call get_command_argument(1, mode)
    if (len_trim(mode) > 0) then
        call cost_comparison(trim(mode))
        call check_summary('FCI reference cost')
        stop
    end if
    call periodic_ladder()
    call mixed_terminal_ladder()
    call interpolated_kernel()
    call invalid_inputs()
    call check_summary('FCI transport reference')
contains
    subroutine from_dense(a, csc, status)
        real(dp), intent(in) :: a(:, :)
        type(csc_t), intent(out) :: csc
        type(fortsparse_status_t), intent(out) :: status
        integer :: rows(count(abs(a) > 0.0_dp)), columns(count(abs(a) > 0.0_dp))
        real(dp) :: values(count(abs(a) > 0.0_dp))
        integer :: i, j, k

        k = 0
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                if (abs(a(i, j)) <= 0.0_dp) cycle
                k = k + 1
                rows(k) = i
                columns(k) = j
                values(k) = a(i, j)
            end do
        end do
        call csc_from_triplet(size(a, 1), size(a, 2), rows, columns, &
            values, csc, status)
    end subroutine from_dense

    subroutine cost_comparison(kind)
        character(*), intent(in) :: kind
        type(fci_transport_reference_t) :: solver
        type(sparse_direct_factor_t) :: factor
        type(sparse_factor_stats_t) :: stats
        type(csc_t) :: gradient, perpendicular, assembled
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: q(:, :), p(:, :), mass(:), weights(:)
        real(dp), allocatable :: exact(:), mode(:), rhs(:), solution(:), a(:, :)
        real(dp) :: eps, h, eigenvalue, setup_seconds, solve_seconds
        real(dp) :: t, error, residual
        integer :: refinement, planes, n, i, j, exponent, repeat, solve_status
        integer, parameter :: repeats = 10
        integer(int64) :: retained_bytes

        call check_condition(kind == 'reference' .or. kind == 'direct', &
            'known comparison method')
        if (kind /= 'reference' .and. kind /= 'direct') return
        write (*, '(a)') &
            'method,n,threads,epsilon,setup_s,solve_s,rel_mass_error,retained_bytes'
        do refinement = 1, 3
            planes = 16*2**(refinement - 1)
            n = 3*planes
            h = 1.0_dp/planes
            call mapped_gradient(planes, .false., q, weights)
            call plane_operator(planes, .false., p, mass)
            call from_dense(q, gradient, status)
            call from_dense(p, perpendicular, status)
            allocate (exact(n), mode(n), rhs(n), solution(n))
            do j = 1, planes
                do i = 1, 3
                    mode(3*(j - 1) + i) = real(i, dp)* &
                        cos(2.0_dp*acos(-1.0_dp)*real(j - 1, dp)/planes)
                end do
            end do
            eigenvalue = 4.0_dp*sin(acos(-1.0_dp)/planes)**2/h
            do exponent = 6, 13, 7
                eps = 10.0_dp**(-exponent)
                exact = 1.0_dp + eps*mode
                rhs = matmul(p, exact) + eigenvalue*mode
                t = omp_get_wtime()
                if (kind == 'reference') then
                    call setup_fci_transport_reference(gradient, weights, mass, &
                        perpendicular, [integer ::], 3, solver, status)
                    call check_condition(status%code == 0, 'reference cost setup')
                    retained_bytes = solver%storage_bytes()
                else
                    a = q
                    do i = 1, size(q, 1)
                        a(i, :) = weights(i)*a(i, :)
                    end do
                    a = p + matmul(transpose(q), a)/eps
                    call from_dense(a, assembled, status)
                    call sparse_direct_factor_csc(factor, n, assembled%col_ptr, &
                        assembled%row_idx, assembled%val, solve_status)
                    call check_condition(solve_status == 0, 'direct cost factor')
                    if (solve_status /= 0) return
                    call sparse_factor_stats(factor, stats, status)
                    retained_bytes = -1
                    if (stats%factor_storage_bytes >= 0 .and. &
                        stats%workspace_bytes >= 0) retained_bytes = &
                        stats%factor_storage_bytes + stats%workspace_bytes
                end if
                setup_seconds = omp_get_wtime() - t
                t = omp_get_wtime()
                do repeat = 1, repeats
                    if (kind == 'reference') then
                        call solve_fci_transport_reference(solver, 0.0_dp, &
                            1.0_dp/eps, 1.0_dp, rhs, [real(dp) ::], &
                            solution, residual, status)
                        call check_condition(status%code == 0, 'reference cost solve')
                    else
                        call sparse_direct_solve_factored(factor, rhs, &
                            solution, solve_status)
                        call check_condition(solve_status == 0, 'direct cost solve')
                    end if
                end do
                solve_seconds = (omp_get_wtime() - t)/repeats
                error = sqrt(sum(mass*(solution - exact)**2)/sum(mass*exact**2))
                if (kind == 'reference') call check_condition(error < 2.0e-11_dp, &
                    'reference cost observable accuracy')
                write (*, '(a,2(",",i0),4(",",es14.6),",",i0)') &
                    kind, n, omp_get_max_threads(), eps, setup_seconds, &
                    solve_seconds, error, retained_bytes
                if (kind == 'direct') call sparse_direct_free(factor)
            end do
            deallocate (exact, mode, rhs, solution, q, p, mass, weights)
        end do
    end subroutine cost_comparison

    subroutine periodic_ladder()
        type(fci_transport_reference_t) :: solver
        type(csc_t) :: gradient, perpendicular
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: q(:, :), p(:, :), mass(:), weights(:)
        real(dp), allocatable :: exact(:), rhs(:), solution(:), mode(:)
        real(dp) :: eps, alpha, h, eigenvalue, error, residual
        integer :: refinement, planes, n, i, j, exponent

        do refinement = 1, 3
            planes = 4*2**(refinement - 1)
            n = 3*planes
            h = 1.0_dp/real(planes, dp)
            call mapped_gradient(planes, .false., q, weights)
            call plane_operator(planes, .false., p, mass)
            call from_dense(q, gradient, status)
            call from_dense(p, perpendicular, status)
            call setup_fci_transport_reference(gradient, weights, mass, &
                perpendicular, [integer ::], 3, solver, status)
            call check_condition(status%code == 0, 'periodic setup')
            if (status%code /= 0) cycle
            call check_condition(solver%nullity() == 3, 'three periodic averages')
            call check_condition(solver%kernel_defect() < 1.0e-12_dp, &
                'FCI kernel satisfies the original support gradient')
            allocate (exact(n), rhs(n), solution(n), mode(n))
            do j = 1, planes
                do i = 1, 3
                    mode(3*(j - 1) + i) = real(i, dp)* &
                        cos(2.0_dp*acos(-1.0_dp)*real(j - 1, dp)/planes)
                end do
            end do
            eigenvalue = 4.0_dp*sin(acos(-1.0_dp)/planes)**2/h
            do exponent = 6, 13
                eps = 10.0_dp**(-exponent)
                exact = 1.0_dp + eps*mode
                do i = 1, 2
                    alpha = real(i - 1, dp)*0.7_dp
                    ! The Fourier eigenvalue is independent of FCI assembly.
                    rhs = alpha*mass*exact + matmul(p, exact) + &
                        eigenvalue*mode
                    call solve_fci_transport_reference(solver, alpha, &
                        1.0_dp/eps, 1.0_dp, rhs, [real(dp) ::], &
                        solution, residual, status)
                    call check_condition(status%code == 0, 'periodic solve')
                    error = sqrt(sum(mass*(solution - exact)**2)/ &
                        sum(mass*exact**2))
                    call check_condition(error < 2.0e-11_dp, &
                        'periodic relative mass norm is uniform in epsilon')
                    call check_condition(residual < 2.0e-12_dp, &
                        'periodic block-scaled residual')
                    call solve_fci_transport_reference(solver, alpha*eps, &
                        1.0_dp, eps, rhs*eps, [real(dp) ::], &
                        solution, residual, status)
                    call check_condition(status%code == 0, 'rescaled periodic solve')
                    call check_condition(maxval(abs(solution - exact)) < 2.0e-11_dp, &
                        'coefficient normalization does not change the solution')
                    if (i == 1) write (*, '(a,i4,1x,es10.2,2(1x,es12.4))') &
                        'periodic ', n, eps, error, residual
                end do
            end do
            deallocate (exact, rhs, solution, mode, q, p, mass, weights)
        end do
    end subroutine periodic_ladder

    subroutine mixed_terminal_ladder()
        type(fci_transport_reference_t) :: solver
        type(csc_t) :: gradient, perpendicular
        type(fortsparse_status_t) :: status
        real(dp), allocatable :: q(:, :), p(:, :), mass(:), weights(:)
        real(dp), allocatable :: exact(:), rhs(:), solution(:), mode(:), base(:)
        real(dp) :: eps, h, error, residual, angle, closed_eigen, open_eigen
        integer :: planes, n, i, j, exponent, refinement

        do refinement = 1, 3
            planes = 4*2**(refinement - 1)
            n = 3*planes + 1
            h = 1.0_dp/planes
            call mapped_gradient(planes, .true., q, weights)
            call plane_operator(planes, .true., p, mass)
            call from_dense(q, gradient, status)
            call from_dense(p, perpendicular, status)
            call setup_fci_transport_reference(gradient, weights, mass, &
                perpendicular, [3, n], 2, solver, status)
            call check_condition(status%code == 0, 'mixed closed/open setup')
            if (status%code /= 0) cycle
            call check_condition(solver%nullity() == 2, &
                'anchored open line contributes no kernel')
            allocate (exact(n), rhs(n), solution(n), mode(n), base(n))
            do j = 1, planes
                angle = acos(-1.0_dp)*real(j - 1, dp)/planes
                do i = 1, 2
                    base(3*(j - 1) + i) = real(i, dp)
                    mode(3*(j - 1) + i) = real(i, dp)*cos(2.0_dp*angle)
                end do
                base(3*j) = 2.0_dp + 3.0_dp*real(j - 1, dp)/planes
                mode(3*j) = sin(angle)
            end do
            base(n) = 5.0_dp
            mode(n) = 0.0_dp
            closed_eigen = 4.0_dp*sin(acos(-1.0_dp)/planes)**2/h
            open_eigen = 4.0_dp*sin(acos(-1.0_dp)/(2*planes))**2/h
            do exponent = 6, 13
                eps = 10.0_dp**(-exponent)
                exact = base + eps*mode
                rhs = matmul(p, exact)
                do j = 1, planes
                    rhs(3*j - 2:3*j - 1) = rhs(3*j - 2:3*j - 1) + &
                        closed_eigen*mode(3*j - 2:3*j - 1)
                    rhs(3*j) = rhs(3*j) + open_eigen*mode(3*j)
                end do
                call solve_fci_transport_reference(solver, 0.0_dp, &
                    1.0_dp/eps, 1.0_dp, rhs, [2.0_dp, 5.0_dp], &
                    solution, residual, status)
                call check_condition(status%code == 0, 'mixed solve')
                error = sqrt(sum(mass*(solution - exact)**2)/sum(mass*exact**2))
                call check_condition(error < 2.0e-11_dp, &
                    'mixed relative mass norm is uniform in epsilon')
                call check_condition(maxval(abs(solution([3, n]) - &
                    [2.0_dp, 5.0_dp])) < 1.0e-14_dp, 'prescribed plate values')
                write (*, '(a,i4,1x,es10.2,2(1x,es12.4))') &
                    'mixed    ', n, eps, error, residual
            end do
            deallocate (exact, rhs, solution, mode, base, q, p, mass, weights)
        end do
    end subroutine mixed_terminal_ladder

    subroutine mapped_gradient(planes, mixed, q, weights)
        integer, intent(in) :: planes
        logical, intent(in) :: mixed
        real(dp), allocatable, intent(out) :: q(:, :), weights(:)
        type(csc_t) :: raw
        type(fortsparse_status_t) :: status
        real(dp) :: forward(3, 3, planes), backward(3, 3, planes)
        real(dp) :: lengths(3, planes)
        integer :: i, column, destination, k, n

        forward = 0.0_dp
        do i = 1, 3
            forward(i, i, :) = 1.0_dp
        end do
        backward = forward
        lengths = 1.0_dp/planes
        call assemble_fci_parallel_gradient_csc( &
            forward, backward, lengths, raw, status)
        call check_condition(status%code == 0, 'existing FCI map consumer')
        n = 3*planes
        if (mixed) n = n + 1
        allocate (q(3*planes, n), weights(3*planes))
        q = 0.0_dp
        weights = 1.0_dp/planes
        do column = 1, raw%ncol
            destination = column
            if (column > 3*planes) then
                destination = column - 3*planes
                if (mixed .and. destination == 3) destination = n
            end if
            do k = raw%col_ptr(column), raw%col_ptr(column + 1) - 1
                q(raw%row_idx(k), destination) = &
                    q(raw%row_idx(k), destination) + raw%val(k)
            end do
        end do
    end subroutine mapped_gradient

    subroutine plane_operator(planes, mixed, p, mass)
        integer, intent(in) :: planes
        logical, intent(in) :: mixed
        real(dp), allocatable, intent(out) :: p(:, :), mass(:)
        integer :: n, i, j, row
        real(dp) :: h, factor

        h = 1.0_dp/planes
        n = 3*planes
        if (mixed) n = n + 1
        allocate (p(n, n), mass(n))
        p = 0.0_dp
        do j = 1, planes
            factor = h*(1.0_dp + 0.2_dp*cos(real(j, dp)))
            do i = 1, 3
                row = 3*(j - 1) + i
                mass(row) = h*(1.0_dp + 0.1_dp*i)
                p(row, row) = 2.0_dp*factor
                if (i == 3) cycle
                p(row, row + 1) = -factor
                p(row + 1, row) = -factor
            end do
        end do
        if (mixed) then
            p(n, n) = h
            mass(n) = h
        end if
    end subroutine plane_operator

    subroutine interpolated_kernel()
        type(fci_transport_reference_t) :: solver
        type(csc_t) :: raw, gradient, perpendicular
        type(fortsparse_status_t) :: status
        real(dp) :: forward(2, 2, 2), backward(2, 2, 2), q(4, 4)
        real(dp) :: lengths(2, 2), p(4, 4), rhs(4), solution(4), exact(4)
        real(dp) :: weights(4), mode(4), residual, eps
        integer :: j, k, target, exponent

        forward(:, :, 1) = reshape([0.8_dp, 0.3_dp, 0.2_dp, 0.7_dp], [2, 2])
        forward(:, :, 2) = forward(:, :, 1)
        backward = 0.0_dp
        backward(1, 1, :) = 1.0_dp
        backward(2, 2, :) = 1.0_dp
        lengths = 1.0_dp
        call assemble_fci_parallel_gradient_csc( &
            forward, backward, lengths, raw, status)
        q = 0.0_dp
        do j = 1, 6
            target = 1 + mod(j - 1, 4)
            do k = raw%col_ptr(j), raw%col_ptr(j + 1) - 1
                q(raw%row_idx(k), target) = q(raw%row_idx(k), target) + raw%val(k)
            end do
        end do
        weights = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
        p = 0.0_dp
        do j = 1, 4
            p(j, j) = real(j, dp)
        end do
        call from_dense(q, gradient, status)
        call from_dense(p, perpendicular, status)
        call setup_fci_transport_reference(gradient, weights, &
            [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], perpendicular, &
            [integer ::], 1, solver, status)
        call check_condition(status%code == 0, 'interpolating periodic map setup')
        if (status%code /= 0) return
        ! Explicit rows provide a second assembly path for the oracle.
        q(1, :) = [-1.0_dp, 0.0_dp, 0.8_dp, 0.2_dp]
        q(2, :) = [0.0_dp, -1.0_dp, 0.3_dp, 0.7_dp]
        q(3, :) = [0.8_dp, 0.2_dp, -1.0_dp, 0.0_dp]
        q(4, :) = [0.3_dp, 0.7_dp, 0.0_dp, -1.0_dp]
        mode = [-0.7_dp, 0.2_dp, 0.9_dp, -0.4_dp]
        do exponent = 6, 13
            eps = 10.0_dp**(-exponent)
            exact = 1.0_dp + eps*mode
            rhs = matmul(p, exact) + matmul(transpose(q), weights*matmul(q, mode))
            call solve_fci_transport_reference(solver, 0.0_dp, 1.0_dp/eps, &
                1.0_dp, rhs, [real(dp) ::], solution, residual, status)
            call check_condition(status%code == 0, 'interpolating map solve')
            call check_condition(maxval(abs(solution - exact)) < 2.0e-12_dp, &
                'off-node map preserves the constant mode at high anisotropy')
        end do
    end subroutine interpolated_kernel

    subroutine invalid_inputs()
        type(fci_transport_reference_t) :: solver
        type(csc_t) :: gradient, perpendicular
        type(fortsparse_status_t) :: status
        real(dp) :: residual, solution(2), nan

        call from_dense(reshape([-1.0_dp, 1.0_dp], [1, 2]), &
            gradient, status)
        call from_dense(reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], &
            [2, 2]), perpendicular, status)
        call setup_fci_transport_reference(gradient, [1.0_dp], [1.0_dp, 1.0_dp], &
            perpendicular, [integer ::], 0, solver, status)
        call check_condition(status%code /= 0, 'reject wrong nullity')
        call setup_fci_transport_reference(gradient, [1.0_dp], [1.0_dp, 1.0_dp], &
            perpendicular, [1, 1], 0, solver, status)
        call check_condition(status%code /= 0, 'reject duplicate anchors')
        call setup_fci_transport_reference(gradient, [-1.0_dp], [1.0_dp, 1.0_dp], &
            perpendicular, [integer ::], 1, solver, status)
        call check_condition(status%code /= 0, 'reject negative support weight')
        call setup_fci_transport_reference(gradient, [1.0_dp], [1.0_dp, 1.0_dp], &
            perpendicular, [integer ::], 1, solver, status)
        call check_condition(status%code == 0, 'valid setup after rejection')
        call solve_fci_transport_reference(solver, 0.0_dp, 1.0_dp, 0.0_dp, &
            [1.0_dp, 1.0_dp], [real(dp) ::], solution, residual, status)
        call check_condition(status%code /= 0, 'reject unregularized kernel')
        nan = ieee_value(0.0_dp, ieee_quiet_nan)
        call solve_fci_transport_reference(solver, nan, 1.0_dp, 1.0_dp, &
            [1.0_dp, 1.0_dp], [real(dp) ::], solution, residual, status)
        call check_condition(status%code /= 0, 'reject nonfinite coefficients')
        call setup_fci_transport_reference(gradient, [1.0_dp], [1.0_dp, 1.0_dp], &
            perpendicular, [1], 0, solver, status)
        call check_condition(status%code == 0, 'fully anchored kernel is empty')
        call solve_fci_transport_reference(solver, 0.0_dp, 1.0_dp, 0.0_dp, &
            [0.0_dp, 0.0_dp], [3.0_dp], solution, residual, status)
        call check_condition(status%code == 0, 'pure parallel anchored solve')
        call check_condition(maxval(abs(solution - 3.0_dp)) < 1.0e-13_dp, &
            'pure parallel constant boundary extension')
    end subroutine invalid_inputs
end program test_fci_transport_reference
