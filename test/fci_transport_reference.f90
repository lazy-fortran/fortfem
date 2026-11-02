module fortfem_fci_transport_reference
    !! Dense, kernel-separated oracle for small FCI transport systems.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortsparse, only: csc_t, csc_is_valid, status_set, fortsparse_status_t, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: fci_transport_reference_t
        private
        integer :: kernel_size = 0, full_size = 0
        real(dp) :: defect = 0.0_dp
        integer, allocatable :: free_dofs(:), fixed_dofs(:)
        real(dp), allocatable :: basis(:, :), sigma(:)
        real(dp), allocatable :: mass_reduced(:, :), perp_reduced(:, :)
        real(dp), allocatable :: perp_boundary(:, :), line_boundary(:, :)
    contains
        procedure :: nullity
        procedure :: kernel_defect
        procedure :: storage_bytes
        procedure, private :: copy_reference
        generic :: assignment(=) => copy_reference
    end type fci_transport_reference_t

    public :: setup_fci_transport_reference, solve_fci_transport_reference

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
    end interface
contains
    subroutine setup_fci_transport_reference(gradient, support_weights, mass, &
            perpendicular, fixed_dofs, expected_nullity, solver, status)
        type(csc_t), intent(in) :: gradient, perpendicular
        real(dp), intent(in) :: support_weights(:), mass(:)
        integer, intent(in) :: fixed_dofs(:), expected_nullity
        type(fci_transport_reference_t), intent(out) :: solver
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: q(:, :), p(:, :), a(:, :), singular(:)
        real(dp), allocatable :: u(:, :), vt(:, :), work(:), product(:, :)
        real(dp) :: query(1), rank_tolerance
        logical, allocatable :: fixed(:)
        integer :: m, n, nf, nb, rank, k, i, j, info

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'FCI reference: incompatible operators or weights')
        if (.not. csc_is_valid(gradient)) return
        if (.not. csc_is_valid(perpendicular)) return
        m = gradient%nrow
        n = gradient%ncol
        nb = size(fixed_dofs)
        nf = n - nb
        if (m < 1 .or. nf < 1) return
        if (perpendicular%nrow /= n .or. perpendicular%ncol /= n) return
        if (size(support_weights) /= m .or. size(mass) /= n) return
        if (expected_nullity < 0 .or. expected_nullity > nf) return
        if (any(.not. ieee_is_finite(gradient%val))) return
        if (any(.not. ieee_is_finite(perpendicular%val))) return
        if (any(.not. ieee_is_finite(support_weights))) return
        if (any(.not. ieee_is_finite(mass))) return
        if (any(support_weights <= 0.0_dp) .or. any(mass <= 0.0_dp)) return
        allocate (fixed(n))
        fixed = .false.
        do i = 1, nb
            j = fixed_dofs(i)
            if (j < 1 .or. j > n) return
            if (fixed(j)) return
            fixed(j) = .true.
        end do
        solver%free_dofs = pack([(i, i=1, n)],.not. fixed)
        solver%fixed_dofs = fixed_dofs
        allocate (q(m, n), p(n, n))
        q = 0.0_dp
        p = 0.0_dp
        do j = 1, n
            do k = gradient%col_ptr(j), gradient%col_ptr(j + 1) - 1
                i = gradient%row_idx(k)
                q(i, j) = q(i, j) + sqrt(support_weights(i))*gradient%val(k)
            end do
            do k = perpendicular%col_ptr(j), perpendicular%col_ptr(j + 1) - 1
                i = perpendicular%row_idx(k)
                p(i, j) = p(i, j) + perpendicular%val(k)
            end do
        end do
        if (any(.not. ieee_is_finite(q))) return
        if (any(.not. ieee_is_finite(p))) return
        if (any(abs(p - transpose(p)) > 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'FCI reference: perpendicular weak matrix must be symmetric')
            return
        end if
        a = q(:, solver%free_dofs)
        do j = 1, nf
            a(:, j) = a(:, j)/sqrt(mass(solver%free_dofs(j)))
        end do
        if (any(.not. ieee_is_finite(a))) return
        allocate (singular(min(m, nf)), u(m, min(m, nf)), vt(nf, nf))
        call dgesvd('S', 'A', m, nf, a, m, singular, u, m, vt, nf, &
            query, -1, info)
        if (info /= 0) return
        allocate (work(max(1, ceiling(query(1)))))
        call dgesvd('S', 'A', m, nf, a, m, singular, u, m, vt, nf, &
            work, size(work), info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'FCI reference: gradient SVD failed')
            return
        end if
        rank_tolerance = 8.0_dp*epsilon(1.0_dp)*max(m, nf)*maxval(singular)
        rank = count(singular > rank_tolerance)
        if (nf - rank /= expected_nullity) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'FCI reference: supplied nullity disagrees with support gradient')
            return
        end if
        solver%kernel_size = expected_nullity
        solver%sigma = singular(:rank)
        allocate (solver%basis(nf, nf))
        solver%basis(:, :expected_nullity) = transpose(vt(rank + 1:, :))
        solver%basis(:, expected_nullity + 1:) = transpose(vt(:rank, :))
        do i = 1, nf
            solver%basis(i, :) = solver%basis(i, :)/ &
                sqrt(mass(solver%free_dofs(i)))
        end do
        if (expected_nullity > 0) then
            product = matmul(q(:, solver%free_dofs), &
                solver%basis(:, :expected_nullity))
            solver%defect = maxval(abs(product))
        end if
        if (allocated(product)) deallocate (product)
        allocate (product(nf, nf), source=0.0_dp)
        ! Column independence keeps the projection work outside serial LAPACK.
        !$omp parallel do private(i) schedule(static)
        do j = 1, nf
            do i = 1, nf
                product(i, j) = mass(solver%free_dofs(i))*solver%basis(i, j)
            end do
        end do
        !$omp end parallel do
        solver%mass_reduced = matmul(transpose(solver%basis), product)
        product = matmul(p(solver%free_dofs, solver%free_dofs), solver%basis)
        solver%perp_reduced = matmul(transpose(solver%basis), product)
        solver%perp_boundary = matmul(transpose(solver%basis), &
            p(solver%free_dofs, fixed_dofs))
        solver%line_boundary = matmul(transpose(u(:, :rank)), &
            q(:, fixed_dofs))
        solver%full_size = n
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine setup_fci_transport_reference

    subroutine solve_fci_transport_reference(solver, alpha, chi_parallel, &
            chi_perpendicular, source, boundary_values, solution, &
            relative_scaled_residual, status)
        type(fci_transport_reference_t), intent(in) :: solver
        real(dp), intent(in) :: alpha, chi_parallel, chi_perpendicular
        real(dp), intent(in) :: source(:), boundary_values(:)
        real(dp), intent(out) :: solution(:), relative_scaled_residual
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: system(:, :), rhs(:), coefficients(:), image(:)
        real(dp) :: slow_scale, row_scale, residual_norm, source_norm
        integer :: n, nz, j, i, info

        solution = 0.0_dp
        relative_scaled_residual = huge(1.0_dp)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            'FCI reference: invalid solve inputs or missing setup')
        if (solver%full_size < 1) return
        if (size(source) /= solver%full_size) return
        if (size(solution) /= solver%full_size) return
        if (size(boundary_values) /= size(solver%fixed_dofs)) return
        if (any(.not. ieee_is_finite(source))) return
        if (any(.not. ieee_is_finite(boundary_values))) return
        if (.not. ieee_is_finite(alpha)) return
        if (.not. ieee_is_finite(chi_parallel)) return
        if (.not. ieee_is_finite(chi_perpendicular)) return
        if (alpha < 0.0_dp .or. chi_parallel <= 0.0_dp) return
        if (chi_perpendicular < 0.0_dp) return
        n = size(solver%free_dofs)
        nz = solver%kernel_size
        slow_scale = max(alpha, chi_perpendicular)
        if (nz > 0 .and. slow_scale <= 0.0_dp) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'FCI reference: field-line averages need a slow operator')
            return
        end if
        allocate (system(n, n), rhs(n), coefficients(n), image(n))
        rhs = matmul(transpose(solver%basis), source(solver%free_dofs))
        ! Scale kernel rows before the large parallel coefficient enters.
        !$omp parallel do private(i, row_scale) schedule(static)
        do j = 1, n
            do i = 1, n
                row_scale = chi_parallel
                if (i <= nz) row_scale = slow_scale
                system(i, j) = (alpha/row_scale)*solver%mass_reduced(i, j) + &
                    (chi_perpendicular/row_scale)*solver%perp_reduced(i, j)
            end do
        end do
        !$omp end parallel do
        do i = 1, n
            row_scale = chi_parallel
            if (i <= nz) row_scale = slow_scale
            rhs(i) = rhs(i)/row_scale - (chi_perpendicular/row_scale)* &
                dot_product(solver%perp_boundary(i, :), boundary_values)
        end do
        do i = nz + 1, n
            system(i, i) = system(i, i) + solver%sigma(i - nz)**2
            rhs(i) = rhs(i) - solver%sigma(i - nz)* &
                dot_product(solver%line_boundary(i - nz, :), boundary_values)
        end do
        if (any(.not. ieee_is_finite(system))) return
        if (any(.not. ieee_is_finite(rhs))) return
        call dense_solve(system, rhs, coefficients, info)
        if (info /= 0) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                'FCI reference: reduced transport solve failed')
            return
        end if
        if (any(.not. ieee_is_finite(coefficients))) return
        solution(solver%free_dofs) = matmul(solver%basis, coefficients)
        solution(solver%fixed_dofs) = boundary_values
        if (any(.not. ieee_is_finite(solution))) return
        image = matmul(system, coefficients) - rhs
        relative_scaled_residual = 0.0_dp
        ! Check both blocks; a small parallel residual cannot mask a bad average.
        if (nz > 0) then
            residual_norm = norm2(image(:nz))
            source_norm = norm2(rhs(:nz))
            if (source_norm > 0.0_dp) residual_norm = residual_norm/source_norm
            relative_scaled_residual = residual_norm
        end if
        if (nz < n) then
            residual_norm = norm2(image(nz + 1:))
            source_norm = norm2(rhs(nz + 1:))
            if (source_norm > 0.0_dp) residual_norm = residual_norm/source_norm
            relative_scaled_residual = max(relative_scaled_residual, residual_norm)
        end if
        call status_set(status, FORTSPARSE_OK, '')
    end subroutine solve_fci_transport_reference

    pure integer function nullity(self)
        class(fci_transport_reference_t), intent(in) :: self
        nullity = self%kernel_size
    end function nullity

    pure real(dp) function kernel_defect(self)
        class(fci_transport_reference_t), intent(in) :: self
        kernel_defect = self%defect
    end function kernel_defect

    pure integer(int64) function storage_bytes(self) result(bytes)
        class(fci_transport_reference_t), intent(in) :: self
        integer(int64) :: real_bytes, integer_bytes

        bytes = 0_int64
        if (self%full_size < 1) return
        real_bytes = storage_size(0.0_dp)/8
        integer_bytes = storage_size(0)/8
        bytes = integer_bytes*(size(self%free_dofs, kind=int64) + &
            size(self%fixed_dofs, kind=int64)) + real_bytes*( &
            size(self%basis, kind=int64) + size(self%sigma, kind=int64) + &
            size(self%mass_reduced, kind=int64) + &
            size(self%perp_reduced, kind=int64) + &
            size(self%perp_boundary, kind=int64) + &
            size(self%line_boundary, kind=int64))
    end function storage_bytes

    subroutine copy_reference(lhs, rhs)
        class(fci_transport_reference_t), intent(out) :: lhs
        class(fci_transport_reference_t), intent(in) :: rhs
        lhs%kernel_size = rhs%kernel_size
        lhs%full_size = rhs%full_size
        lhs%defect = rhs%defect
        if (allocated(rhs%free_dofs)) lhs%free_dofs = rhs%free_dofs
        if (allocated(rhs%fixed_dofs)) lhs%fixed_dofs = rhs%fixed_dofs
        if (allocated(rhs%basis)) lhs%basis = rhs%basis
        if (allocated(rhs%sigma)) lhs%sigma = rhs%sigma
        if (allocated(rhs%mass_reduced)) lhs%mass_reduced = rhs%mass_reduced
        if (allocated(rhs%perp_reduced)) lhs%perp_reduced = rhs%perp_reduced
        if (allocated(rhs%perp_boundary)) lhs%perp_boundary = rhs%perp_boundary
        if (allocated(rhs%line_boundary)) lhs%line_boundary = rhs%line_boundary
    end subroutine copy_reference
end module fortfem_fci_transport_reference
