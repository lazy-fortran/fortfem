module fortfem_spline_energy
    !! Energies sum_points w W(J(c)) of pointwise kernels W over the jets J of
    !! multipatch spline fields, with gradient, assembled Hessian and a
    !! matrix-free Hessian action (doc/design/multipatch_spline_spaces.md,
    !! eqs. 4, 9 and 10).
    !!
    !! The kernel returns densities and their derivatives with respect to
    !! the jet. The library multiplies them by the quadrature weights.
    !! Kernels are called concurrently for different cells and must not
    !! modify shared state.
    use fortfem_kinds, only: dp
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_assemble_blocks
    use fortfem_spline_space, only: spline_space_t, spline_cell_broken
    use fortfem_spline_eval, only: spline_quadrature_t, SPLINE_JET_COUNT, &
        spline_broken_coefficients, spline_cell_tables, spline_cell_xi_jets
    implicit none
    private

    integer, parameter, public :: SPLINE_HESSIAN_DENSE = 0
    integer, parameter, public :: SPLINE_HESSIAN_PACKED = 1

    type, abstract, public :: pointwise_energy_t
        !! field_count fields, jets to order (0, 1 or 2). hessian_storage
        !! selects hessian(ncomp*ncomp, npoint), column-major, or the
        !! upper triangle packed column-major as in LAPACK 'U',
        !! hessian(ncomp*(ncomp + 1)/2, npoint).
        integer :: field_count = 1
        integer :: order = 1
        integer :: hessian_storage = SPLINE_HESSIAN_DENSE
    contains
        procedure(energy_batch), deferred :: evaluate
    end type pointwise_energy_t

    abstract interface
        subroutine energy_batch(self, cell, first_point, npoint, jet, weight, &
                need, w, gradient, hessian, status)
            !! need = 0: w; 1: w and gradient; 2: also hessian. Points
            !! first_point .. first_point + npoint - 1 of the quadrature.
            import :: pointwise_energy_t, dp
            class(pointwise_energy_t), intent(in) :: self
            integer, intent(in) :: cell, first_point, npoint, need
            real(dp), intent(in) :: jet(:, :), weight(:)
            real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
            integer, intent(out) :: status
        end subroutine energy_batch
    end interface

    public :: spline_energy
    public :: spline_energy_hessian_csr
    public :: spline_energy_hessian_apply

contains

    subroutine spline_energy(space, quad, kernel, coef, energy, status, gradient)
        !! Energy and, if present, gradient(field_count, dof_count).
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        class(pointwise_energy_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :)
        real(dp), intent(out) :: energy
        integer, intent(out) :: status
        real(dp), intent(out), optional :: gradient(:, :)

        real(dp), allocatable :: broken(:, :), cell_out(:, :, :), u(:, :)
        real(dp), allocatable :: jet(:, :), g(:, :), h(:, :), w(:), xi_jets(:, :, :)
        real(dp), allocatable :: table1(:, :, :), table2(:, :, :), gx(:, :, :)
        integer, allocatable :: cell_broken(:)
        real(dp) :: total
        integer :: nf, nd, nc, nq, nl, p, cell, q, need, first, ierr, kstat

        status = 1
        if (.not. valid_call(space, kernel, coef)) return
        nf = kernel%field_count
        nd = SPLINE_JET_COUNT(kernel%order)
        nc = nf*nd
        nq = quad%points_per_cell
        p = space%degree
        nl = (p + 1)**2
        need = merge(1, 0, present(gradient))
        allocate(broken(nf, space%broken_count), cell_out(nf, nl, space%cell_count))
        call spline_broken_coefficients(space, coef, broken)
        total = 0.0_dp
        ierr = 0
        !$omp parallel default(none) &
        !$omp shared(space, quad, kernel, broken, cell_out, nf, nd, nc, nq, nl, p, need) &
        !$omp private(u, jet, g, h, w, xi_jets, table1, table2, gx, cell_broken, &
        !$omp cell, q, first, kstat) reduction(+:total) reduction(max:ierr)
        allocate(u(nf, nl), jet(nc, nq), g(nc, nq), h(1, nq), w(nq), &
            xi_jets(nf, 6, nq), table1(p + 1, quad%points_per_direction, 0:2), &
            table2(p + 1, quad%points_per_direction, 0:2), gx(nf, 6, nq), &
            cell_broken(nl))
        !$omp do
        do cell = 1, space%cell_count
            first = quad%point_offset(cell)
            call cell_jets(space, quad, broken, cell, nd, cell_broken, u, &
                table1, table2, xi_jets, jet)
            call kernel%evaluate(cell, first, nq, jet, quad%weight(first:first + nq - 1), &
                need, w, g, h, kstat)
            ierr = max(ierr, kstat)
            total = total + dot_product(w, quad%weight(first:first + nq - 1))
            if (need == 1) then
                do q = 1, nq
                    g(:, q) = quad%weight(first + q - 1)*g(:, q)
                end do
                call jets_transpose(quad, first, nf, nd, g, gx)
                call cell_xi_jets_transpose(table1, table2, gx, u)
                cell_out(:, :, cell) = u
            end if
        end do
        !$omp end do
        !$omp end parallel
        energy = total
        status = 2
        if (ierr /= 0) return
        if (present(gradient)) call scatter(space, cell_out, gradient)
        status = 0
    end subroutine spline_energy

    subroutine spline_energy_hessian_apply(space, quad, kernel, coef, direction, &
            product, status)
        !! product = Hessian(coef) * direction, both (field_count, dof_count),
        !! without forming the matrix.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        class(pointwise_energy_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :), direction(:, :)
        real(dp), intent(out) :: product(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: broken(:, :), broken_v(:, :), cell_out(:, :, :)
        real(dp), allocatable :: u(:, :), jet(:, :), jet_v(:, :), g(:, :), h(:, :)
        real(dp), allocatable :: w(:), xi_jets(:, :, :), table1(:, :, :)
        real(dp), allocatable :: table2(:, :, :), gx(:, :, :), y(:, :)
        integer, allocatable :: cell_broken(:)
        integer :: nf, nd, nc, nq, nl, p, cell, q, first, ierr, kstat, hl

        status = 1
        if (.not. valid_call(space, kernel, coef)) return
        if (any(shape(direction) /= shape(coef))) return
        nf = kernel%field_count
        nd = SPLINE_JET_COUNT(kernel%order)
        nc = nf*nd
        nq = quad%points_per_cell
        p = space%degree
        nl = (p + 1)**2
        hl = hessian_length(kernel%hessian_storage, nc)
        allocate(broken(nf, space%broken_count), broken_v(nf, space%broken_count), &
            cell_out(nf, nl, space%cell_count))
        call spline_broken_coefficients(space, coef, broken)
        call spline_broken_coefficients(space, direction, broken_v)
        ierr = 0
        !$omp parallel default(none) &
        !$omp shared(space, quad, kernel, broken, broken_v, cell_out, nf, nd, nc, nq, &
        !$omp nl, p, hl) &
        !$omp private(u, jet, jet_v, g, h, w, xi_jets, table1, table2, gx, y, &
        !$omp cell_broken, cell, q, first, kstat) reduction(max:ierr)
        allocate(u(nf, nl), jet(nc, nq), jet_v(nc, nq), g(nc, nq), h(hl, nq), &
            w(nq), xi_jets(nf, 6, nq), table1(p + 1, quad%points_per_direction, 0:2), &
            table2(p + 1, quad%points_per_direction, 0:2), gx(nf, 6, nq), &
            y(nc, nq), cell_broken(nl))
        !$omp do
        do cell = 1, space%cell_count
            first = quad%point_offset(cell)
            call cell_jets(space, quad, broken, cell, nd, cell_broken, u, &
                table1, table2, xi_jets, jet)
            call kernel%evaluate(cell, first, nq, jet, quad%weight(first:first + nq - 1), &
                2, w, g, h, kstat)
            ierr = max(ierr, kstat)
            call cell_jets(space, quad, broken_v, cell, nd, cell_broken, u, &
                table1, table2, xi_jets, jet_v)
            do q = 1, nq
                call hessian_times(kernel%hessian_storage, nc, h(:, q), jet_v(:, q), &
                    y(:, q))
                y(:, q) = quad%weight(first + q - 1)*y(:, q)
            end do
            call jets_transpose(quad, first, nf, nd, y, gx)
            call cell_xi_jets_transpose(table1, table2, gx, u)
            cell_out(:, :, cell) = u
        end do
        !$omp end do
        !$omp end parallel
        status = 2
        if (ierr /= 0) return
        call scatter(space, cell_out, product)
        status = 0
    end subroutine spline_energy_hessian_apply

    subroutine spline_energy_hessian_csr(space, quad, kernel, coef, matrix, status)
        !! Assembled Hessian on the unknowns f + field_count*(I - 1). Status 3:
        !! the space has broken functions shared by several global functions
        !! (C1 extraction), which this assembly does not cover yet.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        class(pointwise_energy_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :)
        type(csr_matrix_t), intent(out) :: matrix
        integer, intent(out) :: status

        real(dp), allocatable :: broken(:, :), blocks(:, :, :), u(:, :), jet(:, :)
        real(dp), allocatable :: g(:, :), h(:, :), w(:), xi_jets(:, :, :)
        real(dp), allocatable :: table1(:, :, :), table2(:, :, :), phi(:, :)
        real(dp), allocatable :: hphi(:, :), basis(:, :, :), sgn(:)
        integer, allocatable :: conn(:, :), cell_broken(:)
        integer :: nf, nd, nc, nq, nl, p, cell, q, first, ierr, kstat, hl
        integer :: f, a, j, col, ds

        status = 1
        if (.not. valid_call(space, kernel, coef)) return
        status = 3
        if (any(space%broken_ptr(2:) - space%broken_ptr(:space%broken_count) /= 1)) &
            return
        nf = kernel%field_count
        nd = SPLINE_JET_COUNT(kernel%order)
        nc = nf*nd
        nq = quad%points_per_cell
        p = space%degree
        nl = (p + 1)**2
        hl = hessian_length(kernel%hessian_storage, nc)
        allocate(broken(nf, space%broken_count), conn(nf*nl, space%cell_count), &
            blocks(nf*nl, nf*nl, space%cell_count))
        call spline_broken_coefficients(space, coef, broken)
        ierr = 0
        !$omp parallel default(none) &
        !$omp shared(space, quad, kernel, broken, blocks, conn, nf, nd, nc, nq, nl, &
        !$omp p, hl) &
        !$omp private(u, jet, g, h, w, xi_jets, table1, table2, phi, hphi, basis, &
        !$omp sgn, cell_broken, cell, q, first, kstat, f, a, j, col, ds) &
        !$omp reduction(max:ierr)
        allocate(u(nf, nl), jet(nc, nq), g(nc, nq), h(hl, nq), w(nq), &
            xi_jets(nf, 6, nq), table1(p + 1, quad%points_per_direction, 0:2), &
            table2(p + 1, quad%points_per_direction, 0:2), phi(nc, nf*nl), &
            hphi(nc, nf*nl), basis(nd, nl, nq), sgn(nf*nl), cell_broken(nl))
        !$omp do
        do cell = 1, space%cell_count
            first = quad%point_offset(cell)
            call cell_jets(space, quad, broken, cell, nd, cell_broken, u, &
                table1, table2, xi_jets, jet)
            call kernel%evaluate(cell, first, nq, jet, quad%weight(first:first + nq - 1), &
                2, w, g, h, kstat)
            ierr = max(ierr, kstat)
            call basis_s_jets(quad, first, table1, table2, nd, basis)
            do a = 1, nl
                j = cell_broken(a)
                do f = 1, nf
                    col = f + nf*(a - 1)
                    conn(col, cell) = f + nf*(space%broken_dof(space%broken_ptr(j)) - 1)
                    sgn(col) = space%broken_weight(space%broken_ptr(j))
                end do
            end do
            blocks(:, :, cell) = 0.0_dp
            do q = 1, nq
                phi = 0.0_dp
                do a = 1, nl
                    do f = 1, nf
                        col = f + nf*(a - 1)
                        do ds = 1, nd
                            phi(f + nf*(ds - 1), col) = sgn(col)*basis(ds, a, q)
                        end do
                    end do
                end do
                do col = 1, nf*nl
                    call hessian_times(kernel%hessian_storage, nc, h(:, q), &
                        phi(:, col), hphi(:, col))
                end do
                blocks(:, :, cell) = blocks(:, :, cell) + &
                    quad%weight(first + q - 1)*matmul(transpose(phi), hphi)
            end do
        end do
        !$omp end do
        !$omp end parallel
        status = 2
        if (ierr /= 0) return
        call csr_build_pattern(nf*space%dof_count, conn, .false., matrix, &
            with_slots=.true.)
        call csr_assemble_blocks(matrix, conn, blocks)
        status = 0
    end subroutine spline_energy_hessian_csr

    logical function valid_call(space, kernel, coef) result(valid)
        type(spline_space_t), intent(in) :: space
        class(pointwise_energy_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :)

        valid = kernel%field_count >= 1 .and. kernel%order >= 0 .and. &
            kernel%order <= 2 .and. size(coef, 1) == kernel%field_count .and. &
            size(coef, 2) == space%dof_count .and. &
            (kernel%hessian_storage == SPLINE_HESSIAN_DENSE .or. &
            kernel%hessian_storage == SPLINE_HESSIAN_PACKED)
    end function valid_call

    pure integer function hessian_length(storage, n) result(length)
        integer, intent(in) :: storage, n

        if (storage == SPLINE_HESSIAN_PACKED) then
            length = n*(n + 1)/2
        else
            length = n*n
        end if
    end function hessian_length

    pure subroutine hessian_times(storage, n, h, x, y)
        integer, intent(in) :: storage, n
        real(dp), intent(in) :: h(:), x(:)
        real(dp), intent(out) :: y(:)

        integer :: i, j, k

        if (storage == SPLINE_HESSIAN_DENSE) then
            y = matmul(reshape(h, [n, n]), x(:n))
            return
        end if
        y(:n) = 0.0_dp
        k = 0
        do j = 1, n
            do i = 1, j
                k = k + 1
                y(i) = y(i) + h(k)*x(j)
                if (i /= j) y(j) = y(j) + h(k)*x(i)
            end do
        end do
    end subroutine hessian_times

    subroutine cell_jets(space, quad, broken, cell, nd, cell_broken, u, table1, &
            table2, xi_jets, jet)
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        real(dp), intent(in) :: broken(:, :)
        integer, intent(in) :: cell, nd
        integer, intent(out) :: cell_broken(:)
        real(dp), intent(out) :: u(:, :), table1(:, :, 0:), table2(:, :, 0:)
        real(dp), intent(out) :: xi_jets(:, :, :), jet(:, :)

        real(dp) :: h(2)
        integer :: nf, q, d, c, first

        nf = size(u, 1)
        call spline_cell_broken(space, cell, cell_broken)
        u = broken(:, cell_broken)
        call spline_cell_tables(space, quad, cell, table1, table2, h)
        call spline_cell_xi_jets(table1, table2, u, xi_jets)
        first = quad%point_offset(cell)
        do q = 1, quad%points_per_cell
            do d = 1, nd
                jet(nf*(d - 1) + 1:nf*d, q) = 0.0_dp
                do c = 1, 6
                    if (quad%chain(d, c, first + q - 1) == 0.0_dp) cycle
                    jet(nf*(d - 1) + 1:nf*d, q) = jet(nf*(d - 1) + 1:nf*d, q) + &
                        quad%chain(d, c, first + q - 1)*xi_jets(:, c, q)
                end do
            end do
        end do
    end subroutine cell_jets

    pure subroutine jets_transpose(quad, first, nf, nd, g, gx)
        !! gx(f, c, q) = sum_d chain(d, c) g(f + nf*(d - 1), q).
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: first, nf, nd
        real(dp), intent(in) :: g(:, :)
        real(dp), intent(out) :: gx(:, :, :)

        integer :: q, d, c

        gx = 0.0_dp
        do q = 1, size(g, 2)
            do c = 1, 6
                do d = 1, nd
                    if (quad%chain(d, c, first + q - 1) == 0.0_dp) cycle
                    gx(:, c, q) = gx(:, c, q) + quad%chain(d, c, first + q - 1)* &
                        g(nf*(d - 1) + 1:nf*d, q)
                end do
            end do
        end do
    end subroutine jets_transpose

    pure subroutine cell_xi_jets_transpose(table1, table2, gx, u)
        !! Adjoint of spline_cell_xi_jets: u(f, a + n1*(b - 1)) =
        !! sum_{q, c} N_a^(d1(c))(q1) N_b^(d2(c))(q2) gx(f, c, q).
        real(dp), intent(in) :: table1(:, :, 0:), table2(:, :, 0:), gx(:, :, :)
        real(dp), intent(out) :: u(:, :)

        integer, parameter :: d1(6) = [0, 1, 0, 2, 1, 0], d2(6) = [0, 0, 1, 0, 1, 2]
        real(dp) :: partial(size(u, 1), size(table1, 2), size(table2, 1), 0:2)
        integer :: n1, n2, nq1, nq2, a, b, q1, q2, c

        n1 = size(table1, 1)
        n2 = size(table2, 1)
        nq1 = size(table1, 2)
        nq2 = size(table2, 2)
        partial = 0.0_dp
        do q2 = 1, nq2
            do q1 = 1, nq1
                do c = 1, 6
                    do b = 1, n2
                        partial(:, q1, b, d1(c)) = partial(:, q1, b, d1(c)) + &
                            table2(b, q2, d2(c))*gx(:, c, q1 + nq1*(q2 - 1))
                    end do
                end do
            end do
        end do
        u = 0.0_dp
        do b = 1, n2
            do a = 1, n1
                do c = 0, 2
                    do q1 = 1, nq1
                        u(:, a + n1*(b - 1)) = u(:, a + n1*(b - 1)) + &
                            table1(a, q1, c)*partial(:, q1, b, c)
                    end do
                end do
            end do
        end do
    end subroutine cell_xi_jets_transpose

    pure subroutine basis_s_jets(quad, first, table1, table2, nd, basis)
        !! basis(d, a, q): s-jet component d of cell function a at point q.
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: first, nd
        real(dp), intent(in) :: table1(:, :, 0:), table2(:, :, 0:)
        real(dp), intent(out) :: basis(:, :, :)

        integer, parameter :: d1(6) = [0, 1, 0, 2, 1, 0], d2(6) = [0, 0, 1, 0, 1, 2]
        real(dp) :: xi(6)
        integer :: n1, nq1, a, b, q1, q2, q, c

        n1 = size(table1, 1)
        nq1 = size(table1, 2)
        do q2 = 1, size(table2, 2)
            do q1 = 1, nq1
                q = q1 + nq1*(q2 - 1)
                do b = 1, size(table2, 1)
                    do a = 1, n1
                        do c = 1, 6
                            xi(c) = table1(a, q1, d1(c))*table2(b, q2, d2(c))
                        end do
                        basis(:, a + n1*(b - 1), q) = &
                            matmul(quad%chain(:nd, :, first + q - 1), xi)
                    end do
                end do
            end do
        end do
    end subroutine basis_s_jets

    subroutine scatter(space, cell_out, out)
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: cell_out(:, :, :)
        real(dp), intent(out) :: out(:, :)

        integer, allocatable :: cell_broken(:)
        integer :: cell, a, j, q

        allocate(cell_broken(size(cell_out, 2)))
        out = 0.0_dp
        do cell = 1, space%cell_count
            call spline_cell_broken(space, cell, cell_broken)
            do a = 1, size(cell_broken)
                j = cell_broken(a)
                do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                    out(:, space%broken_dof(q)) = out(:, space%broken_dof(q)) + &
                        space%broken_weight(q)*cell_out(:, a, cell)
                end do
            end do
        end do
    end subroutine scatter

end module fortfem_spline_energy
