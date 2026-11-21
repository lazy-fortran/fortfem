module fortfem_spline_toroidal
    !! Fields u_f(s, t) = sum_m sum_I c(f, m, I) T_m(t) phi_I(s) on a spline
    !! space in the section and a periodic basis T_m in t in [0, 1)
    !! (doc/design/multipatch_spline_spaces.md, eq. 3), with energies summed
    !! over section quadrature points and toroidal planes (eq. 4).
    !!
    !! Coefficients are coef(nfield, nmode, dof_count), the global index of
    !! (f, m, I) is f + nfield (m - 1) + nfield nmode (I - 1). A jet holds per
    !! field the value, the derivatives in (s1, s2, t) and, at order 2, the
    !! second derivatives (11, 21, 31, 22, 32, 33): component f + nfield d.
    !! The mode and plane contraction runs inside the energy calls; the
    !! kernel sees one cell and one plane at a time.
    use fortfem_kinds, only: dp
    use fortfem_gll_line, only: gauss_rule
    use fortfem_interval_bernstein_jets, only: interval_bernstein_jet
    use fortfem_spline_extraction_1d, only: build_periodic_extraction_1d
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_spline_space, only: spline_space_t, spline_cell_broken
    use fortfem_spline_eval, only: spline_quadrature_t, spline_cell_tables, &
        spline_cell_xi_jets, spline_evaluate_at
    implicit none
    private

    integer, parameter, public :: TOROIDAL_FOURIER = 0
    integer, parameter, public :: TOROIDAL_SPLINE = 1
    integer, parameter, public :: SPLINE_JET3_COUNT(0:2) = [1, 4, 10]
    ! 3D component d -> section s-jet component and t-derivative order
    integer, parameter :: c2_of(10) = [1, 2, 3, 1, 4, 5, 2, 6, 3, 1]
    integer, parameter :: dt_of(10) = [0, 0, 0, 1, 0, 0, 1, 0, 1, 2]
    integer, parameter :: section_count(0:2) = [1, 3, 6]

    type, public :: spline_toroidal_t
        integer :: kind = TOROIDAL_FOURIER
        integer :: mode_count = 0
        integer :: plane_count = 0
        integer :: degree = 0 ! periodic splines
        real(dp), allocatable :: t(:), weight(:) ! (plane), weights sum to 1
        real(dp), allocatable :: table(:, :, :) ! (mode, plane, 0:2)
        real(dp), allocatable :: breaks(:), extraction(:, :, :)
    end type spline_toroidal_t

    type, abstract, public :: pointwise_energy3_t
        !! field_count fields, jets to order 0, 1 or 2; hessian is dense,
        !! hessian(ncomp*ncomp, npoint) column-major.
        integer :: field_count = 1
        integer :: order = 1
    contains
        procedure(energy3_batch), deferred :: evaluate
    end type pointwise_energy3_t

    abstract interface
        subroutine energy3_batch(self, cell, plane, first_point, npoint, jet, &
                weight, need, w, gradient, hessian, status)
            !! Points first_point .. first_point + npoint - 1 of the section
            !! quadrature on plane `plane`; weight includes the plane weight.
            import :: pointwise_energy3_t, dp
            class(pointwise_energy3_t), intent(in) :: self
            integer, intent(in) :: cell, plane, first_point, npoint, need
            real(dp), intent(in) :: jet(:, :), weight(:)
            real(dp), intent(out) :: w(:), gradient(:, :), hessian(:, :)
            integer, intent(out) :: status
        end subroutine energy3_batch
    end interface

    public :: spline_toroidal_fourier
    public :: spline_toroidal_periodic
    public :: spline_toroidal_values
    public :: spline_energy3
    public :: spline_energy3_hessian_apply
    public :: spline_energy3_hessian_csr
    public :: spline_evaluate3
    public :: spline_evaluate3_at

contains

    subroutine spline_toroidal_fourier(nmax, plane_count, tor, status)
        !! Modes 1, cos(2 pi m t), sin(2 pi m t), m = 1..nmax, on
        !! plane_count >= 2 nmax + 1 equispaced planes. Status 1: invalid.
        integer, intent(in) :: nmax, plane_count
        type(spline_toroidal_t), intent(out) :: tor
        integer, intent(out) :: status

        integer :: j

        status = 1
        if (nmax < 0 .or. plane_count < 2*nmax + 1) return
        tor%kind = TOROIDAL_FOURIER
        tor%mode_count = 2*nmax + 1
        tor%plane_count = plane_count
        allocate(tor%t(plane_count), tor%weight(plane_count), &
            tor%table(tor%mode_count, plane_count, 0:2))
        do j = 1, plane_count
            tor%t(j) = real(j - 1, dp)/plane_count
            tor%weight(j) = 1.0_dp/plane_count
            call spline_toroidal_values(tor, tor%t(j), tor%table(:, j, :))
        end do
        status = 0
    end subroutine spline_toroidal_fourier

    subroutine spline_toroidal_periodic(degree, element_count, points_per_element, &
            tor, status)
        !! Periodic splines of degree `degree` on element_count uniform
        !! elements, Gauss planes. Status 1: invalid.
        integer, intent(in) :: degree, element_count, points_per_element
        type(spline_toroidal_t), intent(out) :: tor
        integer, intent(out) :: status

        real(dp), allocatable :: knots(:), x(:), w(:)
        integer, allocatable :: first(:)
        integer :: e, q, j

        status = 1
        if (degree < 1 .or. element_count < degree + 1 .or. points_per_element < 1) &
            return
        tor%kind = TOROIDAL_SPLINE
        tor%degree = degree
        tor%mode_count = element_count
        tor%breaks = [(real(e, dp)/element_count, e = 0, element_count)]
        call build_periodic_extraction_1d(degree, tor%breaks, knots, tor%extraction, &
            first, status)
        if (status /= 0) return
        call gauss_rule(points_per_element, x, w, status)
        if (status /= 0) return
        tor%plane_count = element_count*points_per_element
        allocate(tor%t(tor%plane_count), tor%weight(tor%plane_count), &
            tor%table(tor%mode_count, tor%plane_count, 0:2))
        j = 0
        do e = 1, element_count
            do q = 1, points_per_element
                j = j + 1
                tor%t(j) = tor%breaks(e) + 0.5_dp*(x(q) + 1.0_dp)/element_count
                tor%weight(j) = 0.5_dp*w(q)/element_count
                call spline_toroidal_values(tor, tor%t(j), tor%table(:, j, :))
            end do
        end do
        status = 0
    end subroutine spline_toroidal_periodic

    subroutine spline_toroidal_values(tor, t, values)
        !! values(m, d) = d^d T_m/dt^d at t (taken modulo 1).
        type(spline_toroidal_t), intent(in) :: tor
        real(dp), intent(in) :: t
        real(dp), intent(out) :: values(:, 0:)

        real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)
        real(dp) :: x, b(tor%degree + 1, 0:2), h, local(tor%degree + 1)
        integer :: m, e, a, i, d, ierr, n

        values = 0.0_dp
        x = modulo(t, 1.0_dp)
        if (tor%kind == TOROIDAL_FOURIER) then
            values(1, 0) = 1.0_dp
            do m = 1, (tor%mode_count - 1)/2
                values(2*m, 0) = cos(two_pi*m*x)
                values(2*m, 1) = -two_pi*m*sin(two_pi*m*x)
                values(2*m, 2) = -(two_pi*m)**2*cos(two_pi*m*x)
                values(2*m + 1, 0) = sin(two_pi*m*x)
                values(2*m + 1, 1) = two_pi*m*cos(two_pi*m*x)
                values(2*m + 1, 2) = -(two_pi*m)**2*sin(two_pi*m*x)
            end do
            return
        end if
        n = tor%mode_count
        e = min(int(x*n) + 1, n)
        h = 1.0_dp/n
        call interval_bernstein_jet(tor%degree, (x - tor%breaks(e))/h, b(:, 0), &
            b(:, 1), b(:, 2), ierr)
        do d = 0, 2
            local = matmul(tor%extraction(:, :, e), b(:, d))/h**d
            do a = 1, tor%degree + 1
                i = modulo(e + a - 2, n) + 1
                values(i, d) = values(i, d) + local(a)
            end do
        end do
    end subroutine spline_toroidal_values

    subroutine cell_basis(space, quad, cell, nd2, phi, broken)
        !! phi(a, c, q): section s-jet component c of cell function a.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        integer, intent(in) :: cell, nd2
        real(dp), intent(out) :: phi(:, :, :)
        integer, intent(out) :: broken(:)

        real(dp), allocatable :: table1(:, :, :), table2(:, :, :), eye(:, :), xj(:, :, :)
        real(dp) :: h(2)
        integer :: p, nq, nloc, a, q, c, first

        p = space%degree
        nq = quad%points_per_cell
        nloc = (p + 1)**2
        allocate(table1(p + 1, quad%points_per_direction, 0:2), &
            table2(p + 1, quad%points_per_direction, 0:2), eye(nloc, nloc), &
            xj(nloc, 6, nq))
        call spline_cell_tables(space, quad, cell, table1, table2, h)
        eye = 0.0_dp
        do a = 1, nloc
            eye(a, a) = 1.0_dp
        end do
        call spline_cell_xi_jets(table1, table2, eye, xj)
        first = quad%point_offset(cell)
        do q = 1, nq
            do c = 1, nd2
                phi(:, c, q) = matmul(xj(:, :, q), quad%chain(c, :, first + q - 1))
            end do
        end do
        call spline_cell_broken(space, cell, broken)
    end subroutine cell_basis

    subroutine plane_jets(tor, j, nf, nc3, u, jet)
        !! jet(f + nf d, q) on plane j from section jets u(f + nf (m-1), c, q).
        type(spline_toroidal_t), intent(in) :: tor
        integer, intent(in) :: j, nf, nc3
        real(dp), intent(in) :: u(:, :, :)
        real(dp), intent(out) :: jet(:, :)

        integer :: d, m, f, q

        jet = 0.0_dp
        do q = 1, size(u, 3)
            do d = 1, nc3/nf
                do m = 1, tor%mode_count
                    if (tor%table(m, j, dt_of(d)) == 0.0_dp) cycle
                    do f = 1, nf
                        jet(f + nf*(d - 1), q) = jet(f + nf*(d - 1), q) + &
                            tor%table(m, j, dt_of(d))*u(f + nf*(m - 1), c2_of(d), q)
                    end do
                end do
            end do
        end do
    end subroutine plane_jets

    subroutine plane_jets_transpose(tor, j, nf, nc3, g, weight, gu)
        !! gu(f + nf (m-1), c, q) += sum_d T_m^(dt(d)) g(f + nf d, q) weight(q).
        type(spline_toroidal_t), intent(in) :: tor
        integer, intent(in) :: j, nf, nc3
        real(dp), intent(in) :: g(:, :), weight(:)
        real(dp), intent(inout) :: gu(:, :, :)

        integer :: d, m, f, q

        do q = 1, size(gu, 3)
            do d = 1, nc3/nf
                do m = 1, tor%mode_count
                    if (tor%table(m, j, dt_of(d)) == 0.0_dp) cycle
                    do f = 1, nf
                        gu(f + nf*(m - 1), c2_of(d), q) = gu(f + nf*(m - 1), &
                            c2_of(d), q) + tor%table(m, j, dt_of(d))* &
                            g(f + nf*(d - 1), q)*weight(q)
                    end do
                end do
            end do
        end do
    end subroutine plane_jets_transpose

    subroutine broken_of(space, coef, nb, broken)
        !! Broken coefficients (nb, broken) of coef(nb, dof).
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: coef(:, :)
        integer, intent(in) :: nb
        real(dp), intent(out) :: broken(:, :)

        integer :: j, q

        do j = 1, space%broken_count
            broken(:nb, j) = 0.0_dp
            do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                broken(:nb, j) = broken(:nb, j) + &
                    space%broken_weight(q)*coef(:nb, space%broken_dof(q))
            end do
        end do
    end subroutine broken_of

    subroutine gather_global(space, gb, nb, g)
        !! g = E gb for broken values gb(nb, broken).
        type(spline_space_t), intent(in) :: space
        real(dp), intent(in) :: gb(:, :)
        integer, intent(in) :: nb
        real(dp), intent(out) :: g(:, :)

        integer :: j, q

        g = 0.0_dp
        do j = 1, space%broken_count
            do q = space%broken_ptr(j), space%broken_ptr(j + 1) - 1
                g(:nb, space%broken_dof(q)) = g(:nb, space%broken_dof(q)) + &
                    space%broken_weight(q)*gb(:nb, j)
            end do
        end do
    end subroutine gather_global

    subroutine spline_energy3(space, quad, tor, kernel, coef, energy, status, gradient)
        !! Energy and, if present, gradient(nfield, nmode, dof_count).
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(spline_toroidal_t), intent(in) :: tor
        class(pointwise_energy3_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :, :)
        real(dp), intent(out) :: energy
        integer, intent(out) :: status
        real(dp), intent(out), optional :: gradient(:, :, :)

        if (present(gradient)) then
            call energy3_core(space, quad, tor, kernel, coef, 1, energy, status, &
                gradient)
        else
            call energy3_core(space, quad, tor, kernel, coef, 0, energy, status)
        end if
    end subroutine spline_energy3

    subroutine spline_energy3_hessian_apply(space, quad, tor, kernel, coef, &
            direction, product, status)
        !! product = H(coef) direction, all (nfield, nmode, dof_count).
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(spline_toroidal_t), intent(in) :: tor
        class(pointwise_energy3_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :, :), direction(:, :, :)
        real(dp), intent(out) :: product(:, :, :)
        integer, intent(out) :: status

        real(dp) :: energy

        call energy3_core(space, quad, tor, kernel, coef, 2, energy, status, product, &
            direction)
    end subroutine spline_energy3_hessian_apply

    subroutine energy3_core(space, quad, tor, kernel, coef, mode, energy, status, &
            out, direction)
        !! mode 0: energy, 1: gradient, 2: Hessian times direction.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(spline_toroidal_t), intent(in) :: tor
        class(pointwise_energy3_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :, :)
        integer, intent(in) :: mode
        real(dp), intent(out) :: energy
        integer, intent(out) :: status
        real(dp), intent(out), optional :: out(:, :, :)
        real(dp), intent(in), optional :: direction(:, :, :)

        real(dp), allocatable :: bc(:, :), bd(:, :), gb(:, :), phi(:, :, :), u(:, :, :)
        real(dp), allocatable :: v(:, :, :), gu(:, :, :), jet(:, :), vjet(:, :)
        real(dp), allocatable :: g(:, :), h(:, :), w(:), wt(:), hv(:, :), go(:, :)
        integer, allocatable :: broken(:)
        integer :: nf, nm, nb, nc3, nd2, nq, nloc, cell, j, q, first, a, need, k, ierr

        status = 1
        nf = kernel%field_count
        nm = tor%mode_count
        nb = nf*nm
        if (kernel%order < 0 .or. kernel%order > 2) return
        if (size(coef, 1) /= nf .or. size(coef, 2) /= nm .or. &
            size(coef, 3) /= space%dof_count) return
        nc3 = nf*SPLINE_JET3_COUNT(kernel%order)
        nd2 = section_count(kernel%order)
        nq = quad%points_per_cell
        nloc = (space%degree + 1)**2
        need = merge(0, merge(1, 2, mode == 1), mode == 0)
        allocate(bc(nb, space%broken_count))
        call broken_of(space, reshape(coef, [nb, space%dof_count]), nb, bc)
        allocate(gb(nb, merge(space%broken_count, 0, mode >= 1)))
        gb = 0.0_dp
        if (mode == 2) then
            allocate(bd(nb, space%broken_count))
            call broken_of(space, reshape(direction, [nb, space%dof_count]), nb, bd)
        end if
        energy = 0.0_dp
        ierr = 0
        !$omp parallel default(shared) private(cell, first, phi, broken, u, jet, g, h, &
        !$omp& w, wt, gu, v, vjet, hv, j, q, a, k, status) reduction(+:energy, gb)
        allocate(phi(nloc, nd2, nq), broken(nloc), u(nb, nd2, nq), jet(nc3, nq), &
            g(nc3, nq), h(nc3*nc3, nq), w(nq), wt(nq), gu(nb, nd2, nq), &
            v(nb, nd2, nq), vjet(nc3, nq), hv(nc3, nq))
        !$omp do schedule(dynamic)
        do cell = 1, space%cell_count
            if (ierr /= 0) cycle
            first = quad%point_offset(cell)
            call cell_basis(space, quad, cell, nd2, phi, broken)
            call section_jets(bc(:, broken), phi, u)
            if (mode == 2) call section_jets(bd(:, broken), phi, v)
            if (mode >= 1) gu = 0.0_dp
            do j = 1, tor%plane_count
                wt = quad%weight(first:first + nq - 1)*tor%weight(j)
                call plane_jets(tor, j, nf, nc3, u, jet)
                call kernel%evaluate(cell, j, first, nq, jet, wt, need, w, g, h, status)
                if (status /= 0) then
                    !$omp atomic write
                    ierr = status
                end if
                energy = energy + sum(wt*w)
                if (mode == 1) then
                    call plane_jets_transpose(tor, j, nf, nc3, g, wt, gu)
                else if (mode == 2) then
                    call plane_jets(tor, j, nf, nc3, v, vjet)
                    do q = 1, nq
                        hv(:, q) = matmul(reshape(h(:, q), [nc3, nc3]), vjet(:, q))
                    end do
                    call plane_jets_transpose(tor, j, nf, nc3, hv, wt, gu)
                end if
            end do
            if (mode >= 1) then
                do a = 1, nloc
                    do q = 1, nq
                        do k = 1, nd2
                            gb(:, broken(a)) = gb(:, broken(a)) + gu(:, k, q)*phi(a, k, q)
                        end do
                    end do
                end do
            end if
        end do
        !$omp end do
        !$omp end parallel
        if (ierr /= 0) then
            status = ierr
            return
        end if
        if (mode >= 1) then
            allocate(go(nb, space%dof_count))
            call gather_global(space, gb, nb, go)
            out = reshape(go, [nf, nm, space%dof_count])
        end if
        status = 0
    end subroutine energy3_core

    pure subroutine section_jets(uloc, phi, u)
        !! u(i, c, q) = sum_a uloc(i, a) phi(a, c, q).
        real(dp), intent(in) :: uloc(:, :), phi(:, :, :)
        real(dp), intent(out) :: u(:, :, :)

        integer :: q

        do q = 1, size(phi, 3)
            u(:, :, q) = matmul(uloc, phi(:, :, q))
        end do
    end subroutine section_jets

    subroutine spline_energy3_hessian_csr(space, quad, tor, kernel, coef, matrix, &
            status)
        !! Assembled Hessian over the global index f + nf (m-1) + nf nm (I-1).
        !! Only mode pairs that share a plane are stored (all pairs for
        !! Fourier modes, a band of width 2q + 1 for periodic splines). Cells
        !! are processed in parallel; the kernel must be thread safe.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(spline_toroidal_t), intent(in) :: tor
        class(pointwise_energy3_t), intent(in) :: kernel
        real(dp), intent(in) :: coef(:, :, :)
        type(csr_matrix_t), intent(out) :: matrix
        integer, intent(out) :: status

        real(dp), allocatable :: bc(:, :), phi(:, :, :), u(:, :, :), jet(:, :)
        real(dp), allocatable :: g(:, :), w(:), m4(:, :, :, :), hall(:, :, :)
        real(dp), allocatable :: wall(:, :), lp(:, :, :), blocks(:, :)
        integer, allocatable :: broken(:), rowptr(:), col(:), pos(:)
        integer, allocatable :: pptr(:), pcol(:)
        logical, allocatable :: couple(:, :)
        integer :: nf, nm, nb, nc3, nd2, nq, nloc, cell, j, q, first, a, a2, d, d2
        integer :: f, f2, m, m2, i, i2, c, c2, pa, pb, ia, ib, r, k, nnz, np, pp
        integer :: ierr
        real(dp) :: tw, wa, wb, wab

        status = 1
        nf = kernel%field_count
        nm = tor%mode_count
        nb = nf*nm
        if (kernel%order < 0 .or. kernel%order > 2) return
        if (size(coef, 1) /= nf .or. size(coef, 2) /= nm .or. &
            size(coef, 3) /= space%dof_count) return
        nc3 = nf*SPLINE_JET3_COUNT(kernel%order)
        nd2 = section_count(kernel%order)
        nq = quad%points_per_cell
        nloc = (space%degree + 1)**2
        call section_pattern(space, rowptr, col)
        nnz = size(col)
        ! coupled (i, i2) pairs: modes that are both nonzero on some plane
        allocate(couple(nm, nm))
        couple = .false.
        do j = 1, tor%plane_count
            do m2 = 1, nm
                if (all(tor%table(m2, j, :) == 0.0_dp)) cycle
                do m = 1, nm
                    if (any(tor%table(m, j, :) /= 0.0_dp)) couple(m, m2) = .true.
                end do
            end do
        end do
        allocate(pptr(nb + 1))
        pptr(1) = 1
        do i = 1, nb
            m = (i - 1)/nf + 1
            pptr(i + 1) = pptr(i) + nf*count(couple(m, :))
        end do
        np = pptr(nb + 1) - 1
        allocate(pcol(np))
        do i = 1, nb
            m = (i - 1)/nf + 1
            k = pptr(i)
            do i2 = 1, nb
                if (.not. couple(m, (i2 - 1)/nf + 1)) cycle
                pcol(k) = i2
                k = k + 1
            end do
        end do
        allocate(blocks(np, nnz), bc(nb, space%broken_count))
        blocks = 0.0_dp
        call broken_of(space, reshape(coef, [nb, space%dof_count]), nb, bc)
        ierr = 0
        !$omp parallel default(shared) private(cell, first, phi, broken, u, jet, g, &
        !$omp& w, hall, wall, m4, lp, j, q, d, d2, f, f2, c, c2, tw, m, m2, wa, wb, &
        !$omp& i, pp, a, a2, pa, pb, ia, ib, r, wab, pos, status)
        allocate(phi(nloc, nd2, nq), broken(nloc), u(nb, nd2, nq), jet(nc3, nq), &
            g(nc3, nq), w(nq), m4(nb, nd2, nb, nd2), lp(nloc, nloc, np), &
            hall(nc3*nc3, nq, tor%plane_count), wall(nq, tor%plane_count), &
            pos(space%dof_count))
        pos = 0
        !$omp do schedule(dynamic)
        do cell = 1, space%cell_count
            if (ierr /= 0) cycle
            first = quad%point_offset(cell)
            call cell_basis(space, quad, cell, nd2, phi, broken)
            call section_jets(bc(:, broken), phi, u)
            do j = 1, tor%plane_count
                wall(:, j) = quad%weight(first:first + nq - 1)*tor%weight(j)
                call plane_jets(tor, j, nf, nc3, u, jet)
                call kernel%evaluate(cell, j, first, nq, jet, wall(:, j), 2, w, g, &
                    hall(:, :, j), status)
                if (status /= 0) then
                    !$omp atomic write
                    ierr = status
                end if
            end do
            lp = 0.0_dp
            do q = 1, nq
                m4 = 0.0_dp
                do j = 1, tor%plane_count
                    do d2 = 1, nc3/nf
                        do f2 = 1, nf
                            c2 = f2 + nf*(d2 - 1)
                            do d = 1, nc3/nf
                                do f = 1, nf
                                    c = f + nf*(d - 1)
                                    tw = hall(c + nc3*(c2 - 1), q, j)*wall(q, j)
                                    if (tw == 0.0_dp) cycle
                                    do m2 = 1, nm
                                        wb = tor%table(m2, j, dt_of(d2))
                                        if (wb == 0.0_dp) cycle
                                        do m = 1, nm
                                            wa = tor%table(m, j, dt_of(d))
                                            if (wa == 0.0_dp) cycle
                                            m4(f + nf*(m - 1), c2_of(d), &
                                                f2 + nf*(m2 - 1), c2_of(d2)) = &
                                                m4(f + nf*(m - 1), c2_of(d), &
                                                f2 + nf*(m2 - 1), c2_of(d2)) + wa*tw*wb
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
                do i = 1, nb
                    do pp = pptr(i), pptr(i + 1) - 1
                        lp(:, :, pp) = lp(:, :, pp) + matmul(phi(:, :, q), &
                            matmul(m4(i, :, pcol(pp), :), transpose(phi(:, :, q))))
                    end do
                end do
            end do
            !$omp critical (spline_energy3_scatter)
            do a = 1, nloc
                do pa = space%broken_ptr(broken(a)), space%broken_ptr(broken(a) + 1) - 1
                    ia = space%broken_dof(pa)
                    do r = rowptr(ia), rowptr(ia + 1) - 1
                        pos(col(r)) = r
                    end do
                    do a2 = 1, nloc
                        do pb = space%broken_ptr(broken(a2)), &
                                space%broken_ptr(broken(a2) + 1) - 1
                            ib = space%broken_dof(pb)
                            wab = space%broken_weight(pa)*space%broken_weight(pb)
                            blocks(:, pos(ib)) = blocks(:, pos(ib)) + wab*lp(a, a2, :)
                        end do
                    end do
                end do
            end do
            !$omp end critical (spline_energy3_scatter)
        end do
        !$omp end do
        !$omp end parallel
        if (ierr /= 0) then
            status = ierr
            return
        end if
        matrix%n = nb*space%dof_count
        allocate(matrix%rowptr(matrix%n + 1), matrix%colidx(nnz*np), &
            matrix%val(nnz*np))
        k = 0
        matrix%rowptr(1) = 1
        do ia = 1, space%dof_count
            do i = 1, nb
                do r = rowptr(ia), rowptr(ia + 1) - 1
                    do pp = pptr(i), pptr(i + 1) - 1
                        k = k + 1
                        matrix%colidx(k) = pcol(pp) + nb*(col(r) - 1)
                        matrix%val(k) = blocks(pp, r)
                    end do
                end do
                matrix%rowptr(i + nb*(ia - 1) + 1) = k + 1
            end do
        end do
        status = 0
    end subroutine spline_energy3_hessian_csr

    subroutine section_pattern(space, rowptr, col)
        !! Sorted section dof graph: I, J coupled when both are nonzero on a
        !! common cell.
        type(spline_space_t), intent(in) :: space
        integer, allocatable, intent(out) :: rowptr(:), col(:)

        integer, allocatable :: cptr(:), cdof(:), dptr(:), dcell(:), mark(:), next(:)
        integer :: nloc, cell, a, p, i, n, r, c, s, pass, x, y
        integer :: broken((space%degree + 1)**2)

        nloc = (space%degree + 1)**2
        allocate(cptr(space%cell_count + 1), mark(space%dof_count))
        mark = 0
        do pass = 1, 2
            n = 0
            do cell = 1, space%cell_count
                if (pass == 1) cptr(cell) = n + 1
                call spline_cell_broken(space, cell, broken)
                do a = 1, nloc
                    do p = space%broken_ptr(broken(a)), space%broken_ptr(broken(a) + 1) - 1
                        i = space%broken_dof(p)
                        if (mark(i) == cell + (pass - 1)*space%cell_count) cycle
                        mark(i) = cell + (pass - 1)*space%cell_count
                        n = n + 1
                        if (pass == 2) cdof(n) = i
                    end do
                end do
            end do
            if (pass == 1) then
                cptr(space%cell_count + 1) = n + 1
                allocate(cdof(n))
            end if
        end do
        allocate(dptr(space%dof_count + 1), next(space%dof_count + 1))
        dptr = 0
        do r = 1, size(cdof)
            dptr(cdof(r) + 1) = dptr(cdof(r) + 1) + 1
        end do
        dptr(1) = 1
        do i = 2, space%dof_count + 1
            dptr(i) = dptr(i) + dptr(i - 1)
        end do
        allocate(dcell(size(cdof)))
        next = dptr
        do cell = 1, space%cell_count
            do r = cptr(cell), cptr(cell + 1) - 1
                dcell(next(cdof(r))) = cell
                next(cdof(r)) = next(cdof(r)) + 1
            end do
        end do
        mark = 0
        allocate(rowptr(space%dof_count + 1))
        do pass = 1, 2
            n = 0
            do i = 1, space%dof_count
                if (pass == 1) rowptr(i) = n + 1
                s = n
                do r = dptr(i), dptr(i + 1) - 1
                    c = dcell(r)
                    do p = cptr(c), cptr(c + 1) - 1
                        if (mark(cdof(p)) == i + (pass - 1)*space%dof_count) cycle
                        mark(cdof(p)) = i + (pass - 1)*space%dof_count
                        n = n + 1
                        if (pass == 2) col(n) = cdof(p)
                    end do
                end do
                if (pass == 2) then
                    do x = s + 2, n
                        y = col(x)
                        p = x - 1
                        do while (p > s)
                            if (col(p) <= y) exit
                            col(p + 1) = col(p)
                            p = p - 1
                        end do
                        col(p + 1) = y
                    end do
                end if
            end do
            if (pass == 1) then
                rowptr(space%dof_count + 1) = n + 1
                allocate(col(n))
            end if
        end do
    end subroutine section_pattern

    subroutine spline_evaluate3(space, quad, tor, coef, order, jet, status)
        !! jet(nfield*count(order), point, plane) at all points and planes.
        type(spline_space_t), intent(in) :: space
        type(spline_quadrature_t), intent(in) :: quad
        type(spline_toroidal_t), intent(in) :: tor
        real(dp), intent(in) :: coef(:, :, :)
        integer, intent(in) :: order
        real(dp), intent(out) :: jet(:, :, :)
        integer, intent(out) :: status

        real(dp), allocatable :: bc(:, :), phi(:, :, :), u(:, :, :)
        integer, allocatable :: broken(:)
        integer :: nf, nm, nb, nc3, nd2, nq, nloc, cell, j, first

        status = 1
        nf = size(coef, 1)
        nm = tor%mode_count
        nb = nf*nm
        if (order < 0 .or. order > 2 .or. size(coef, 2) /= nm) return
        nc3 = nf*SPLINE_JET3_COUNT(order)
        nd2 = section_count(order)
        nq = quad%points_per_cell
        nloc = (space%degree + 1)**2
        allocate(bc(nb, space%broken_count), phi(nloc, nd2, nq), broken(nloc), &
            u(nb, nd2, nq))
        call broken_of(space, reshape(coef, [nb, space%dof_count]), nb, bc)
        do cell = 1, space%cell_count
            first = quad%point_offset(cell)
            call cell_basis(space, quad, cell, nd2, phi, broken)
            call section_jets(bc(:, broken), phi, u)
            do j = 1, tor%plane_count
                call plane_jets(tor, j, nf, nc3, u, jet(:, first:first + nq - 1, j))
            end do
        end do
        status = 0
    end subroutine spline_evaluate3

    subroutine spline_evaluate3_at(space, tor, coef, block, xi, t, order, jet, status)
        !! Jet at (block, xi) and t.
        type(spline_space_t), intent(in) :: space
        type(spline_toroidal_t), intent(in) :: tor
        real(dp), intent(in) :: coef(:, :, :), xi(2), t
        integer, intent(in) :: block, order
        real(dp), intent(out) :: jet(:)
        integer, intent(out) :: status

        real(dp), allocatable :: section(:), values(:, :), u(:, :, :)
        integer :: nf, nm, nb, nd2, c, d, f, m

        status = 1
        nf = size(coef, 1)
        nm = tor%mode_count
        nb = nf*nm
        if (order < 0 .or. order > 2) return
        nd2 = section_count(order)
        allocate(section(nb*nd2), values(nm, 0:2), u(nb, nd2, 1))
        call spline_evaluate_at(space, reshape(coef, [nb, space%dof_count]), block, xi, &
            order, section, status)
        if (status /= 0) return
        do c = 1, nd2
            u(:, c, 1) = section(nb*(c - 1) + 1:nb*c)
        end do
        call spline_toroidal_values(tor, t, values)
        jet = 0.0_dp
        do d = 1, SPLINE_JET3_COUNT(order)
            do m = 1, nm
                do f = 1, nf
                    jet(f + nf*(d - 1)) = jet(f + nf*(d - 1)) + &
                        values(m, dt_of(d))*u(f + nf*(m - 1), c2_of(d), 1)
                end do
            end do
        end do
        status = 0
    end subroutine spline_evaluate3_at

end module fortfem_spline_toroidal
