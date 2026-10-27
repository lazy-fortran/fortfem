program test_element_to_csr
    !! Element blocks to CSR on curved isoparametric P1/P2/P3 triangle meshes:
    !! pattern symmetry, CSR product vs element-by-element scatter apply, exchange
    !! form, colouring, AD products, thread invariance, the cached slots with the
    !! symmetric column path, and long (heapsorted) rows. Scaling is measured by
    !! benchmark/fortfem_csr_spmv_scaling.f90.
    use, intrinsic :: iso_fortran_env, only: int64
    use check, only: check_condition, check_summary
    use fortfem_element_to_csr, only: csr_assemble_blocks, csr_build_pattern, &
        csr_element_colouring, csr_find_entry
    use fortfem_feec, only: evaluate_triangle_lagrange_basis, &
        initialize_triangle_lagrange_basis, triangle_lagrange_basis_t, &
        triangle_lagrange_nodes
    use fortfem_kinds, only: dp
    use fortfem_sparse_matrix, only: csr_matrix_t, spmv, spmv_jvp, spmv_vjp
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    integer, parameter :: thread_counts(3) = [1, 2, 8]
    integer :: degree, hardware_threads

    hardware_threads = max_threads()
    do degree = 1, 3
        call run_degree(degree, 24)
    end do
    call check_long_rows()
    call check_summary("element_to_csr")

contains

    subroutine run_degree(degree, nx)
        integer, intent(in) :: degree, nx

        type(csr_matrix_t) :: exchange, full, tangent
        integer, allocatable :: colour(:), conn(:, :)
        real(dp), allocatable :: blocks_full(:, :, :), blocks_stiff(:, :, :)
        real(dp), allocatable :: reference(:), x(:), y(:), xdot(:), vdot(:)
        real(dp), allocatable :: vbar(:), xbar(:), ybar(:), ydot(:)
        real(dp) :: defect, lhs, rhs, scale
        integer :: e, n_colours, n_dofs, pass
        character(len=24) :: tag
        logical :: ok

        write (tag, "(a,i0)") "P", degree
        call build_mesh(degree, nx, conn, blocks_stiff, blocks_full)
        n_dofs = maxval(conn)

        call csr_build_pattern(n_dofs, conn, .false., full)
        call csr_build_pattern(n_dofs, conn, .true., exchange)
        call csr_assemble_blocks(full, conn, blocks_full)
        call csr_assemble_blocks(exchange, conn, blocks_stiff, defect)

        call check_condition(pattern_ok(full) .and. pattern_ok(exchange), &
            trim(tag)//": pattern symmetric, sorted, tpos consistent")
        call check_condition(covers_elements(full, conn) .and. &
            covers_elements(exchange, conn), &
            trim(tag)//": pattern contains every dof pair sharing an element")
        call check_condition(size(full%val) == size(exchange%val) + n_dofs, &
            trim(tag)//": exchange form drops exactly the diagonal")

        allocate(x(n_dofs), y(n_dofs), reference(n_dofs))
        call random_vector(x, 17)
        call element_apply(conn, blocks_full, x, reference)
        call spmv(full, x, y)
        call check_condition(maxval(abs(y - reference)) <= &
            1.0e-14_dp*maxval(abs(reference)), &
            trim(tag)//": full CSR apply equals element apply to 1e-14")
        call element_apply(conn, blocks_stiff, x, reference)
        call spmv(exchange, x, y)
        call check_condition(maxval(abs(y - reference)) <= &
            1.0e-14_dp*maxval(abs(reference)), &
            trim(tag)//": exchange CSR apply equals element apply to 1e-14")
        x = 1.0_dp
        call spmv(exchange, x, y)
        call check_condition(all(y == 0.0_dp) .and. defect < 1.0e-14_dp, &
            trim(tag)//": exchange form has exactly zero row sums")

        allocate(xdot(n_dofs), ydot(n_dofs), ybar(n_dofs), xbar(n_dofs))
        ok = .true.
        do pass = 1, 2
            if (pass == 1) then
                call copy_matrix(full, tangent)
            else
                call copy_matrix(exchange, tangent)
            end if
            allocate(vdot(size(tangent%val)), vbar(size(tangent%val)))
            call random_vector(x, 3)
            call random_vector(xdot, 5)
            call random_vector(ybar, 7)
            call random_vector(vdot, 11)
            call spmv_jvp(tangent, x, vdot, xdot, ydot)
            call spmv_vjp(tangent, x, ybar, vbar, xbar)
            lhs = dot_product(ybar, ydot)
            rhs = dot_product(vbar, vdot) + dot_product(xbar, xdot)
            scale = abs(lhs) + abs(rhs) + tiny(1.0_dp)
            ok = ok .and. abs(lhs - rhs) <= 1.0e-12_dp*scale
            deallocate(vdot, vbar)
        end do
        call check_condition(ok, trim(tag)//": JVP and VJP satisfy the adjoint identity")

        call csr_element_colouring(full, conn, colour, n_colours)
        ok = n_colours <= 16
        do e = 1, size(conn, 2)
            ok = ok .and. colour(e) >= 1
        end do
        call check_condition(ok .and. colouring_valid(full, conn, colour), &
            trim(tag)//": element colouring has no shared dof within a colour")

        call check_thread_invariance(conn, blocks_full, blocks_stiff, tag)
        call check_slots(conn, blocks_full, blocks_stiff, tag)
    end subroutine run_degree

    subroutine check_slots(conn, blocks_full, blocks_stiff, tag)
        !! The cached slots give the values of the search path bitwise, and with
        !! exactly symmetric blocks the symmetric (column) path does too.
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks_full(:, :, :), blocks_stiff(:, :, :)
        character(len=*), intent(in) :: tag

        type(csr_matrix_t) :: reference, slotted, ex_reference, ex_slotted
        real(dp), allocatable :: sym_full(:, :, :), sym_stiff(:, :, :)
        integer :: e, n_dofs
        logical :: same

        n_dofs = maxval(conn)
        sym_full = blocks_full
        sym_stiff = blocks_stiff
        do e = 1, size(conn, 2)
            sym_full(:, :, e) = 0.5_dp*(blocks_full(:, :, e) &
                + transpose(blocks_full(:, :, e)))
            sym_stiff(:, :, e) = 0.5_dp*(blocks_stiff(:, :, e) &
                + transpose(blocks_stiff(:, :, e)))
        end do
        call csr_build_pattern(n_dofs, conn, .false., reference)
        call csr_build_pattern(n_dofs, conn, .false., slotted, with_slots=.true.)
        call csr_build_pattern(n_dofs, conn, .true., ex_reference)
        call csr_build_pattern(n_dofs, conn, .true., ex_slotted, with_slots=.true.)
        call csr_assemble_blocks(reference, conn, blocks_full)
        call csr_assemble_blocks(slotted, conn, blocks_full)
        same = all(slotted%colidx == reference%colidx) .and. &
            all(slotted%val == reference%val)
        call csr_assemble_blocks(reference, conn, sym_full)
        call csr_assemble_blocks(slotted, conn, sym_full, symmetric=.true.)
        same = same .and. all(slotted%val == reference%val)
        call csr_assemble_blocks(ex_reference, conn, sym_stiff)
        call csr_assemble_blocks(ex_slotted, conn, sym_stiff, symmetric=.true.)
        same = same .and. all(ex_slotted%val == ex_reference%val)
        call check_condition(same, trim(tag)// &
            ": slots and the symmetric column path bitwise equal to the search path")
    end subroutine check_slots

    subroutine check_long_rows()
        !! Rows longer than the insertion-sort cutoff (heapsort path) are sorted,
        !! unique and hold exactly the dofs that share an element.
        integer, parameter :: n_dofs = 300, nel = 40, nloc = 60
        type(csr_matrix_t) :: matrix
        integer :: conn(nloc, nel), a, b, e, i, k, count_row
        logical :: member(n_dofs), ok
        real(dp) :: r(nloc*nel)

        call random_vector(r, 29)
        do e = 1, nel
            member = .false.
            a = 0
            k = 0
            do while (a < nloc)
                k = k + 1
                i = 1 + mod(int(abs(r(1 + mod(k + nloc*(e - 1), nloc*nel)))*1.0e6_dp) &
                    + 37*k + 11*e, n_dofs)
                if (member(i)) cycle
                member(i) = .true.
                a = a + 1
                conn(a, e) = i
            end do
        end do
        call csr_build_pattern(n_dofs, conn, .false., matrix, with_slots=.true.)
        ok = pattern_ok(matrix) .and. covers_elements(matrix, conn)
        do i = 1, n_dofs
            member = .false.
            do e = 1, nel
                if (.not. any(conn(:, e) == i)) cycle
                do b = 1, nloc
                    member(conn(b, e)) = .true.
                end do
            end do
            member(i) = .true.
            count_row = matrix%rowptr(i + 1) - matrix%rowptr(i)
            if (count_row /= count(member)) ok = .false.
        end do
        ok = ok .and. maxval(matrix%rowptr(2:) - matrix%rowptr(:n_dofs)) > 32
        call check_condition(ok, "long rows: heapsorted pattern complete and sorted")
    end subroutine check_long_rows

    subroutine check_thread_invariance(conn, blocks_full, blocks_stiff, tag)
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks_full(:, :, :), blocks_stiff(:, :, :)
        character(len=*), intent(in) :: tag

        type(csr_matrix_t) :: exchange, full
        real(dp), allocatable :: x(:), y(:), y_ref(:), z(:), z_ref(:)
        real(dp), allocatable :: val_ref(:), val_ex_ref(:)
        integer, allocatable :: col_ref(:)
        integer :: it, n_dofs
        logical :: same

        n_dofs = maxval(conn)
        allocate(x(n_dofs), y(n_dofs), y_ref(n_dofs), z(n_dofs), z_ref(n_dofs))
        call random_vector(x, 23)
        same = .true.
        do it = 1, size(thread_counts)
            call set_threads(thread_counts(it))
            call csr_build_pattern(n_dofs, conn, .false., full)
            call csr_build_pattern(n_dofs, conn, .true., exchange)
            call csr_assemble_blocks(full, conn, blocks_full)
            call csr_assemble_blocks(exchange, conn, blocks_stiff)
            call spmv(full, x, y)
            call spmv(exchange, x, z)
            if (it == 1) then
                y_ref = y
                z_ref = z
                val_ref = full%val
                val_ex_ref = exchange%val
                col_ref = full%colidx
            else
                same = same .and. all(full%colidx == col_ref) .and. &
                    all(full%val == val_ref) .and. all(exchange%val == val_ex_ref) &
                    .and. all(y == y_ref) .and. all(z == z_ref)
            end if
        end do
        call set_threads(hardware_threads)
        call check_condition(same, trim(tag)// &
            ": pattern, assembly and SpMV bitwise equal for 1, 2, 8 threads")
    end subroutine check_thread_invariance

    ! ---------------------------------------------------------------- mesh

    subroutine build_mesh(degree, nx, conn, blocks_stiff, blocks_full)
        !! Unit square split into 2*nx^2 triangles, nodes on a lattice of spacing
        !! 1/(degree*nx) pushed through a smooth non-affine map. Each element uses
        !! the isoparametric degree-p map, so its edges curve for p >= 2.
        integer, intent(in) :: degree, nx
        integer, allocatable, intent(out) :: conn(:, :)
        real(dp), allocatable, intent(out) :: blocks_stiff(:, :, :)
        real(dp), allocatable, intent(out) :: blocks_full(:, :, :)

        type(triangle_lagrange_basis_t) :: basis
        real(dp), allocatable :: coords(:, :, :), gradients(:, :, :, :)
        real(dp), allocatable :: nodes(:, :), shape_values(:, :), weights(:)
        real(dp), allocatable :: gauss_x(:), gauss_w(:), values(:), grads(:, :)
        real(dp) :: det, h, inverse(2, 2), jac(2, 2), px, py, xi, eta
        integer :: a, b, e, half, i, ia, ib, j, iq, jq, m, nloc, nq, q, side, status
        integer :: lattice

        nloc = (degree + 1)*(degree + 2)/2
        call initialize_triangle_lagrange_basis(degree, basis, status)
        call triangle_lagrange_nodes(basis, nodes)
        lattice = degree*nx + 1
        h = 1.0_dp/real(degree*nx, dp)

        ! Collapsed Gauss-Legendre rule on the reference triangle.
        call gauss_legendre(6, gauss_x, gauss_w)
        nq = 36
        allocate(weights(nq), shape_values(nloc, nq), gradients(2, nloc, nq, 1))
        allocate(values(nloc), grads(2, nloc))
        q = 0
        do iq = 1, 6
            do jq = 1, 6
                q = q + 1
                xi = gauss_x(iq)
                eta = gauss_x(jq)*(1.0_dp - xi)
                weights(q) = gauss_w(iq)*gauss_w(jq)*(1.0_dp - xi)
                call evaluate_triangle_lagrange_basis( &
                    basis, xi, eta, values, grads, status)
                shape_values(:, q) = values
                gradients(:, :, q, 1) = grads
            end do
        end do

        allocate(conn(nloc, 2*nx*nx), coords(2, nloc, 2*nx*nx))
        allocate(blocks_stiff(nloc, nloc, 2*nx*nx), blocks_full(nloc, nloc, 2*nx*nx))
        e = 0
        do j = 0, nx - 1
            do i = 0, nx - 1
                do half = 0, 1
                    e = e + 1
                    side = 1 - 2*half
                    do m = 1, nloc
                        ia = i*degree + half*degree + side*nint(nodes(1, m)*degree)
                        ib = j*degree + half*degree + side*nint(nodes(2, m)*degree)
                        conn(m, e) = ia + lattice*ib + 1
                        call curved_map(real(ia, dp)*h, real(ib, dp)*h, px, py)
                        coords(:, m, e) = [px, py]
                    end do
                end do
            end do
        end do
        do e = 1, size(conn, 2)
            blocks_stiff(:, :, e) = 0.0_dp
            blocks_full(:, :, e) = 0.0_dp
            do q = 1, nq
                jac = 0.0_dp
                do m = 1, nloc
                    jac(1, 1) = jac(1, 1) + coords(1, m, e)*gradients(1, m, q, 1)
                    jac(1, 2) = jac(1, 2) + coords(1, m, e)*gradients(2, m, q, 1)
                    jac(2, 1) = jac(2, 1) + coords(2, m, e)*gradients(1, m, q, 1)
                    jac(2, 2) = jac(2, 2) + coords(2, m, e)*gradients(2, m, q, 1)
                end do
                det = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
                if (det <= 0.0_dp) error stop "test mesh: inverted element"
                inverse(1, :) = [jac(2, 2), -jac(1, 2)]/det
                inverse(2, :) = [-jac(2, 1), jac(1, 1)]/det
                do a = 1, nloc
                    grads(:, a) = [ &
                        inverse(1, 1)*gradients(1, a, q, 1) + &
                        inverse(2, 1)*gradients(2, a, q, 1), &
                        inverse(1, 2)*gradients(1, a, q, 1) + &
                        inverse(2, 2)*gradients(2, a, q, 1)]
                end do
                do b = 1, nloc
                    do a = 1, nloc
                        blocks_stiff(a, b, e) = blocks_stiff(a, b, e) + &
                            weights(q)*det*dot_product(grads(:, a), grads(:, b))
                        blocks_full(a, b, e) = blocks_full(a, b, e) + &
                            weights(q)*det*(dot_product(grads(:, a), grads(:, b)) &
                            + shape_values(a, q)*shape_values(b, q))
                    end do
                end do
            end do
        end do
    end subroutine build_mesh

    pure subroutine curved_map(u, v, px, py)
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: px, py

        real(dp), parameter :: two_pi = 6.283185307179586_dp

        px = u + 0.08_dp*sin(two_pi*v) + 0.05_dp*u*v
        py = v + 0.08_dp*sin(two_pi*u) - 0.04_dp*u*u
    end subroutine curved_map

    subroutine gauss_legendre(n, x, w)
        !! Gauss-Legendre points and weights on [0, 1].
        integer, intent(in) :: n
        real(dp), allocatable, intent(out) :: x(:), w(:)

        real(dp), parameter :: pi = 3.141592653589793_dp
        real(dp) :: dp0, p0, p1, p2, t
        integer :: i, iter, k

        allocate(x(n), w(n))
        do i = 1, n
            t = cos(pi*(real(i, dp) - 0.25_dp)/(real(n, dp) + 0.5_dp))
            do iter = 1, 100
                p0 = 1.0_dp
                p1 = t
                do k = 2, n
                    p2 = ((2*k - 1)*t*p1 - (k - 1)*p0)/k
                    p0 = p1
                    p1 = p2
                end do
                dp0 = n*(t*p1 - p0)/(t*t - 1.0_dp)
                t = t - p1/dp0
                if (abs(p1/dp0) < 1.0e-15_dp) exit
            end do
            x(i) = 0.5_dp*(1.0_dp + t)
            w(i) = 1.0_dp/((1.0_dp - t*t)*dp0*dp0)
        end do
    end subroutine gauss_legendre

    ! ------------------------------------------------------------- helpers

    subroutine element_apply(conn, blocks, x, y)
        !! Reference apply: serial element-by-element scatter.
        integer, intent(in) :: conn(:, :)
        real(dp), intent(in) :: blocks(:, :, :), x(:)
        real(dp), intent(out) :: y(:)

        integer :: a, b, e

        y = 0.0_dp
        do e = 1, size(conn, 2)
            do b = 1, size(conn, 1)
                do a = 1, size(conn, 1)
                    y(conn(a, e)) = y(conn(a, e)) + blocks(a, b, e)*x(conn(b, e))
                end do
            end do
        end do
    end subroutine element_apply

    logical function pattern_ok(matrix)
        type(csr_matrix_t), intent(in) :: matrix

        integer :: i, k, p

        pattern_ok = .true.
        do i = 1, matrix%n
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                if (k > matrix%rowptr(i)) then
                    if (matrix%colidx(k) <= matrix%colidx(k - 1)) pattern_ok = .false.
                end if
                if (matrix%exchange .and. matrix%colidx(k) == i) pattern_ok = .false.
                p = csr_find_entry(matrix, matrix%colidx(k), i)
                if (p < 1) then
                    pattern_ok = .false.
                else
                    if (matrix%tpos(k) /= p) pattern_ok = .false.
                end if
            end do
            if (.not. matrix%exchange) then
                if (csr_find_entry(matrix, i, i) < 1) pattern_ok = .false.
            end if
        end do
    end function pattern_ok

    logical function covers_elements(matrix, conn)
        type(csr_matrix_t), intent(in) :: matrix
        integer, intent(in) :: conn(:, :)

        integer :: a, b, e

        covers_elements = .true.
        do e = 1, size(conn, 2)
            do b = 1, size(conn, 1)
                do a = 1, size(conn, 1)
                    if (matrix%exchange .and. conn(a, e) == conn(b, e)) cycle
                    if (csr_find_entry(matrix, conn(a, e), conn(b, e)) < 1) then
                        covers_elements = .false.
                    end if
                end do
            end do
        end do
    end function covers_elements

    logical function colouring_valid(matrix, conn, colour)
        type(csr_matrix_t), intent(in) :: matrix
        integer, intent(in) :: conn(:, :), colour(:)

        integer :: dof, k, l

        colouring_valid = .true.
        do dof = 1, matrix%n
            do k = matrix%elem_ptr(dof), matrix%elem_ptr(dof + 1) - 1
                do l = k + 1, matrix%elem_ptr(dof + 1) - 1
                    if (colour(matrix%elem_idx(k)) == colour(matrix%elem_idx(l))) then
                        colouring_valid = .false.
                    end if
                end do
            end do
        end do
        if (size(conn) < 0) colouring_valid = .false.
    end function colouring_valid

    subroutine copy_matrix(source, target)
        type(csr_matrix_t), intent(in) :: source
        type(csr_matrix_t), intent(out) :: target

        integer :: k

        target = source
        do k = 1, size(source%val)
            target%val(k) = source%val(k)*(1.0_dp + 0.25_dp*cos(real(k, dp)))
        end do
    end subroutine copy_matrix

    subroutine random_vector(v, seed)
        real(dp), intent(out) :: v(:)
        integer, intent(in) :: seed

        integer(int64) :: state
        integer :: i

        state = 12345_int64 + 7919_int64*seed
        do i = 1, size(v)
            state = mod(state*6364136223846793005_int64 + 1442695040888963407_int64, &
                2_int64**62)
            v(i) = real(mod(state/1024_int64, 2000001_int64), dp)/1.0e6_dp - 1.0_dp
        end do
    end subroutine random_vector

    subroutine set_threads(n)
        integer, intent(in) :: n

        if (n < 1) error stop "set_threads: invalid thread count"
        !$ call omp_set_num_threads(n)
    end subroutine set_threads

    integer function max_threads()
        max_threads = 1
        !$ max_threads = omp_get_max_threads()
    end function max_threads

end program test_element_to_csr
