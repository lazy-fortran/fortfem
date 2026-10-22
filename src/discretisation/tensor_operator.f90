module fortfem_tensor_operator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_generated_tensor_determinant, only: generated_tensor_determinant
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    use fortfem_gll_line, only: build_gll_line
    use fortfem_line_mesh, only: line_geometry_count, line_geometry_tables, unit_line, &
                                 line_matrix_t
    use fortfem_tensor_contractions, only: tensor_contract
    use fortfem_generated_tensor_metric, only: generated_tensor_metric
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
                            csr_find_entry, build_node_to_element, csr_element_colouring
    implicit none
    private
    public :: tensor_operator_t, build_tensor_operator, apply_tensor_operator
    public :: assemble_tensor_operator, tensor_diagonal, tensor_load, tensor_l2_error
    public :: diffusion_coefficient, scalar_function, metric_at, determinant_at
    public :: tensor_values, tensor_boundary_projection
    type :: line_transpose_t
        real(dp), allocatable :: bt(:, :, :), dt(:, :, :) ! (nloc, nq, ncell)
    end type
    ! One operator for every tensor family: per-direction line tables (GLL,
    ! B-spline, Fourier or mixed), cell(:, e) selects each direction's 1D cell.
    type :: tensor_operator_t
        type(tensor_mesh_t) :: mesh
        type(tensor_line_tables_t) :: line(3)
        type(line_transpose_t) :: trans(3)
        integer, allocatable :: cell(:, :)
        integer :: ni(3) = 1, no(3) = 1, ncolours = 0
        logical :: precompute = .true.
        real(dp), allocatable :: jac(:, :, :, :), coef(:, :, :, :), reaction(:, :)
        real(dp), allocatable :: weight(:, :), points(:, :, :), volume(:, :)
        real(dp), allocatable :: metric(:, :, :, :), mass(:, :)
        integer, allocatable :: colour_ptr(:), colour_cells(:)
    end type
    abstract interface
        subroutine diffusion_coefficient(x, k, c)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: k(:, :), c
        end subroutine
        function scalar_function(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function
    end interface
contains
    pure function determinant_at(j) result(det)
        real(dp), intent(in) :: j(3, 3)
        real(dp) :: det
        call generated_tensor_determinant(j(1, 1), j(2, 1), j(3, 1), &
                              j(1, 2), j(2, 2), j(3, 2), j(1, 3), j(2, 3), j(3, 3), det)
    end function

    subroutine metric_at(j, k, w, c, g, m, det)
        real(dp), intent(in) :: j(3, 3), k(3, 3), w, c
        real(dp), intent(out) :: g(3, 3), m, det
        det = determinant_at(j)
        g = 0; m = 0
        if (.not. ieee_is_finite(det)) return
        if (det <= 0) return
        call generated_tensor_metric(j(1, 1), j(2, 1), j(3, 1), &
                                 j(1, 2), j(2, 2), j(3, 2), j(1, 3), j(2, 3), j(3, 3), &
                                 k(1, 1), k(2, 1), k(3, 1), k(1, 2), k(2, 2), k(3, 2), &
                           k(1, 3), k(2, 3), k(3, 3), w, c, g(1, 1), g(2, 1), g(3, 1), &
                           g(1, 2), g(2, 2), g(3, 2), g(1, 3), g(2, 3), g(3, 3), m, det)
    end subroutine

    subroutine build_tensor_operator(mesh, coefficient, op, status, precompute, nq)
        ! GLL p-nodal meshes build GLL lines (nq Gauss points, default p+1) on the
        ! reference cell; line-table meshes carry their own lines and rules.
        type(tensor_mesh_t), intent(in) :: mesh
        procedure(diffusion_coefficient) :: coefficient
        type(tensor_operator_t), intent(out) :: op
        integer, intent(out) :: status
        logical, intent(in), optional :: precompute
        integer, intent(in), optional :: nq
        type(csr_matrix_t) :: adjacency
        type(line_matrix_t), allocatable :: gb(:, :), gd(:, :)
        real(dp), allocatable :: t1(:), t2(:), values(:)
        integer, allocatable :: colour(:), fill(:)
        real(dp) :: g(3, 3), m, det, kmat(3, 3), c
        integer :: ng(3), mq, e, comp, dir, i, j, l, point, q, d, ce(3), ws
        status = 1
        if (mesh%dim < 2 .or. mesh%dim > 3) return
        op%mesh = mesh
        if (present(precompute)) op%precompute = precompute
        q = mesh%geometry_degree
        if (allocated(mesh%cell)) then
            if (present(nq)) return
            op%line = mesh%line
            op%cell = mesh%cell
        else
            mq = mesh%degree + 1
            if (present(nq)) mq = nq
            do d = 1, mesh%dim
                call build_gll_line(mesh%degree, 1, op%line(d), status, nq=mq, &
                                    x0=-1.0_dp, x1=1.0_dp)
                if (status /= 0) return
            end do
            if (mesh%dim == 2) call unit_line(op%line(3))
            allocate (op%cell(3, mesh%nelem))
            op%cell = 1
        end if
        status = 1
        op%ni = op%line%nloc
        op%no = op%line%nq
        do d = 1, 3
            ng(d) = line_geometry_count(op%line(d), q)
            allocate (op%trans(d)%bt(op%ni(d), op%no(d), op%line(d)%nelem), &
                      op%trans(d)%dt(op%ni(d), op%no(d), op%line(d)%nelem))
            do i = 1, op%line(d)%nelem
                op%trans(d)%bt(:, :, i) = transpose(op%line(d)%b(:, :, i))
                op%trans(d)%dt(:, :, i) = transpose(op%line(d)%d(:, :, i))
            end do
        end do
        if (size(mesh%geometry, 2) /= product(ng)) return
        allocate (gb(maxval(op%line%nelem), 3), gd(maxval(op%line%nelem), 3))
        do d = 1, 3
            do i = 1, op%line(d)%nelem
             call line_geometry_tables(op%line(d), i, q, gb(i, d)%a, gd(i, d)%a, status)
                if (status /= 0) return
            end do
        end do
        mq = product(op%no)
        allocate (op%jac(3, 3, mq, mesh%nelem), &
                  op%coef(3, 3, mq, mesh%nelem), op%reaction(mq, mesh%nelem), &
                  op%points(mesh%dim, mq, mesh%nelem), op%weight(mq, mesh%nelem), &
                  op%volume(mq, mesh%nelem))
        if (op%precompute) allocate (op%metric(3, 3, mq, mesh%nelem), &
                                     op%mass(mq, mesh%nelem))
        ws = maxval([ng, op%no])**3
        allocate (t1(ws), t2(ws), values(mq))
        op%jac = 0
        if (mesh%dim == 2) op%jac(3, 3, :, :) = 1
        do e = 1, mesh%nelem
            ce = op%cell(:, e)
            do l = 1, op%no(3)
                do j = 1, op%no(2)
                    do i = 1, op%no(1)
                        point = i + op%no(1)*(j - 1 + op%no(2)*(l - 1))
                        op%weight(point, e) = op%line(1)%wq(i, ce(1))* &
                                         op%line(2)%wq(j, ce(2))*op%line(3)%wq(l, ce(3))
                    end do
                end do
            end do
            do comp = 1, mesh%dim
                call tensor_contract(ng, op%no, gb(ce(1), 1)%a, gb(ce(2), 2)%a, &
                              gb(ce(3), 3)%a, mesh%geometry(comp, :, e), values, t1, t2)
                op%points(comp, :, e) = values
                do dir = 1, mesh%dim
                    select case (dir)
                    case (1)
                       call tensor_contract(ng, op%no, gd(ce(1), 1)%a, gb(ce(2), 2)%a, &
                                    gb(ce(3), 3)%a, mesh%geometry(comp, :, e), values, &
                                             t1, t2)
                    case (2)
                       call tensor_contract(ng, op%no, gb(ce(1), 1)%a, gd(ce(2), 2)%a, &
                                    gb(ce(3), 3)%a, mesh%geometry(comp, :, e), values, &
                                             t1, t2)
                    case default
                       call tensor_contract(ng, op%no, gb(ce(1), 1)%a, gb(ce(2), 2)%a, &
                                    gd(ce(3), 3)%a, mesh%geometry(comp, :, e), values, &
                                             t1, t2)
                    end select
                    op%jac(comp, dir, :, e) = values
                end do
            end do
            do point = 1, mq
                kmat = 0
                kmat(3, 3) = 1
                call coefficient(op%points(:, point, e), &
                                 kmat(1:mesh%dim, 1:mesh%dim), c)
                if (.not. all(ieee_is_finite(kmat))) return
                if (.not. ieee_is_finite(c)) return
                if (c < 0) return
                if (maxval(abs(kmat - transpose(kmat))) > &
                    1.e-13_dp*maxval(abs(kmat))) return
                if (kmat(1, 1) <= 0) return
                if (kmat(1, 1)*kmat(2, 2) - kmat(1, 2)**2 <= 0) return
                if (determinant_at(kmat) <= 0) return
                op%coef(:, :, point, e) = kmat
                op%reaction(point, e) = c
                call metric_at(op%jac(:, :, point, e), kmat, op%weight(point, e), &
                               c, g, m, det)
                if (.not. ieee_is_finite(det) .or. det <= 0) then
                    status = 2
                    return
                end if
                op%volume(point, e) = op%weight(point, e)*det
                if (op%precompute) then
                    op%metric(:, :, point, e) = g
                    op%mass(point, e) = m
                end if
            end do
        end do
        ! Colour from node-to-cell adjacency only, without a dense CSR pattern.
        call build_node_to_element(mesh%ndof, mesh%dof, &
                                   adjacency%elem_ptr, adjacency%elem_idx)
        call csr_element_colouring(adjacency, mesh%dof, colour, op%ncolours)
        allocate (op%colour_ptr(op%ncolours + 1), &
                  op%colour_cells(mesh%nelem), fill(op%ncolours))
        op%colour_ptr = 0
        do e = 1, mesh%nelem
            op%colour_ptr(colour(e) + 1) = op%colour_ptr(colour(e) + 1) + 1
        end do
        op%colour_ptr(1) = 1
        do i = 1, op%ncolours
            op%colour_ptr(i + 1) = op%colour_ptr(i + 1) + op%colour_ptr(i)
        end do
        fill = op%colour_ptr(1:op%ncolours)
        do e = 1, mesh%nelem
            i = colour(e)
            op%colour_cells(fill(i)) = e
            fill(i) = fill(i) + 1
        end do
        status = 0
    end subroutine

    pure subroutine cell_values(op, c, kind, u, v, t1, t2)
        ! Sum-factorised interpolation of cell coefficients u to quadrature:
        ! values (kind 0) or the parameter derivative along direction `kind`.
        type(tensor_operator_t), intent(in) :: op
        integer, intent(in) :: c(3), kind
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        real(dp), intent(inout) :: t1(:), t2(:)
        select case (kind)
        case (1)
            call tensor_contract(op%ni, op%no, op%line(1)%d(:, :, c(1)), &
                                 op%line(2)%b(:, :, c(2)), op%line(3)%b(:, :, c(3)), &
                                 u, v, t1, t2)
        case (2)
            call tensor_contract(op%ni, op%no, op%line(1)%b(:, :, c(1)), &
                                 op%line(2)%d(:, :, c(2)), op%line(3)%b(:, :, c(3)), &
                                 u, v, t1, t2)
        case (3)
            call tensor_contract(op%ni, op%no, op%line(1)%b(:, :, c(1)), &
                                 op%line(2)%b(:, :, c(2)), op%line(3)%d(:, :, c(3)), &
                                 u, v, t1, t2)
        case default
            call tensor_contract(op%ni, op%no, op%line(1)%b(:, :, c(1)), &
                                 op%line(2)%b(:, :, c(2)), op%line(3)%b(:, :, c(3)), &
                                 u, v, t1, t2)
        end select
    end subroutine

    pure subroutine cell_tests(op, c, kind, f, v, t1, t2)
        ! Transposed contraction: quadrature data f against test functions.
        type(tensor_operator_t), intent(in) :: op
        integer, intent(in) :: c(3), kind
        real(dp), intent(in) :: f(:)
        real(dp), intent(out) :: v(:)
        real(dp), intent(inout) :: t1(:), t2(:)
        select case (kind)
        case (1)
            call tensor_contract(op%no, op%ni, op%trans(1)%dt(:, :, c(1)), &
                               op%trans(2)%bt(:, :, c(2)), op%trans(3)%bt(:, :, c(3)), &
                                 f, v, t1, t2)
        case (2)
            call tensor_contract(op%no, op%ni, op%trans(1)%bt(:, :, c(1)), &
                               op%trans(2)%dt(:, :, c(2)), op%trans(3)%bt(:, :, c(3)), &
                                 f, v, t1, t2)
        case (3)
            call tensor_contract(op%no, op%ni, op%trans(1)%bt(:, :, c(1)), &
                               op%trans(2)%bt(:, :, c(2)), op%trans(3)%dt(:, :, c(3)), &
                                 f, v, t1, t2)
        case default
            call tensor_contract(op%no, op%ni, op%trans(1)%bt(:, :, c(1)), &
                               op%trans(2)%bt(:, :, c(2)), op%trans(3)%bt(:, :, c(3)), &
                                 f, v, t1, t2)
        end select
    end subroutine

    subroutine apply_tensor_operator(op, x, y, eliminate)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        logical, intent(in), optional :: eliminate
        logical :: mask
        integer :: nl, mq, ws, colour, ce, e, a, dir, q, dim
        mask = .true.
        if (present(eliminate)) mask = eliminate
        dim = op%mesh%dim
        nl = product(op%ni); mq = product(op%no)
        ws = maxval([op%ni, op%no])**3
        y = 0
        !$omp parallel default(shared) private(colour,ce,e,a,dir,q)
        block
            real(dp), allocatable :: u(:), v(:), tmp(:), grad(:, :), flux(:, :), &
                                     t1(:), t2(:)
            real(dp) :: g(3, 3), mass, det
            allocate (u(nl), v(nl), tmp(max(nl, mq)), grad(mq, 4), flux(mq, 4), &
                      t1(ws), t2(ws))
            do colour = 1, op%ncolours
                !$omp do schedule(static)
                do ce = op%colour_ptr(colour), op%colour_ptr(colour + 1) - 1
                    e = op%colour_cells(ce)
                    do a = 1, nl
                        u(a) = x(op%mesh%dof(a, e))
                        if (mask) then
                            if (op%mesh%boundary(op%mesh%dof(a, e))) u(a) = 0
                        end if
                    end do
                    do dir = 1, dim + 1
                        call cell_values(op, op%cell(:, e), modulo(dir, dim + 1), u, &
                                         grad(:, dir), t1, t2)
                    end do
                    do q = 1, mq
                        if (op%precompute) then
                            g = op%metric(:, :, q, e); mass = op%mass(q, e)
                        else
                            call metric_at(op%jac(:, :, q, e), op%coef(:, :, q, e), &
                                       op%weight(q, e), op%reaction(q, e), g, mass, det)
                        end if
                        do dir = 1, dim
                            flux(q, dir) = dot_product(g(dir, 1:dim), grad(q, 1:dim))
                        end do
                        flux(q, dim + 1) = mass*grad(q, dim + 1)
                    end do
                    v = 0
                    do dir = 1, dim + 1
                        call cell_tests(op, op%cell(:, e), modulo(dir, dim + 1), &
                                        flux(:, dir), tmp, t1, t2)
                        v = v + tmp(1:nl)
                    end do
                    do a = 1, nl
                        y(op%mesh%dof(a, e)) = y(op%mesh%dof(a, e)) + v(a)
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        if (mask) then
            where (op%mesh%boundary) y = x
        end if
    end subroutine

    pure function decode(a, n) result(index)
        integer, intent(in) :: a, n(3)
        integer :: index(3)
        index = [modulo(a - 1, n(1)), modulo((a - 1)/n(1), n(2)), &
                 (a - 1)/(n(1)*n(2))] + 1
    end function

    subroutine assemble_tensor_operator(op, matrix)
        type(tensor_operator_t), intent(in) :: op
        type(csr_matrix_t), intent(out) :: matrix
        real(dp), allocatable :: basis(:, :), gradient(:, :, :), element(:, :)
        real(dp) :: g(3, 3), mass, det, val
        integer :: a, b, q, e, i, j, k, l, dir, ia(3), iq(3), nl, mq, entry, c(3)
        nl = product(op%ni); mq = product(op%no)
        allocate (basis(mq, nl), gradient(3, mq, nl), element(nl, nl))
        gradient = 0
        call csr_build_pattern(op%mesh%ndof, op%mesh%dof, .false., matrix)
        do e = 1, op%mesh%nelem
            c = op%cell(:, e)
            ! Deliberately unfactored quadrature oracle, independent of MF schedules.
            do a = 1, nl
                ia = decode(a, op%ni)
                do q = 1, mq
                    iq = decode(q, op%no)
                    val = 1
                    do dir = 1, 3
                        val = val*op%line(dir)%b(iq(dir), ia(dir), c(dir))
                    end do
                    basis(q, a) = val
                    do dir = 1, op%mesh%dim
                        val = 1
                        do l = 1, 3
                            if (l == dir) then
                                val = val*op%line(l)%d(iq(l), ia(l), c(l))
                            else
                                val = val*op%line(l)%b(iq(l), ia(l), c(l))
                            end if
                        end do
                        gradient(dir, q, a) = val
                    end do
                end do
            end do
            element = 0
            do q = 1, mq
                call metric_at(op%jac(:, :, q, e), op%coef(:, :, q, e), &
                               op%weight(q, e), op%reaction(q, e), g, mass, det)
                do b = 1, nl
                    do a = 1, b
                        element(a, b) = element(a, b) + mass*basis(q, a)*basis(q, b) + &
                            dot_product(gradient(:, q, a), matmul(g, gradient(:, q, b)))
                    end do
                end do
            end do
            ! FortSparse requires canonical, exactly symmetric CSC input.
            ! Compute each symmetric bilinear-form entry once, then mirror it.
            do b = 1, nl
                do a = 1, b - 1
                    element(b, a) = element(a, b)
                end do
            end do
            do b = 1, nl
                j = op%mesh%dof(b, e)
                do a = 1, nl
                    i = op%mesh%dof(a, e)
                    entry = csr_find_entry(matrix, i, j)
                    matrix%val(entry) = matrix%val(entry) + element(a, b)
                end do
            end do
        end do
        do i = 1, matrix%n
            do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                j = matrix%colidx(k)
                if (op%mesh%boundary(i) .or. op%mesh%boundary(j)) then
                    matrix%val(k) = 0
                    if (i == j) matrix%val(k) = 1
                end if
            end do
        end do
    end subroutine

    subroutine tensor_diagonal(op, diagonal)
        type(tensor_operator_t), intent(in) :: op
        real(dp), allocatable, intent(out) :: diagonal(:)
        integer :: nl, mq, ws, colour, ce, e, a, r, s, d, q, dim
        dim = op%mesh%dim
        nl = product(op%ni); mq = product(op%no); ws = maxval([op%ni, op%no])**3
        allocate (diagonal(op%mesh%ndof))
        diagonal = 0
        !$omp parallel default(shared) private(colour,ce,e,a,r,s,d,q)
        block
            real(dp), allocatable :: t1(:), t2(:), local(:), term(:), values(:)
            type(line_matrix_t) :: table(3)
            real(dp) :: g(3, 3), mass, det
            integer :: c(3)
            allocate (t1(ws), t2(ws), local(nl), term(nl), values(mq))
            do colour = 1, op%ncolours
                !$omp do schedule(static)
                do ce = op%colour_ptr(colour), op%colour_ptr(colour + 1) - 1
                    e = op%colour_cells(ce)
                    c = op%cell(:, e)
                    local = 0
                    do r = 1, dim + 1
                        do s = 1, dim + 1
                            if (max(r, s) == dim + 1 .and. r /= s) cycle
                            do d = 1, 3
                                associate (b => op%line(d)%b(:, :, c(d)), &
                                           dd => op%line(d)%d(:, :, c(d)))
                                    if (r == d .and. s == d .and. d <= dim) then
                                        table(d)%a = transpose(dd**2)
                                    else if ((r == d .or. s == d) .and. d <= dim) then
                                        table(d)%a = transpose(b*dd)
                                    else
                                        table(d)%a = transpose(b**2)
                                    end if
                                end associate
                            end do
                            do q = 1, mq
                               call metric_at(op%jac(:, :, q, e), op%coef(:, :, q, e), &
                                       op%weight(q, e), op%reaction(q, e), g, mass, det)
                                if (r == dim + 1) then
                                    values(q) = mass
                                else
                                    values(q) = g(r, s)
                                end if
                            end do
                            call tensor_contract(op%no, op%ni, table(1)%a, table(2)%a, &
                                                 table(3)%a, values, term, t1, t2)
                            local = local + term
                        end do
                    end do
                    do a = 1, nl
                        diagonal(op%mesh%dof(a, e)) = &
                            diagonal(op%mesh%dof(a, e)) + local(a)
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        where (op%mesh%boundary) diagonal = 1
    end subroutine

    subroutine tensor_load(op, forcing, rhs)
        type(tensor_operator_t), intent(in) :: op
        procedure(scalar_function) :: forcing
        real(dp), allocatable, intent(out) :: rhs(:)
        real(dp), allocatable :: values(:), local(:), t1(:), t2(:)
        integer :: e, q, a, ws
        ws = maxval([op%ni, op%no])**3
        allocate (rhs(op%mesh%ndof), values(product(op%no)), &
                  local(product(op%ni)), t1(ws), t2(ws))
        rhs = 0
        do e = 1, op%mesh%nelem
            do q = 1, product(op%no)
                values(q) = op%volume(q, e)*forcing(op%points(:, q, e))
            end do
            call cell_tests(op, op%cell(:, e), 0, values, local, t1, t2)
            do a = 1, product(op%ni)
                rhs(op%mesh%dof(a, e)) = rhs(op%mesh%dof(a, e)) + local(a)
            end do
        end do
        where (op%mesh%boundary) rhs = 0
    end subroutine

    subroutine tensor_values(op, x, values)
        ! Discrete field at the operator's quadrature points, (nq, nelem).
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), allocatable, intent(out) :: values(:, :)
        real(dp), allocatable :: local(:), t1(:), t2(:)
        integer :: e, ws
        ws = maxval([op%ni, op%no])**3
        allocate (values(product(op%no), op%mesh%nelem))
        !$omp parallel default(shared) private(e, local, t1, t2)
        allocate (local(product(op%ni)), t1(ws), t2(ws))
        !$omp do schedule(static)
        do e = 1, op%mesh%nelem
            local = x(op%mesh%dof(:, e))
            call cell_values(op, op%cell(:, e), 0, local, values(:, e), t1, t2)
        end do
        !$omp end do
        !$omp end parallel
    end subroutine

    function tensor_l2_error(op, x, exact) result(error)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        procedure(scalar_function) :: exact
        real(dp) :: error
        real(dp), allocatable :: values(:, :)
        integer :: e, q
        call tensor_values(op, x, values)
        error = 0
        do e = 1, op%mesh%nelem
            do q = 1, product(op%no)
           error = error + op%volume(q, e)*(values(q, e) - exact(op%points(:, q, e)))**2
            end do
        end do
        error = sqrt(error)
    end function

    subroutine tensor_boundary_projection(op, g, x, status, iterations)
        ! L2 projection of Dirichlet data g onto the trace of the discrete space,
        ! for line families whose end functions interpolate (GLL, open B-spline).
        ! Returns x with boundary coefficients set and interior entries zero.
        ! Boundary faces are cell faces whose dofs are all Dirichlet dofs.
        type(tensor_operator_t), intent(in) :: op
        procedure(scalar_function) :: g
        real(dp), allocatable, intent(out) :: x(:)
        integer, intent(out) :: status
        integer, intent(out), optional :: iterations
        type(line_matrix_t) :: gb(3), gd(3)
        type(line_matrix_t), allocatable :: matrices(:)
        real(dp), allocatable :: tangent(:, :, :), xs(:, :), values(:), t1(:), t2(:), &
                                 basis(:, :), local(:), rhs(:), r(:), z(:), p(:), &
                                 ap(:), diagonal(:), measure(:)
        integer, allocatable :: faces(:, :), dofs(:, :), compact(:)
        integer :: ng(3), nof(3), e, d, side, a, i, j, k, nface, nf, nfaces, nb, &
                   comp, t, other(2), ia(3), iq(3), ws, it, q, dim, ce(3)
        real(dp) :: rho, previous, alpha, normal(3), initial
        logical :: face
        status = 1
        dim = op%mesh%dim
        nb = count(op%mesh%boundary)
        allocate (x(op%mesh%ndof), compact(op%mesh%ndof))
        x = 0
        compact = 0
        j = 0
        do i = 1, op%mesh%ndof
            if (.not. op%mesh%boundary(i)) cycle
            j = j + 1
            compact(i) = j
        end do
        do d = 1, 3
            ng(d) = line_geometry_count(op%line(d), op%mesh%geometry_degree)
        end do
        ! Collect boundary faces (cell, direction, side).
        allocate (faces(3, 6*op%mesh%nelem))
        nfaces = 0
        do e = 1, op%mesh%nelem
            do d = 1, dim
                if (op%line(d)%periodic) cycle
                do side = 0, 1
                    face = .true.
                    do a = 1, product(op%ni)
                        ia = decode(a, op%ni)
                        if (ia(d) /= merge(1, op%ni(d), side == 0)) cycle
                        if (.not. op%mesh%boundary(op%mesh%dof(a, e))) face = .false.
                    end do
                    if (.not. face) cycle
                    nfaces = nfaces + 1
                    faces(:, nfaces) = [e, d, side]
                end do
            end do
        end do
        if (nb > 0 .and. nfaces == 0) return
        nf = maxval([(product(op%ni)/op%ni(d), d=1, dim)])
        allocate (matrices(nfaces), dofs(nf, nfaces), rhs(nb))
        dofs = 0
        rhs = 0
        ws = maxval([ng, op%no])**3
        allocate (t1(ws), t2(ws))
        do t = 1, nfaces
            e = faces(1, t); d = faces(2, t); side = faces(3, t)
            ce = op%cell(:, e)
            nof = op%no
            nof(d) = 1
            other = pack([1, 2, 3], [1, 2, 3] /= d)
            do i = 1, 3
                if (i == d) then
                 call line_geometry_tables(op%line(i), ce(i), op%mesh%geometry_degree, &
                                              gb(i)%a, gd(i)%a, status, &
                                              face=merge(-1.0_dp, 1.0_dp, side == 0))
                else
                 call line_geometry_tables(op%line(i), ce(i), op%mesh%geometry_degree, &
                                              gb(i)%a, gd(i)%a, status)
                end if
                if (status /= 0) return
            end do
            status = 1
            q = product(nof)
            if (allocated(xs)) deallocate (xs, tangent, values, measure)
            allocate (xs(dim, q), tangent(3, q, 2), values(q), measure(q))
            tangent = 0
            do comp = 1, dim
                call tensor_contract(ng, nof, gb(1)%a, gb(2)%a, gb(3)%a, &
                                     op%mesh%geometry(comp, :, e), values, t1, t2)
                xs(comp, :) = values
                do k = 1, dim - 1
                    select case (other(k))
                    case (1)
                        call tensor_contract(ng, nof, gd(1)%a, gb(2)%a, gb(3)%a, &
                                           op%mesh%geometry(comp, :, e), values, t1, t2)
                    case (2)
                        call tensor_contract(ng, nof, gb(1)%a, gd(2)%a, gb(3)%a, &
                                           op%mesh%geometry(comp, :, e), values, t1, t2)
                    case default
                        call tensor_contract(ng, nof, gb(1)%a, gb(2)%a, gd(3)%a, &
                                           op%mesh%geometry(comp, :, e), values, t1, t2)
                    end select
                    tangent(comp, :, k) = values
                end do
            end do
            ! Face dofs, their trace basis at the face points and the measure.
            nface = product(op%ni)/op%ni(d)
            if (allocated(basis)) deallocate (basis, local)
            allocate (basis(q, nface), local(nface))
            j = 0
            do a = 1, product(op%ni)
                ia = decode(a, op%ni)
                if (ia(d) /= merge(1, op%ni(d), side == 0)) cycle
                j = j + 1
                dofs(j, t) = op%mesh%dof(a, e)
                do i = 1, q
                    iq = decode(i, nof)
                    basis(i, j) = op%line(other(1))%b(iq(other(1)), ia(other(1)), &
                                                      ce(other(1)))* &
                                  op%line(other(2))%b(iq(other(2)), ia(other(2)), &
                                                      ce(other(2)))
                end do
            end do
            do i = 1, q
                iq = decode(i, nof)
                if (dim == 2) then
                    measure(i) = norm2(tangent(:, i, 1))
                else
                    normal = [tangent(2, i, 1)*tangent(3, i, 2) - &
                              tangent(3, i, 1)*tangent(2, i, 2), &
                              tangent(3, i, 1)*tangent(1, i, 2) - &
                              tangent(1, i, 1)*tangent(3, i, 2), &
                              tangent(1, i, 1)*tangent(2, i, 2) - &
                              tangent(2, i, 1)*tangent(1, i, 2)]
                    measure(i) = norm2(normal)
                end if
                measure(i) = measure(i)* &
                             op%line(other(1))%wq(iq(other(1)), ce(other(1)))* &
                             op%line(other(2))%wq(iq(other(2)), ce(other(2)))
            end do
            allocate (matrices(t)%a(nface, nface))
            matrices(t)%a = matmul(transpose(basis), spread(measure, 2, nface)*basis)
            do i = 1, q
                values(i) = measure(i)*g(xs(:, i))
            end do
            local = matmul(values, basis)
            do j = 1, nface
                rhs(compact(dofs(j, t))) = rhs(compact(dofs(j, t))) + local(j)
            end do
        end do
        ! Jacobi-PCG on the face-wise (matrix-free) boundary mass.
        allocate (r(nb), z(nb), p(nb), ap(nb), diagonal(nb), local(nb))
        diagonal = 0
        do t = 1, nfaces
            do j = 1, size(matrices(t)%a, 1)
                i = compact(dofs(j, t))
                diagonal(i) = diagonal(i) + matrices(t)%a(j, j)
            end do
        end do
        if (any(diagonal <= 0)) return
        local = 0
        r = rhs
        initial = max(norm2(rhs), tiny(1.0_dp))
        previous = 1
        do it = 1, 20*nb + 100
            if (norm2(r) <= 1.e-14_dp*initial) exit
            z = r/diagonal
            rho = dot_product(r, z)
            if (it == 1) then
                p = z
            else
                p = z + (rho/previous)*p
            end if
            call boundary_mass(p, ap)
            alpha = rho/dot_product(p, ap)
            local = local + alpha*p
            r = r - alpha*ap
            previous = rho
        end do
        if (present(iterations)) iterations = it
        if (norm2(r) > 1.e-12_dp*initial) return
        do i = 1, op%mesh%ndof
            if (compact(i) > 0) x(i) = local(compact(i))
        end do
        status = 0
    contains
        subroutine boundary_mass(u, v)
            real(dp), intent(in) :: u(:)
            real(dp), intent(out) :: v(:)
            integer :: tt, n
            v = 0
            do tt = 1, nfaces
                n = size(matrices(tt)%a, 1)
                v(compact(dofs(1:n, tt))) = v(compact(dofs(1:n, tt))) + &
                                       matmul(matrices(tt)%a, u(compact(dofs(1:n, tt))))
            end do
        end subroutine
    end subroutine
end module
