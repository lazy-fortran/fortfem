module fortfem_tensor_seam
    !! Periodic seams whose two faces do not match node for node.
    !!
    !! A field-aligned mesh closes a period with a poloidal twist that is in
    !! general not a multiple of the cell size. The slave face DOFs are then
    !! eliminated by an affine constraint u_s = sum_m C_sm u_m on the primary
    !! (all other) DOFs, and every operator acts as P^T A P with the global
    !! prolongation P = [I; C]. Two constraints are available:
    !!
    !! - SEAM_MORTAR: the standard mortar condition
    !!   integral psi_i (u_slave - u_master o T) = 0 for slave test functions
    !!   psi_i whose degree is lowered by one in elements along a Dirichlet rim.
    !!   Overlaps are integrated by a Gauss fan rule on the clipped curved
    !!   polygon, so no master or slave kink lies inside a quadrature cell.
    !! - SEAM_INTERPOLATION: slave nodes interpolate the master trace. This
    !!   is the pointwise-matching comparison; it is not optimal for p >= 2.
    !!
    !! A seam whose faces match node for node is closed by DOF identification
    !! instead, which is bit-for-bit the ordinary periodic gather map.
    !! The seam map T is a rigid motion y = Q x + t from slave to master
    !! coordinates. Optional translations replicate the master face for a
    !! second periodic direction, for example a periodic slab coordinate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_gll_line, only: gll_rule, gauss_rule, lagrange_table
    use fortfem_tensor_mesh, only: tensor_mesh_t, build_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, apply_tensor_operator, &
                                       assemble_tensor_operator, tensor_diagonal
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_generated_tensor_surface, only: generated_tensor_surface_jacobian
    implicit none
    private
    public :: tensor_seam_t, close_tensor_seam, identity_tensor_seam
    public :: coarsen_tensor_seam, seam_expand, seam_fold, apply_seam_operator
    public :: seam_diagonal, assemble_seam_operator
    public :: SEAM_MORTAR, SEAM_INTERPOLATION
    integer, parameter :: SEAM_MORTAR = 1, SEAM_INTERPOLATION = 2
    real(dp), parameter :: rigid_identity(3, 4) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1, &
                                                           0, 0, 0]*1.0_dp, [3, 4])
    type :: seam_cubature_t
        !! Master face, slave and master reference points and physical weights.
        integer, allocatable :: master(:)
        real(dp), allocatable :: slave_ab(:, :), master_ab(:, :), weight(:)
    end type
    type :: tensor_seam_t
        !! No default initialisation: identity_tensor_seam, store_definition
        !! and build_constraint set every scalar.
        integer :: method, nfull, nreduced, nslave
        logical :: conforming
        !! reduced(i) > 0: primary index of full DOF i; < 0: minus slave row.
        integer, allocatable :: reduced(:), primary(:), slave(:)
        integer, allocatable :: ptr(:), col(:), tptr(:), trow(:)
        real(dp), allocatable :: val(:), tval(:)
        logical, allocatable :: boundary(:)
        integer, allocatable :: master_faces(:, :), slave_faces(:, :)
        real(dp) :: motion(3, 4)
        real(dp), allocatable :: images(:, :)
        !! Largest relative defect of the overlap areas covering a slave cell.
        real(dp) :: area_defect
        !! Mortar overlap cubature per slave face; geometry only, so coarser
        !! p-levels on the same geometry reuse it.
        type(seam_cubature_t), allocatable :: cubature(:)
    end type
    type :: seam_face_t
        integer :: element = 0, dir = 0, side = 0
        real(dp), allocatable :: geometry(:, :)
        integer, allocatable :: dof(:)
        real(dp) :: center(3) = 0, radius = 0
        !! A part of a parent face: parent reference = origin + scale*own.
        real(dp) :: origin(2) = 0, scale = 1
    end type
    type :: face_split_t
        type(seam_face_t), allocatable :: part(:)
    end type
contains
    subroutine identity_tensor_seam(mesh, seam)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_seam_t), intent(out) :: seam
        integer :: i
        seam%method = SEAM_MORTAR
        seam%nfull = mesh%ndof
        seam%nreduced = mesh%ndof
        seam%nslave = 0
        seam%conforming = .true.
        seam%motion = rigid_identity
        seam%area_defect = 0
        allocate (seam%reduced(mesh%ndof), seam%primary(mesh%ndof), seam%slave(0))
        do i = 1, mesh%ndof
            seam%reduced(i) = i
            seam%primary(i) = i
        end do
        allocate (seam%ptr(1), seam%col(0), seam%val(0))
        seam%ptr = 1
        allocate (seam%tptr(mesh%ndof + 1), seam%trow(0), seam%tval(0))
        seam%tptr = 1
        seam%boundary = mesh%boundary
    end subroutine

    subroutine close_tensor_seam(mesh, master, slave, closed, seam, status, &
                                 motion, images, method)
        !! Close a periodic seam. master/slave list faces as (element,
        !! direction, side) with side 0 at reference -1 and 1 at +1.
        type(tensor_mesh_t), intent(in) :: mesh
        integer, intent(in) :: master(:, :), slave(:, :)
        type(tensor_mesh_t), intent(out) :: closed
        type(tensor_seam_t), intent(out) :: seam
        integer, intent(out) :: status
        real(dp), intent(in), optional :: motion(3, 4), images(:, :)
        integer, intent(in), optional :: method
        type(tensor_seam_t) :: setup
        integer, allocatable :: match(:)
        logical :: matched
        status = 1
        call store_definition(setup, master, slave, motion, images, method, status)
        if (status /= 0) return
        if (mesh%dim /= 3) then
            status = 1
            return
        end if
        call match_faces(mesh, setup, match, matched, status)
        if (status /= 0) return
        if (matched) then
            call merge_seam(mesh, match, closed, status)
            if (status /= 0) return
            call identity_tensor_seam(closed, seam)
        else
            closed = mesh
            call build_constraint(mesh, setup, seam, status)
        end if
    end subroutine

    subroutine coarsen_tensor_seam(fine, mesh, seam, status)
        !! The same seam on a mesh with the same cells and another degree.
        type(tensor_seam_t), intent(in) :: fine
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_seam_t), intent(out) :: seam
        integer, intent(out) :: status
        type(tensor_seam_t) :: setup
        status = 0
        if (fine%conforming) then
            call identity_tensor_seam(mesh, seam)
            return
        end if
        call store_definition(setup, fine%master_faces, fine%slave_faces, &
                              fine%motion, fine%images, fine%method, status)
        if (status /= 0) return
        if (allocated(fine%cubature)) setup%cubature = fine%cubature
        call build_constraint(mesh, setup, seam, status)
    end subroutine

    subroutine store_definition(seam, master, slave, motion, images, method, status)
        type(tensor_seam_t), intent(inout) :: seam
        integer, intent(in) :: master(:, :), slave(:, :)
        real(dp), intent(in), optional :: motion(3, 4), images(:, :)
        integer, intent(in), optional :: method
        integer, intent(out) :: status
        real(dp) :: q(3, 3)
        integer :: i
        status = 1
        if (size(master, 1) /= 3 .or. size(slave, 1) /= 3) return
        if (size(master, 2) < 1 .or. size(slave, 2) < 1) return
        if (any(master(2, :) < 1) .or. any(master(2, :) > 3)) return
        if (any(slave(2, :) < 1) .or. any(slave(2, :) > 3)) return
        if (any(master(3, :) < 0) .or. any(master(3, :) > 1)) return
        if (any(slave(3, :) < 0) .or. any(slave(3, :) > 1)) return
        seam%master_faces = master
        seam%slave_faces = slave
        seam%motion = rigid_identity
        if (present(motion)) seam%motion = motion
        q = matmul(transpose(seam%motion(:, 1:3)), seam%motion(:, 1:3))
        do i = 1, 3
            q(i, i) = q(i, i) - 1
        end do
        if (maxval(abs(q)) > 1.e-12_dp) return
        if (present(images)) then
            if (size(images, 1) /= 3) return
            seam%images = images
        else
            allocate (seam%images(3, 1))
            seam%images = 0
        end if
        seam%method = SEAM_MORTAR
        if (present(method)) seam%method = method
        if (seam%method /= SEAM_MORTAR .and. seam%method /= SEAM_INTERPOLATION) return
        status = 0
    end subroutine

    pure subroutine face_nodes(n, dir, side, idx)
        !! Volume-local indices of the n*n face nodes, first tangent fastest.
        integer, intent(in) :: n, dir, side
        integer, intent(out) :: idx(:)
        integer :: i, j, pos(3), t(2), k
        k = 0
        do i = 1, 3
            if (i == dir) cycle
            k = k + 1
            t(k) = i
        end do
        do j = 0, n - 1
            do i = 0, n - 1
                pos(dir) = side*(n - 1)
                pos(t(1)) = i
                pos(t(2)) = j
                idx(1 + i + n*j) = 1 + pos(1) + n*(pos(2) + n*pos(3))
            end do
        end do
    end subroutine

    subroutine load_faces(mesh, list, shift, faces, status)
        type(tensor_mesh_t), intent(in) :: mesh
        integer, intent(in) :: list(:, :)
        real(dp), intent(in) :: shift(3)
        type(seam_face_t), allocatable, intent(out) :: faces(:)
        integer, intent(out) :: status
        integer, allocatable :: gidx(:), nidx(:)
        integer :: f, e, nq, np, a
        status = 1
        nq = mesh%geometry_degree + 1
        np = mesh%degree + 1
        allocate (faces(size(list, 2)), gidx(nq*nq), nidx(np*np))
        do f = 1, size(list, 2)
            e = list(1, f)
            if (e < 1 .or. e > mesh%nelem) return
            faces(f)%element = e
            faces(f)%dir = list(2, f)
            faces(f)%side = list(3, f)
            call face_nodes(nq, list(2, f), list(3, f), gidx)
            call face_nodes(np, list(2, f), list(3, f), nidx)
            faces(f)%geometry = mesh%geometry(:, gidx, e)
            do a = 1, nq*nq
                faces(f)%geometry(:, a) = faces(f)%geometry(:, a) + shift
            end do
            faces(f)%dof = mesh%dof(nidx, e)
            faces(f)%center = sum(faces(f)%geometry, 2)/real(nq*nq, dp)
            faces(f)%radius = 0
            do a = 1, nq*nq
                faces(f)%radius = max(faces(f)%radius, &
                                      norm2(faces(f)%geometry(:, a) - faces(f)%center))
            end do
        end do
        status = 0
    end subroutine

    subroutine match_faces(mesh, seam, match, matched, status)
        !! match(i) is the master DOF at the image of slave-face DOF i.
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_seam_t), intent(in) :: seam
        integer, allocatable, intent(out) :: match(:)
        logical, intent(out) :: matched
        integer, intent(out) :: status
        type(seam_face_t), allocatable :: master(:), slave(:)
        logical, allocatable :: on_master(:)
        real(dp) :: y(3), scale, tol
        integer :: f, a, m, g, i, k, found
        allocate (match(mesh%ndof), on_master(mesh%ndof))
        match = 0
        on_master = .false.
        matched = .false.
        call load_faces(mesh, seam%slave_faces, [0.0_dp, 0.0_dp, 0.0_dp], slave, &
                        status)
        if (status /= 0) return
        do g = 1, size(seam%images, 2)
            call load_faces(mesh, seam%master_faces, seam%images(:, g), master, &
                            status)
            if (status /= 0) return
            do f = 1, size(master)
                on_master(master(f)%dof) = .true.
            end do
        end do
        scale = maxval(abs(mesh%coordinates)) + 1
        tol = 1.e-10_dp*scale
        do f = 1, size(slave)
            do a = 1, size(slave(f)%dof)
                i = slave(f)%dof(a)
                if (match(i) /= 0) cycle
                if (on_master(i)) then
                    match(i) = i
                    cycle
                end if
                y = matmul(seam%motion(:, 1:3), mesh%coordinates(:, i)) + &
                    seam%motion(:, 4)
                found = 0
                images: do g = 1, size(seam%images, 2)
                    do m = 1, mesh%ndof
                        if (.not. on_master(m)) cycle
                        if (maxval(abs(mesh%coordinates(:, m) + seam%images(:, g) - &
                                       y)) <= tol) then
                            found = m
                            exit images
                        end if
                    end do
                end do images
                if (found == 0) return
                match(i) = found
            end do
        end do
        ! A matching node set is necessary; corner-to-vertex agreement is
        ! checked when the identified mesh is built.
        k = 0
        do i = 1, mesh%ndof
            if (match(i) /= 0 .and. match(i) /= i) k = k + 1
        end do
        matched = k > 0
    end subroutine

    subroutine merge_seam(mesh, match, closed, status)
        type(tensor_mesh_t), intent(in) :: mesh
        integer, intent(in) :: match(:)
        type(tensor_mesh_t), intent(out) :: closed
        integer, intent(out) :: status
        integer, allocatable :: map(:), newid(:), vertex_of(:), vmap(:), vid(:), &
                                gather(:, :), vertices(:, :)
        logical, allocatable :: boundary(:)
        integer :: i, e, c, bits(3), idx, nv, n, k
        status = 1
        n = mesh%degree + 1
        allocate (map(mesh%ndof), newid(mesh%ndof))
        do i = 1, mesh%ndof
            map(i) = i
            if (match(i) /= 0) map(i) = match(i)
        end do
        k = 0
        newid = 0
        do i = 1, mesh%ndof
            if (map(i) /= i) cycle
            k = k + 1
            newid(i) = k
        end do
        nv = maxval(mesh%vertices)
        allocate (vertex_of(mesh%ndof), vmap(nv), vid(nv))
        vertex_of = 0
        do e = 1, mesh%nelem
            do c = 1, 2**mesh%dim
                call corner_index(c, n, bits, idx)
                vertex_of(mesh%dof(idx, e)) = mesh%vertices(c, e)
            end do
        end do
        do i = 1, nv
            vmap(i) = i
        end do
        do e = 1, mesh%nelem
            do c = 1, 2**mesh%dim
                call corner_index(c, n, bits, idx)
                i = mesh%dof(idx, e)
                if (map(i) == i) cycle
                if (vertex_of(map(i)) == 0) return
                vmap(mesh%vertices(c, e)) = vertex_of(map(i))
            end do
        end do
        k = 0
        vid = 0
        do i = 1, nv
            if (vmap(i) /= i) cycle
            k = k + 1
            vid(i) = k
        end do
        allocate (gather(size(mesh%dof, 1), mesh%nelem), &
                  vertices(size(mesh%vertices, 1), mesh%nelem))
        do e = 1, mesh%nelem
            gather(:, e) = newid(map(mesh%dof(:, e)))
            vertices(:, e) = vid(vmap(mesh%vertices(:, e)))
        end do
        allocate (boundary(maxval(newid)))
        do i = 1, mesh%ndof
            if (newid(i) > 0) boundary(newid(i)) = mesh%boundary(i)
        end do
        call build_tensor_mesh(mesh%dim, mesh%degree, mesh%geometry_degree, &
                               vertices, mesh%geometry, closed, status, gather=gather)
        if (status /= 0) return
        closed%boundary = boundary
    end subroutine

    pure subroutine corner_index(c, n, bits, idx)
        integer, intent(in) :: c, n
        integer, intent(out) :: bits(3), idx
        bits(1) = modulo(c - 1, 2)
        bits(2) = modulo((c - 1)/2, 2)
        bits(3) = (c - 1)/4
        idx = 1 + bits(1)*(n - 1) + n*(bits(2)*(n - 1) + n*bits(3)*(n - 1))
    end subroutine

    subroutine face_map(g, zq, a, b, x, xa, xb)
        real(dp), intent(in) :: g(:, :), zq(:), a, b
        real(dp), intent(out) :: x(3), xa(3), xb(3)
        real(dp), allocatable :: l(:, :), dl(:, :)
        integer :: r, s, nq, k
        nq = size(zq)
        call lagrange_table(zq, [a, b], l, dl)
        x = 0; xa = 0; xb = 0
        do s = 1, nq
            do r = 1, nq
                k = r + nq*(s - 1)
                x = x + l(1, r)*l(2, s)*g(:, k)
                xa = xa + dl(1, r)*l(2, s)*g(:, k)
                xb = xb + l(1, r)*dl(2, s)*g(:, k)
            end do
        end do
    end subroutine

    subroutine invert_face(g, zq, y, ab, ok, reach)
        !! Damped Gauss-Newton inverse of a face map, also beyond its reference
        !! square. A limit counts only if the residual is normal to the face,
        !! so a fold of the polynomial extension is never taken for a preimage.
        !! Far from the square a curved face's extension folds and has spurious
        !! preimages, so only preimages within reach (max norm) count; callers
        !! treat a point without one as outside the cell, and cubature nodes
        !! are checked to lie in both cells. See reach_for.
        real(dp), intent(in) :: g(:, :), zq(:), y(3)
        real(dp), intent(inout) :: ab(2)
        logical, intent(out) :: ok
        real(dp), intent(in) :: reach
        real(dp), parameter :: starts(2, 9) = reshape([0.0_dp, 0.0_dp, &
                   -0.8_dp, -0.8_dp, 0.8_dp, -0.8_dp, 0.8_dp, 0.8_dp, -0.8_dp, 0.8_dp, &
              0.0_dp, -0.8_dp, 0.8_dp, 0.0_dp, 0.0_dp, 0.8_dp, -0.8_dp, 0.0_dp], [2, 9])
        real(dp) :: x(3), xa(3), xb(3), m(2, 2), rhs(2), step(2), det, r(3), &
                    normal(3), res, trial, scale, guess(2), c(2)
        integer :: it, start, half, far
        guess = ab
        ok = .false.
        far = 0
        do start = 0, 9
            c = guess
            if (start > 0) c = starts(:, start)
            call face_map(g, zq, c(1), c(2), x, xa, xb)
            res = norm2(x - y)
            do it = 1, 60
                m(1, 1) = dot_product(xa, xa)
                m(1, 2) = dot_product(xa, xb)
                m(2, 2) = dot_product(xb, xb)
                rhs(1) = dot_product(xa, y - x)
                rhs(2) = dot_product(xb, y - x)
                det = m(1, 1)*m(2, 2) - m(1, 2)**2
                if (det <= 0) exit
                step(1) = (m(2, 2)*rhs(1) - m(1, 2)*rhs(2))/det
                step(2) = (m(1, 1)*rhs(2) - m(1, 2)*rhs(1))/det
                do half = 0, 30
                    call face_map(g, zq, c(1) + step(1), c(2) + step(2), x, xa, xb)
                    trial = norm2(x - y)
                    if (trial <= res) exit
                    step = step/2
                end do
                ! A stalled line search is a limit point: test it below.
                if (trial > res) exit
                c = c + step
                res = trial
                if (maxval(abs(c)) > 1.e2_dp) exit
                if (maxval(abs(step)) <= 4*epsilon(1.0_dp)*(1 + maxval(abs(c)))) exit
            end do
            call face_map(g, zq, c(1), c(2), x, xa, xb)
            normal = [xa(2)*xb(3) - xa(3)*xb(2), xa(3)*xb(1) - xa(1)*xb(3), &
                      xa(1)*xb(2) - xa(2)*xb(1)]
            if (norm2(normal) <= 0) cycle
            normal = normal/norm2(normal)
            r = x - y
            r = r - normal*dot_product(normal, r)
            scale = norm2(xa) + norm2(xb)
            if (norm2(r) <= 1.e-12_dp*scale*(1 + maxval(abs(c)))**size(zq)) then
                if (maxval(abs(c)) <= reach) then
                    ab = c
                    ok = .true.
                    return
                end if
                ! Two starts converging beyond reach: the point lies outside.
                far = far + 1
                if (far == 2) return
            end if
        end do
    end subroutine

    pure real(dp) function reach_for(small, large)
        !! Reach in the chart of the larger cell: it must cover the smaller cell
        !! wherever that can meet the square, else the clipper sees a spurious
        !! crossing off its clip line, yet stay below the fold branches of a
        !! curved cell, which lie about two cells beyond its square.
        type(seam_face_t), intent(in) :: small, large
        reach_for = 1.5_dp + 4*small%radius/large%radius
    end function

    subroutine rim_tests(n, x, remove, coef)
        !! psi_i = phi_i + sum_r coef(i,r) phi_r has degree n-1-size(remove).
        !! Leading coefficients of Lagrange phi_j are lambda_j, and the k-th
        !! power-sum condition sum lambda_j x_j^k = 0 removes the top terms.
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:)
        logical, intent(in) :: remove(:)
        real(dp), intent(out) :: coef(:, :)
        real(dp) :: lambda(n), m(2, 2), r(2), det
        integer :: i, j, k, rr(2), nr
        coef = 0
        do i = 1, n
            coef(i, i) = 1
        end do
        nr = 0
        do j = 1, n
            if (remove(j)) then
                nr = nr + 1
                rr(nr) = j
            end if
        end do
        if (nr == 0) return
        do j = 1, n
            lambda(j) = 1
            do k = 1, n
                if (k /= j) lambda(j) = lambda(j)/(x(j) - x(k))
            end do
        end do
        do i = 1, n
            if (remove(i)) then
                coef(i, :) = 0
                cycle
            end if
            if (nr == 1) then
                coef(i, rr(1)) = -lambda(i)/lambda(rr(1))
            else
                do k = 1, 2
                    m(k, 1) = lambda(rr(1))*x(rr(1))**(k - 1)
                    m(k, 2) = lambda(rr(2))*x(rr(2))**(k - 1)
                    r(k) = -lambda(i)*x(i)**(k - 1)
                end do
                det = m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1)
                coef(i, rr(1)) = (m(2, 2)*r(1) - m(1, 2)*r(2))/det
                coef(i, rr(2)) = (m(1, 1)*r(2) - m(2, 1)*r(1))/det
            end if
        end do
    end subroutine

    subroutine build_constraint(mesh, setup, seam, status)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_seam_t), intent(in) :: setup
        type(tensor_seam_t), intent(out) :: seam
        integer, intent(out) :: status
        type(seam_face_t), allocatable :: slave(:), master(:), copies(:)
        type(face_split_t), allocatable :: sparts(:), mparts(:)
        type(seam_cubature_t), allocatable :: cub(:)
        real(dp), allocatable :: zq(:), zp(:), w(:), d(:, :), b(:, :), test(:, :, :, :)
        real(dp), allocatable :: overlap(:), area(:)
        logical, allocatable :: on_master(:), remove(:)
        integer, allocatable :: srow(:), ccol(:), cols(:)
        integer :: n, nf, f, g, i, k, a, ns, nc, nm, side, info
        real(dp) :: defect, small
        status = 1
        seam%method = setup%method
        seam%conforming = .false.
        seam%area_defect = 0
        seam%master_faces = setup%master_faces
        seam%slave_faces = setup%slave_faces
        seam%motion = setup%motion
        seam%images = setup%images
        n = mesh%degree + 1
        call gll_rule(mesh%geometry_degree, zq, w, status)
        if (status /= 0) return
        call gll_rule(mesh%degree, zp, w, status)
        if (status /= 0) return
        call load_faces(mesh, setup%slave_faces, [0.0_dp, 0.0_dp, 0.0_dp], slave, &
                        status)
        if (status /= 0) return
        nm = size(setup%master_faces, 2)
        allocate (master(nm*size(setup%images, 2)))
        do g = 1, size(setup%images, 2)
            call load_faces(mesh, setup%master_faces, setup%images(:, g), copies, &
                            status)
            if (status /= 0) return
            master(1 + nm*(g - 1):nm*g) = copies
        end do
        allocate (on_master(mesh%ndof), srow(mesh%ndof), ccol(mesh%ndof))
        on_master = .false.
        do f = 1, size(master)
            on_master(master(f)%dof) = .true.
        end do
        ! Slave rows: unconstrained slave-face DOFs that are not Dirichlet.
        srow = 0
        ns = 0
        do f = 1, size(slave)
            do a = 1, n*n
                i = slave(f)%dof(a)
                if (mesh%boundary(i) .or. on_master(i) .or. srow(i) /= 0) cycle
                ns = ns + 1
                srow(i) = ns
            end do
        end do
        ! Compact columns: every primary DOF a slave row can depend on.
        ccol = 0
        nc = 0
        allocate (cols(mesh%ndof))
        do f = 1, size(master)
            call add_columns(master(f)%dof)
        end do
        do f = 1, size(slave)
            call add_columns(slave(f)%dof)
        end do
        allocate (b(ns, nc), test(n, n, 2, size(slave)), remove(n))
        b = 0
        ! Rim-modified one-dimensional test bases per slave face and tangent.
        do f = 1, size(slave)
            do k = 1, 2
                remove = .false.
                do side = 0, 1
                    remove(1 + side*(n - 1)) = rim_side(slave(f), k, side)
                end do
                call rim_tests(n, zp, remove, test(:, :, k, f))
            end do
        end do
        nf = size(slave)
        allocate (area(nf), overlap(nf))
        area = 0
        overlap = 0
        if (setup%method == SEAM_INTERPOLATION) then
            call interpolation_rows(status)
            if (status /= 0) return
        else
            allocate (d(ns, ns))
            d = 0
            if (allocated(setup%cubature)) then
                cub = setup%cubature
            else
                allocate (cub(size(slave)), sparts(size(slave)), mparts(size(master)))
            end if
            ! The clipper inverts a slave part's boundary in a master part's
            ! chart, which is reliable only near its square and away from the
            ! folds of a curved cell: master parts are nearly bilinear and
            ! slave parts small against them.
            !$omp parallel do default(shared) private(f) schedule(dynamic)
            do f = 1, merge(0, size(master), allocated(setup%cubature))
                call split_face(master(f), zq, 0.05_dp, mparts(f))
            end do
            !$omp end parallel do
            if (.not. allocated(setup%cubature)) then
                small = huge(1.0_dp)
                do f = 1, size(master)
                    small = min(small, 0.25_dp*minval(mparts(f)%part%radius))
                end do
            end if
            !$omp parallel do default(shared) private(f) schedule(dynamic)
            do f = 1, merge(0, size(slave), allocated(setup%cubature))
                call split_face(slave(f), zq, 0.01_dp, sparts(f), small)
            end do
            !$omp end parallel do
            call mortar_rows(status)
            if (status /= 0) return
            call dense_solve(ns, nc, d, b, info)
            if (info /= 0) then
                status = 3
                return
            end if
            defect = 0
            do f = 1, nf
                if (area(f) > 0) defect = max(defect, abs(overlap(f) - area(f))/area(f))
            end do
            seam%area_defect = defect
            call move_alloc(cub, seam%cubature)
        end if
        call finish_constraint(mesh, srow, ns, cols, nc, b, seam)
        status = 0
    contains
        subroutine add_columns(list)
            integer, intent(in) :: list(:)
            integer :: j
            do j = 1, size(list)
                if (srow(list(j)) /= 0 .or. ccol(list(j)) /= 0) cycle
                nc = nc + 1
                ccol(list(j)) = nc
                cols(nc) = list(j)
            end do
        end subroutine

        logical function rim_side(face, k, side)
            type(seam_face_t), intent(in) :: face
            integer, intent(in) :: k, side
            integer :: j, idx
            rim_side = .true.
            do j = 0, n - 1
                if (k == 1) then
                    idx = 1 + side*(n - 1) + n*j
                else
                    idx = 1 + j + n*side*(n - 1)
                end if
                if (.not. mesh%boundary(face%dof(idx))) rim_side = .false.
            end do
        end function

        subroutine interpolation_rows(info)
            integer, intent(out) :: info
            real(dp), allocatable :: phi(:)
            real(dp) :: y(3), ab(2), excess, best, bab(2)
            integer :: fs, a, i, m, bm
            logical :: ok
            logical, allocatable :: done(:)
            allocate (phi(n*n), done(ns))
            done = .false.
            info = 0
            do fs = 1, size(slave)
                do a = 1, n*n
                    i = slave(fs)%dof(a)
                    if (srow(i) == 0) cycle
                    if (done(srow(i))) cycle
                    done(srow(i)) = .true.
                    y = matmul(setup%motion(:, 1:3), mesh%coordinates(:, i)) + &
                        setup%motion(:, 4)
                    best = huge(1.0_dp)
                    bm = 0
                    do m = 1, size(master)
                        if (norm2(y - master(m)%center) > 1.5_dp*master(m)%radius) &
                            cycle
                        ab = 0
                        call invert_face(master(m)%geometry, zq, y, ab, ok, 1.5_dp)
                        if (.not. ok) cycle
                        excess = maxval(abs(ab)) - 1
                        if (excess < best) then
                            best = excess
                            bm = m
                            bab = ab
                        end if
                        if (excess <= 1.e-10_dp) exit
                    end do
                    if (bm == 0 .or. best > 1.e-6_dp) then
                        info = 4
                        return
                    end if
                    call face_basis(zp, bab, phi)
                    do m = 1, n*n
                        b(srow(i), ccol(master(bm)%dof(m))) = &
                            b(srow(i), ccol(master(bm)%dof(m))) + phi(m)
                    end do
                end do
            end do
        end subroutine

        subroutine mortar_rows(info)
            integer, intent(out) :: info
            integer :: fs, fail, cell
            real(dp), allocatable :: dl(:, :), bl(:, :, :)
            integer, allocatable :: hits(:)
            info = 0
            fail = 0
            !$omp parallel do default(shared) private(fs,dl,bl,hits,cell) &
            !$omp schedule(dynamic) reduction(max:fail)
            do fs = 1, size(slave)
                cell = 0
                if (.not. allocated(setup%cubature)) call overlap_points(fs, cell)
                fail = max(fail, cell)
                if (cell /= 0) cycle
                call slave_cell(fs, dl, bl, hits)
                !$omp critical (seam_scatter)
                call scatter_cell(fs, dl, bl, hits)
                !$omp end critical (seam_scatter)
            end do
            !$omp end parallel do
            info = fail
        end subroutine

        subroutine slave_cell(fs, dl, bl, hits)
            !! Local slave mass and the cross mass with every overlapping master.
            integer, intent(in) :: fs
            real(dp), allocatable, intent(out) :: dl(:, :), bl(:, :, :)
            integer, allocatable, intent(out) :: hits(:)
            real(dp), allocatable :: xg(:), wg(:), psi(:), phi(:), work(:, :, :)
            integer, allocatable :: slot(:)
            real(dp) :: x(3), xa(3), xb(3), jac
            integer :: i, j, m, q, nh, st
            allocate (psi(n*n), phi(n*n), dl(n*n, n*n), work(n*n, n*n, size(master)), &
                      hits(size(master)), slot(size(master)))
            dl = 0
            nh = 0
            slot = 0
            call gauss_rule(mesh%degree + mesh%geometry_degree + 2, xg, wg, st)
            do j = 1, size(xg)
                do i = 1, size(xg)
                    call face_map(slave(fs)%geometry, zq, xg(i), xg(j), x, xa, xb)
                    call generated_tensor_surface_jacobian(xa(1), xa(2), xa(3), &
                                                           xb(1), xb(2), xb(3), jac)
                    call test_basis(fs, [xg(i), xg(j)], psi)
                    call face_basis(zp, [xg(i), xg(j)], phi)
                    dl = dl + wg(i)*wg(j)*jac*spread(psi, 2, n*n)*spread(phi, 1, n*n)
                    area(fs) = area(fs) + wg(i)*wg(j)*jac
                end do
            end do
            associate (cb => cub(fs))
                do q = 1, size(cb%weight)
                    m = cb%master(q)
                    if (slot(m) == 0) then
                        nh = nh + 1
                        hits(nh) = m
                        slot(m) = nh
                        work(:, :, nh) = 0
                    end if
                    call test_basis(fs, cb%slave_ab(:, q), psi)
                    call face_basis(zp, cb%master_ab(:, q), phi)
                    work(:, :, slot(m)) = work(:, :, slot(m)) + cb%weight(q)* &
                                          spread(psi, 2, n*n)*spread(phi, 1, n*n)
                    overlap(fs) = overlap(fs) + cb%weight(q)
                end do
            end associate
            bl = work(:, :, 1:nh)
            hits = hits(1:nh)
        end subroutine

        subroutine overlap_points(fs, fail)
            !! Overlap cubature of slave face fs with every master face.
            integer, intent(in) :: fs
            integer, intent(inout) :: fail
            real(dp), allocatable :: pts(:, :), wts(:), sq(:, :), mq(:, :), wq(:), &
                                     t2(:, :), t1(:)
            integer, allocatable :: owner(:), ti(:)
            real(dp) :: x(3), xa(3), xb(3), y(3), ab(2), jac, c(3), cs(3), sab(2), &
                        weight, inverse(3, 4)
            integer :: m, q, np, st, sp, mp
            logical :: ok, swap, rim(4), open(4)
            real(dp), parameter :: edge = 1.e-12_dp
            allocate (sq(2, 64), mq(2, 64), wq(64), owner(64))
            np = 0
            inverse(:, 1:3) = transpose(setup%motion(:, 1:3))
            inverse(:, 4) = -matmul(inverse(:, 1:3), setup%motion(:, 4))
            c = matmul(setup%motion(:, 1:3), slave(fs)%center) + setup%motion(:, 4)
            do m = 1, size(master)
                if (norm2(c - master(m)%center) > &
                    1.5_dp*(slave(fs)%radius + master(m)%radius)) cycle
                ! The discrete rims of both faces differ by the geometry error;
                ! slave slivers beyond a master rim see its polynomial extension.
                rim = [rim_side(master(m), 1, 1), rim_side(master(m), 1, 0), &
                       rim_side(master(m), 2, 1), rim_side(master(m), 2, 0)]
                do sp = 1, size(sparts(fs)%part)
                    associate (s => sparts(fs)%part(sp))
                        cs = matmul(setup%motion(:, 1:3), s%center) + setup%motion(:, 4)
                        do mp = 1, size(mparts(m)%part)
                            associate (t => mparts(m)%part(mp))
                                if (norm2(cs - t%center) > 1.5_dp*(s%radius + &
                                                                   t%radius)) cycle
                                ! Only part sides on the parent's rim are open.
                                open = rim .and. [t%origin(1) + t%scale >= 1 - edge, &
                                                  t%origin(1) - t%scale <= edge - 1, &
                                                  t%origin(2) + t%scale >= 1 - edge, &
                                                  t%origin(2) - t%scale <= edge - 1]
                                if (.not. cells_overlap(s, t, setup%motion, zq, &
                                                        open)) cycle
                                ! Clip the smaller part in the chart of the larger
                                ! one, whose extension then stays near its square.
                                swap = t%radius > s%radius .or. any(open)
                                if (swap) then
                                    call overlap_cubature(t, s, inverse, zq, &
                                             2*mesh%degree + mesh%geometry_degree + 4, &
                                                          pts, wts, st, open=open)
                                else
                                    call overlap_cubature(s, t, setup%motion, zq, &
                                             2*mesh%degree + mesh%geometry_degree + 4, &
                                                          pts, wts, st)
                                end if
                                if (st /= 0) then
                                    fail = 2
                                    return
                                end if
                                do q = 1, size(wts)
                                    if (swap) then
                                        ab = pts(:, q)
                                        call face_map(t%geometry, zq, ab(1), ab(2), &
                                                      y, xa, xb)
                                        call generated_tensor_surface_jacobian( &
                                            xa(1), xa(2), xa(3), xb(1), xb(2), xb(3), &
                                            jac)
                                        weight = wts(q)*jac
                                        ab = t%origin + t%scale*ab
                                        x = matmul(inverse(:, 1:3), y) + inverse(:, 4)
                                        sab = 0
                                        call invert_face(s%geometry, zq, x, sab, ok, &
                                                         1.5_dp)
                                        ! Nodes lie in the overlap, hence in s.
                                        if (maxval(abs(sab)) > 1 + 1.e-8_dp) &
                                            ok = .false.
                                        sab = s%origin + s%scale*sab
                                    else
                                        sab = s%origin + s%scale*pts(:, q)
                                        call face_map(slave(fs)%geometry, zq, sab(1), &
                                                      sab(2), x, xa, xb)
                                        call generated_tensor_surface_jacobian( &
                                            xa(1), xa(2), xa(3), xb(1), xb(2), xb(3), &
                                            jac)
                                        weight = wts(q)*s%scale**2*jac
                                        y = matmul(setup%motion(:, 1:3), x) + &
                                            setup%motion(:, 4)
                                        ab = 0
                                        call invert_face(t%geometry, zq, y, ab, ok, &
                                                         1.5_dp)
                                        if (maxval(abs(ab)) > 1 + 1.e-8_dp) ok = .false.
                                        ab = t%origin + t%scale*ab
                                    end if
                                    if (.not. ok) then
                                        fail = 2
                                        return
                                    end if
                                    if (np == size(wq)) then
                                        allocate (t2(2, 2*np), t1(2*np), ti(2*np))
                                        t2(:, 1:np) = sq
                                        call move_alloc(t2, sq)
                                        allocate (t2(2, 2*np))
                                        t2(:, 1:np) = mq
                                        call move_alloc(t2, mq)
                                        t1(1:np) = wq
                                        call move_alloc(t1, wq)
                                        ti(1:np) = owner
                                        call move_alloc(ti, owner)
                                    end if
                                    np = np + 1
                                    sq(:, np) = sab
                                    mq(:, np) = ab
                                    wq(np) = weight
                                    owner(np) = m
                                end do
                            end associate
                        end do
                    end associate
                end do
            end do
            cub(fs)%master = owner(1:np)
            cub(fs)%slave_ab = sq(:, 1:np)
            cub(fs)%master_ab = mq(:, 1:np)
            cub(fs)%weight = wq(1:np)
        end subroutine

        subroutine scatter_cell(fs, dl, bl, hits)
            integer, intent(in) :: fs, hits(:)
            real(dp), intent(in) :: dl(:, :), bl(:, :, :)
            integer :: a, c, h, r, j
            do a = 1, n*n
                r = srow(slave(fs)%dof(a))
                if (r == 0) cycle
                do c = 1, n*n
                    j = slave(fs)%dof(c)
                    if (srow(j) /= 0) then
                        d(r, srow(j)) = d(r, srow(j)) + dl(a, c)
                    else
                        b(r, ccol(j)) = b(r, ccol(j)) - dl(a, c)
                    end if
                end do
                do h = 1, size(hits)
                    do c = 1, n*n
                        j = ccol(master(hits(h))%dof(c))
                        b(r, j) = b(r, j) + bl(a, c, h)
                    end do
                end do
            end do
        end subroutine

        subroutine test_basis(fs, ab, psi)
            integer, intent(in) :: fs
            real(dp), intent(in) :: ab(2)
            real(dp), intent(out) :: psi(:)
            real(dp), allocatable :: l(:, :), dl(:, :)
            real(dp) :: ta(n), tb(n)
            integer :: i, j
            call lagrange_table(zp, ab, l, dl)
            ta = matmul(test(:, :, 1, fs), l(1, :))
            tb = matmul(test(:, :, 2, fs), l(2, :))
            do j = 1, n
                do i = 1, n
                    psi(i + n*(j - 1)) = ta(i)*tb(j)
                end do
            end do
        end subroutine
    end subroutine

    subroutine face_basis(zp, ab, phi)
        real(dp), intent(in) :: zp(:), ab(2)
        real(dp), intent(out) :: phi(:)
        real(dp), allocatable :: l(:, :), dl(:, :)
        integer :: i, j, n
        n = size(zp)
        call lagrange_table(zp, ab, l, dl)
        do j = 1, n
            do i = 1, n
                phi(i + n*(j - 1)) = l(1, i)*l(2, j)
            end do
        end do
    end subroutine

    subroutine finish_constraint(mesh, srow, ns, cols, nc, c, seam)
        type(tensor_mesh_t), intent(in) :: mesh
        integer, intent(in) :: srow(:), ns, cols(:), nc
        real(dp), intent(in) :: c(:, :)
        type(tensor_seam_t), intent(inout) :: seam
        integer :: i, k, r, j, nz
        integer, allocatable :: fill(:)
        seam%nfull = mesh%ndof
        seam%nslave = ns
        allocate (seam%reduced(mesh%ndof), seam%slave(ns))
        k = 0
        do i = 1, mesh%ndof
            if (srow(i) /= 0) then
                seam%reduced(i) = -srow(i)
                seam%slave(srow(i)) = i
            else
                k = k + 1
                seam%reduced(i) = k
            end if
        end do
        seam%nreduced = k
        allocate (seam%primary(k), seam%boundary(k))
        do i = 1, mesh%ndof
            if (seam%reduced(i) > 0) then
                seam%primary(seam%reduced(i)) = i
                seam%boundary(seam%reduced(i)) = mesh%boundary(i)
            end if
        end do
        nz = count(abs(c(1:ns, 1:nc)) > 0)
        allocate (seam%ptr(ns + 1), seam%col(nz), seam%val(nz))
        seam%ptr(1) = 1
        nz = 0
        do r = 1, ns
            do j = 1, nc
                if (.not. abs(c(r, j)) > 0) cycle
                nz = nz + 1
                seam%col(nz) = seam%reduced(cols(j))
                seam%val(nz) = c(r, j)
            end do
            seam%ptr(r + 1) = nz + 1
        end do
        allocate (seam%tptr(k + 1), seam%trow(nz), seam%tval(nz), fill(k))
        seam%tptr = 0
        do i = 1, nz
            seam%tptr(seam%col(i) + 1) = seam%tptr(seam%col(i) + 1) + 1
        end do
        seam%tptr(1) = 1
        do i = 1, k
            seam%tptr(i + 1) = seam%tptr(i + 1) + seam%tptr(i)
        end do
        fill = seam%tptr(1:k)
        do r = 1, ns
            do i = seam%ptr(r), seam%ptr(r + 1) - 1
                j = seam%col(i)
                seam%trow(fill(j)) = r
                seam%tval(fill(j)) = seam%val(i)
                fill(j) = fill(j) + 1
            end do
        end do
    end subroutine

    subroutine dense_solve(n, m, a, b, info)
        !! Overwrite b with a^{-1} b by partial-pivoting LU (setup only).
        integer, intent(in) :: n, m
        real(dp), intent(inout) :: a(:, :), b(:, :)
        integer, intent(out) :: info
        integer :: i, j, k, piv
        real(dp) :: tmp
        real(dp), allocatable :: row(:)
        info = 0
        allocate (row(max(n, m)))
        do k = 1, n
            piv = k - 1 + maxloc(abs(a(k:n, k)), 1)
            if (abs(a(piv, k)) <= tiny(1.0_dp)) then
                info = k
                return
            end if
            if (piv /= k) then
                row(1:n) = a(k, 1:n); a(k, 1:n) = a(piv, 1:n); a(piv, 1:n) = row(1:n)
                row(1:m) = b(k, 1:m); b(k, 1:m) = b(piv, 1:m); b(piv, 1:m) = row(1:m)
            end if
            a(k + 1:n, k) = a(k + 1:n, k)/a(k, k)
            !$omp parallel do default(shared) private(j) schedule(static)
            do j = k + 1, n
                a(k + 1:n, j) = a(k + 1:n, j) - a(k + 1:n, k)*a(k, j)
            end do
            !$omp end parallel do
        end do
        !$omp parallel do default(shared) private(j,i,tmp) schedule(static)
        do j = 1, m
            do i = 2, n
                tmp = dot_product(a(i, 1:i - 1), b(1:i - 1, j))
                b(i, j) = b(i, j) - tmp
            end do
            do i = n, 1, -1
                tmp = 0
                if (i < n) tmp = dot_product(a(i, i + 1:n), b(i + 1:n, j))
                b(i, j) = (b(i, j) - tmp)/a(i, i)
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine split_face(face, zq, limit, split, radius)
        !! Halve a curved face until every part deviates from its bilinear
        !! corner interpolant by at most limit times its diameter; the chart
        !! of a nearly affine part is then injective well beyond its square.
        !! Parts are also halved until their half diagonals are at most radius.
        type(seam_face_t), intent(in) :: face
        real(dp), intent(in) :: zq(:), limit
        type(face_split_t), intent(out) :: split
        real(dp), intent(in), optional :: radius
        real(dp) :: x(3), xa(3), xb(3), corner(3, 4), a0, b0, h, worst, dev, diam, &
                    u, v, largest
        integer :: k, i, j, r, t, nq, part, c
        nq = size(zq)
        k = 1
        do
            worst = 0
            largest = 0
            do j = 0, k - 1
                do i = 0, k - 1
                    h = 1.0_dp/k
                    a0 = -1 + (2*i + 1)*h
                    b0 = -1 + (2*j + 1)*h
                    do c = 1, 4
                        call face_map(face%geometry, zq, a0 + h*merge(-1, 1, c == 1 &
                                            .or. c == 4), b0 + h*merge(-1, 1, c <= 2), &
                                      corner(:, c), xa, xb)
                    end do
                    diam = max(norm2(corner(:, 1) - corner(:, 3)), &
                               norm2(corner(:, 2) - corner(:, 4)))
                    largest = max(largest, diam/2)
                    do t = 0, 4
                        do r = 0, 4
                            u = -1 + r/2.0_dp
                            v = -1 + t/2.0_dp
                            call face_map(face%geometry, zq, a0 + h*u, b0 + h*v, x, &
                                          xa, xb)
                            dev = norm2(x - ((1 - u)*(1 - v)*corner(:, 1) + &
                                             (1 + u)*(1 - v)*corner(:, 2) + &
                                             (1 + u)*(1 + v)*corner(:, 3) + &
                                             (1 - u)*(1 + v)*corner(:, 4))/4)
                            worst = max(worst, dev/diam)
                        end do
                    end do
                end do
            end do
            if (worst <= limit .or. k >= 32) then
                if (.not. present(radius)) exit
                if (largest <= radius .or. k >= 128) exit
            end if
            k = 2*k
        end do
        allocate (split%part(k*k))
        part = 0
        do j = 0, k - 1
            do i = 0, k - 1
                part = part + 1
                h = 1.0_dp/k
                split%part(part)%origin = [-1 + (2*i + 1)*h, -1 + (2*j + 1)*h]
                split%part(part)%scale = h
                allocate (split%part(part)%geometry(3, nq*nq))
                do t = 1, nq
                    do r = 1, nq
                        call face_map(face%geometry, zq, split%part(part)%origin(1) + &
                                      h*zq(r), split%part(part)%origin(2) + h*zq(t), &
                                      x, xa, xb)
                        split%part(part)%geometry(:, r + nq*(t - 1)) = x
                    end do
                end do
                split%part(part)%center = sum(split%part(part)%geometry, 2)/(nq*nq)
                split%part(part)%radius = 0
                do r = 1, nq*nq
                    split%part(part)%radius = max(split%part(part)%radius, &
                       norm2(split%part(part)%geometry(:, r) - split%part(part)%center))
                end do
            end do
        end do
    end subroutine

    logical function cells_overlap(slave, master, motion, zq, open)
        !! Sampled test whether two face cells overlap under the seam map: the
        !! smaller cell is sampled in the chart of the larger one, where it is
        !! nearly a parallelogram, and its padded bounding box meets the square.
        !! open(c) drops side c of a master chart, as in overlap_cubature.
        type(seam_face_t), intent(in) :: slave, master
        real(dp), intent(in) :: motion(3, 4), zq(:)
        logical, intent(in), optional :: open (4)
        real(dp) :: x(3), xa(3), xb(3), y(3), ab(2), lo(2), hi(2), pad
        integer :: i, j, m
        logical :: ok, large_master
        m = 2*size(zq)
        large_master = master%radius >= slave%radius
        cells_overlap = .false.
        ab = 0
        ! A much smaller cell whose centre is beyond reach cannot overlap.
        if (large_master .and. slave%radius < 0.25_dp*master%radius) then
            y = matmul(motion(:, 1:3), slave%center) + motion(:, 4)
            call invert_face(master%geometry, zq, y, ab, ok, reach_for(slave, master))
            if (.not. ok) return
        else if (master%radius < 0.25_dp*slave%radius) then
            x = matmul(transpose(motion(:, 1:3)), master%center - motion(:, 4))
            call invert_face(slave%geometry, zq, x, ab, ok, reach_for(master, slave))
            if (.not. ok) return
        end if
        lo = huge(1.0_dp)
        hi = -huge(1.0_dp)
        do j = 0, m
            do i = 0, m
                if (large_master) then
                    call face_map(slave%geometry, zq, -1 + 2*real(i, dp)/m, &
                                  -1 + 2*real(j, dp)/m, x, xa, xb)
                    y = matmul(motion(:, 1:3), x) + motion(:, 4)
                    call invert_face(master%geometry, zq, y, ab, ok, &
                                     reach_for(slave, master))
                else
                    call face_map(master%geometry, zq, -1 + 2*real(i, dp)/m, &
                                  -1 + 2*real(j, dp)/m, y, xa, xb)
                    x = matmul(transpose(motion(:, 1:3)), y - motion(:, 4))
                    call invert_face(slave%geometry, zq, x, ab, ok, &
                                     reach_for(master, slave))
                end if
                ! A sample without a preimage near the square lies outside.
                if (.not. ok) then
                    ab = 0
                    cycle
                end if
                lo = min(lo, ab)
                hi = max(hi, ab)
            end do
        end do
        cells_overlap = .false.
        if (any(hi < lo)) return
        pad = 0.05_dp*maxval(hi - lo) + 1.e-9_dp
        if (large_master .and. present(open)) then
            where (open ([1, 3])) lo = -huge(1.0_dp)
            where (open ([2, 4])) hi = huge(1.0_dp)
        end if
        cells_overlap = all(lo - pad <= 1) .and. all(hi + pad >= -1)
    end function

    ! ---------------------------------------------------------------------
    ! Overlap cubature of slave cell S and master cell M in slave reference
    ! coordinates: the image of M's boundary is clipped against [-1,1]^2, and
    ! every boundary piece g(u) spans a curved triangle with a centre c,
    ! integral f = sum of integral f(c + t (g - c)) t (g - c) x g' dt du.
    ! The signed fan is exact for any c; for the nearly convex overlaps of a
    ! small cell with a square its nodes stay inside the polygon.
    ! ---------------------------------------------------------------------
    subroutine overlap_cubature(slave, master, motion, zq, order, pts, wts, status, &
                                open)
        !! open(c): half-plane c of the chart cell is not clipped, so the cell
        !! extends polynomially across that edge, e.g. a Dirichlet rim.
        type(seam_face_t), intent(in) :: slave, master
        real(dp), intent(in) :: motion(3, 4), zq(:)
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: pts(:, :), wts(:)
        integer, intent(out) :: status
        logical, intent(in), optional :: open (4)
        integer, allocatable :: kind(:)
        real(dp), allocatable :: piece(:, :), xg(:), wg(:)
        real(dp) :: p0(2), dp0(2), centre(2), ab(2), dab(2), total, t
        integer :: clip, i, j, k, nq, np
        logical :: ok
        status = 0
        allocate (kind(4), piece(4, 4))
        do k = 1, 4
            kind(k) = 1
            piece(:, k) = [real(k, dp), -1.0_dp, 1.0_dp, 0.0_dp]
        end do
        do clip = 1, 4
            if (present(open)) then
                if (open (clip)) cycle
            end if
            call clip_polygon(clip, kind, piece, ok)
            if (.not. ok) then
                status = 2
                return
            end if
            if (size(kind) == 0) exit
        end do
        allocate (pts(2, 0), wts(0))
        if (size(kind) == 0) return
        centre = 0
        do k = 1, size(kind)
            call piece_point(kind(k), piece(:, k), 0.0_dp, p0, dp0, ok)
            if (.not. ok) then
                status = 2
                return
            end if
            centre = centre + p0/size(kind)
        end do
        call gauss_rule(order, xg, wg, i)
        nq = size(xg)
        np = 0
        deallocate (pts, wts)
        allocate (pts(2, size(kind)*nq*nq), wts(size(kind)*nq*nq))
        do k = 1, size(kind)
            do i = 1, nq
                call piece_point(kind(k), piece(:, k), (xg(i) + 1)/2, ab, dab, ok)
                if (.not. ok) then
                    status = 2
                    return
                end if
                do j = 1, nq
                    np = np + 1
                    t = (xg(j) + 1)/2
                    pts(:, np) = centre + t*(ab - centre)
                    wts(np) = wg(i)/2*wg(j)/2*t*((ab(1) - centre(1))*dab(2) - &
                                                 (ab(2) - centre(2))*dab(1))
                end do
            end do
        end do
        pts = pts(:, 1:np)
        wts = wts(1:np)
        total = sum(wts)
        ! A reversed parametrisation orientation flips every fan weight.
        if (total < 0) wts = -wts
        if (abs(total) <= 1.e-14_dp) then
            deallocate (pts, wts)
            allocate (pts(2, 0), wts(0))
        end if
    contains
        subroutine piece_point(kd, dat, u, ab, dab, ok)
            !! kd 1: master edge dat(1) over t in [dat(2),dat(3)];
            !! kd 2: straight segment from dat(1:2) to dat(3:4).
            integer, intent(in) :: kd
            real(dp), intent(in) :: dat(4), u
            real(dp), intent(out) :: ab(2), dab(2)
            logical, intent(out) :: ok
            real(dp) :: t, m(2), dm(2), x(3), xa(3), xb(3), xs(3), dx(3), &
                        g(2, 2), rhs(2), det
            if (kd == 2) then
                ab = dat(1:2) + u*(dat(3:4) - dat(1:2))
                dab = dat(3:4) - dat(1:2)
                ok = .true.
                return
            end if
            t = dat(2) + u*(dat(3) - dat(2))
            select case (nint(dat(1)))
            case (1)
                m = [t, -1.0_dp]; dm = [1.0_dp, 0.0_dp]
            case (2)
                m = [1.0_dp, t]; dm = [0.0_dp, 1.0_dp]
            case (3)
                m = [-t, 1.0_dp]; dm = [-1.0_dp, 0.0_dp]
            case default
                m = [-1.0_dp, -t]; dm = [0.0_dp, -1.0_dp]
            end select
            call face_map(master%geometry, zq, m(1), m(2), x, xa, xb)
            dx = (xa*dm(1) + xb*dm(2))*(dat(3) - dat(2))
            xs = matmul(transpose(motion(:, 1:3)), x - motion(:, 4))
            dx = matmul(transpose(motion(:, 1:3)), dx)
            ! Always from the centre, so a point keeps one preimage throughout.
            ab = 0
            call invert_face(slave%geometry, zq, xs, ab, ok, reach_for(master, slave))
            if (.not. ok) return
            call face_map(slave%geometry, zq, ab(1), ab(2), x, xa, xb)
            g(1, 1) = dot_product(xa, xa)
            g(1, 2) = dot_product(xa, xb)
            g(2, 2) = dot_product(xb, xb)
            rhs = [dot_product(xa, dx), dot_product(xb, dx)]
            det = g(1, 1)*g(2, 2) - g(1, 2)**2
            dab(1) = (g(2, 2)*rhs(1) - g(1, 2)*rhs(2))/det
            dab(2) = (g(1, 1)*rhs(2) - g(1, 2)*rhs(1))/det
        end subroutine

        real(dp) function level(c, ab)
            !! Half-plane c of the reference square: level <= 0 is inside.
            integer, intent(in) :: c
            real(dp), intent(in) :: ab(2)
            select case (c)
            case (1)
                level = ab(1) - 1
            case (2)
                level = -ab(1) - 1
            case (3)
                level = ab(2) - 1
            case default
                level = -ab(2) - 1
            end select
            level = level - 1.e-13_dp
        end function

        real(dp) function signed_level(c, kd, dat, u)
            !! A master boundary point without a slave-chart preimage lies
            !! far outside the slave cell and counts as outside.
            integer, intent(in) :: c, kd
            real(dp), intent(in) :: dat(4), u
            real(dp) :: ab(2), dab(2)
            logical :: okp
            call piece_point(kd, dat, u, ab, dab, okp)
            signed_level = 1
            if (okp) signed_level = level(c, ab)
        end function

        subroutine sub_piece(kd, dat, u0, u1, out)
            integer, intent(in) :: kd
            real(dp), intent(in) :: dat(4), u0, u1
            real(dp), intent(out) :: out(4)
            real(dp) :: a0(2), a1(2), d(2)
            logical :: ok
            if (kd == 1) then
                out = [dat(1), dat(2) + u0*(dat(3) - dat(2)), &
                       dat(2) + u1*(dat(3) - dat(2)), 0.0_dp]
            else
                call piece_point(kd, dat, u0, a0, d, ok)
                call piece_point(kd, dat, u1, a1, d, ok)
                out = [a0, a1]
            end if
        end subroutine

        subroutine clip_polygon(c, kind, piece, ok)
            !! Sutherland-Hodgman with curved pieces; bridges along the clip
            !! line are straight and cancel in the loop integral if doubled.
            integer, intent(in) :: c
            integer, allocatable, intent(inout) :: kind(:)
            real(dp), allocatable, intent(inout) :: piece(:, :)
            logical, intent(out) :: ok
            integer, parameter :: ns = 8
            integer, allocatable :: sk(:), ok_kind(:)
            real(dp), allocatable :: sp(:, :), op(:, :)
            logical, allocatable :: inside(:)
            real(dp) :: u(0:ns), g(0:ns), roots(ns + 2), ab(2), dab(2), ua, ub, &
                        ga, gb, um, gm, exit_point(2), entry_point(2), dat(4)
            integer :: k, j, nr, it, nsub, start, idx, kept
            logical :: pending, any_in, okp
            ok = .true.
            allocate (sk(0), sp(4, 0), inside(0))
            do k = 1, size(kind)
                do j = 0, ns
                    u(j) = real(j, dp)/ns
                    g(j) = signed_level(c, kind(k), piece(:, k), u(j))
                end do
                nr = 1
                roots(1) = 0
                do j = 1, ns
                    if ((g(j - 1) <= 0) .eqv. (g(j) <= 0)) cycle
                    ua = u(j - 1); ub = u(j); ga = g(j - 1); gb = g(j)
                    ! Illinois false position, every fourth step a bisection;
                    ! the inside/outside bracket is kept as in plain bisection.
                    kept = 0
                    do it = 1, 200
                        um = (ua + ub)/2
                        if (mod(it, 4) /= 0 .and. abs(gb - ga) > 0) then
                            um = (ua*gb - ub*ga)/(gb - ga)
                            if (.not. (um > ua .and. um < ub)) um = (ua + ub)/2
                        end if
                        gm = signed_level(c, kind(k), piece(:, k), um)
                        if (abs(gm) <= 4*epsilon(1.0_dp)) then
                            ua = um
                            ub = um
                            exit
                        end if
                        if ((gm <= 0) .eqv. (ga <= 0)) then
                            ua = um; ga = gm
                            if (kept == 1) gb = gb/2
                            kept = 1
                        else
                            ub = um; gb = gm
                            if (kept == -1) ga = ga/2
                            kept = -1
                        end if
                        if (ub - ua <= 4*epsilon(1.0_dp)) exit
                    end do
                    nr = nr + 1
                    roots(nr) = (ua + ub)/2
                end do
                nr = nr + 1
                roots(nr) = 1
                do j = 1, nr - 1
                    if (roots(j + 1) <= roots(j)) cycle
                    um = (roots(j) + roots(j + 1))/2
                    call sub_piece(kind(k), piece(:, k), roots(j), roots(j + 1), dat)
                    sk = [sk, kind(k)]
                    sp = reshape([sp, dat], [4, size(sk)])
                    inside = [inside, signed_level(c, kind(k), piece(:, k), um) <= 0]
                end do
            end do
            nsub = size(sk)
            any_in = any(inside)
            if (.not. any_in) then
                deallocate (kind, piece)
                allocate (kind(0), piece(4, 0))
                return
            end if
            if (all(inside)) then
                kind = sk
                piece = sp
                return
            end if
            ! Start at an inside piece that follows an outside piece.
            start = 1
            do k = 1, nsub
                j = modulo(k - 2, nsub) + 1
                if (inside(k) .and. .not. inside(j)) then
                    start = k
                    exit
                end if
            end do
            allocate (ok_kind(0), op(4, 0))
            pending = .false.
            call piece_point(sk(start), sp(:, start), 0.0_dp, entry_point, dab, okp)
            exit_point = entry_point
            do j = 0, nsub - 1
                idx = modulo(start - 1 + j, nsub) + 1
                if (inside(idx)) then
                    if (pending) then
                        call piece_point(sk(idx), sp(:, idx), 0.0_dp, ab, dab, okp)
                        ok_kind = [ok_kind, 2]
                        op = reshape([op, exit_point, ab], [4, size(ok_kind)])
                        pending = .false.
                    end if
                    ok_kind = [ok_kind, sk(idx)]
                    op = reshape([op, sp(:, idx)], [4, size(ok_kind)])
                else if (.not. pending) then
                    call piece_point(sk(idx), sp(:, idx), 0.0_dp, exit_point, dab, okp)
                    pending = .true.
                end if
            end do
            if (pending) then
                ok_kind = [ok_kind, 2]
                op = reshape([op, exit_point, entry_point], [4, size(ok_kind)])
            end if
            call move_alloc(ok_kind, kind)
            call move_alloc(op, piece)
        end subroutine
    end subroutine

    ! ---------------------------------------------------------------------
    ! Constrained vectors and operators.
    ! ---------------------------------------------------------------------
    subroutine seam_expand(seam, reduced, full)
        !! full = P reduced.
        type(tensor_seam_t), intent(in) :: seam
        real(dp), intent(in) :: reduced(:)
        real(dp), intent(out) :: full(:)
        integer :: k, r, i
        real(dp) :: s
        !$omp parallel default(shared) private(k,r,i,s)
        !$omp do schedule(static)
        do k = 1, seam%nreduced
            full(seam%primary(k)) = reduced(k)
        end do
        !$omp end do nowait
        !$omp do schedule(static)
        do r = 1, seam%nslave
            s = 0
            do i = seam%ptr(r), seam%ptr(r + 1) - 1
                s = s + seam%val(i)*reduced(seam%col(i))
            end do
            full(seam%slave(r)) = s
        end do
        !$omp end do
        !$omp end parallel
    end subroutine

    subroutine seam_fold(seam, full, reduced)
        !! reduced = P^T full.
        type(tensor_seam_t), intent(in) :: seam
        real(dp), intent(in) :: full(:)
        real(dp), intent(out) :: reduced(:)
        integer :: k, i
        real(dp) :: s
        !$omp parallel do default(shared) private(k,i,s) schedule(static)
        do k = 1, seam%nreduced
            s = full(seam%primary(k))
            do i = seam%tptr(k), seam%tptr(k + 1) - 1
                s = s + seam%tval(i)*full(seam%slave(seam%trow(i)))
            end do
            reduced(k) = s
        end do
        !$omp end parallel do
    end subroutine

    subroutine apply_seam_operator(op, seam, x, y, eliminate)
        !! y = P^T A P x on reduced DOFs; Dirichlet rows as in the full operator.
        type(tensor_operator_t), intent(in) :: op
        type(tensor_seam_t), intent(in) :: seam
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        logical, intent(in), optional :: eliminate
        real(dp), allocatable :: xm(:), xf(:), yf(:)
        logical :: mask
        mask = .true.
        if (present(eliminate)) mask = eliminate
        if (seam%nslave == 0) then
            call apply_tensor_operator(op, x, y, eliminate=mask)
            return
        end if
        allocate (xm(seam%nreduced), xf(seam%nfull), yf(seam%nfull))
        xm = x
        if (mask) then
            where (seam%boundary) xm = 0
        end if
        call seam_expand(seam, xm, xf)
        call apply_tensor_operator(op, xf, yf, eliminate=.false.)
        call seam_fold(seam, yf, y)
        if (mask) then
            where (seam%boundary) y = x
        end if
    end subroutine

    subroutine seam_diagonal(op, seam, diagonal)
        !! Jacobi diagonal of P^T A P without slave-slave cross terms.
        type(tensor_operator_t), intent(in) :: op
        type(tensor_seam_t), intent(in) :: seam
        real(dp), allocatable, intent(out) :: diagonal(:)
        real(dp), allocatable :: full(:)
        integer :: k, i
        call tensor_diagonal(op, full)
        if (seam%nslave == 0) then
            call move_alloc(full, diagonal)
            return
        end if
        allocate (diagonal(seam%nreduced))
        !$omp parallel do default(shared) private(k,i) schedule(static)
        do k = 1, seam%nreduced
            diagonal(k) = full(seam%primary(k))
            do i = seam%tptr(k), seam%tptr(k + 1) - 1
                diagonal(k) = diagonal(k) + seam%tval(i)**2* &
                              full(seam%slave(seam%trow(i)))
            end do
        end do
        !$omp end parallel do
        where (seam%boundary) diagonal = 1
    end subroutine

    subroutine assemble_seam_operator(op, seam, matrix)
        !! Assembled P^T A P with exact symmetry, the oracle and coarse matrix.
        type(tensor_operator_t), intent(in) :: op
        type(tensor_seam_t), intent(in) :: seam
        type(csr_matrix_t), intent(out) :: matrix
        type(csr_matrix_t) :: full
        integer, allocatable :: aptr(:), acol(:), mark(:), list(:), bptr(:), bcol(:)
        real(dp), allocatable :: aval(:), acc(:), bval(:)
        integer :: i, j, k, r, t, nl, nz, n, cap, s, e
        real(dp) :: avg
        if (seam%nslave == 0) then
            call assemble_tensor_operator(op, matrix)
            return
        end if
        call assemble_tensor_operator(op, full, eliminate=.false.)
        n = seam%nreduced
        allocate (mark(n), list(n), acc(n))
        mark = 0
        acc = 0
        ! AP = A P, rows in full numbering.
        allocate (aptr(seam%nfull + 1))
        cap = size(full%val)
        allocate (acol(cap), aval(cap))
        aptr(1) = 1
        nz = 0
        do i = 1, seam%nfull
            nl = 0
            do k = full%rowptr(i), full%rowptr(i + 1) - 1
                j = full%colidx(k)
                if (seam%reduced(j) > 0) then
                    call push(seam%reduced(j), full%val(k))
                else
                    r = -seam%reduced(j)
                    do t = seam%ptr(r), seam%ptr(r + 1) - 1
                        call push(seam%col(t), full%val(k)*seam%val(t))
                    end do
                end if
            end do
            call flush_row(nz, acol, aval)
            aptr(i + 1) = nz + 1
        end do
        ! B = P^T (AP), rows in reduced numbering.
        allocate (bptr(n + 1), bcol(cap), bval(cap))
        bptr(1) = 1
        nz = 0
        do j = 1, n
            nl = 0
            if (.not. seam%boundary(j)) then
                call add_row(seam%primary(j), 1.0_dp)
                do t = seam%tptr(j), seam%tptr(j + 1) - 1
                    call add_row(seam%slave(seam%trow(t)), seam%tval(t))
                end do
                do k = 1, nl
                    if (seam%boundary(list(k))) acc(list(k)) = 0
                end do
            else
                call push(j, 0.0_dp)
                acc(j) = 1
            end if
            call flush_row(nz, bcol, bval)
            bptr(j + 1) = nz + 1
        end do
        ! Exact symmetry for the symmetric coarse factorisation.
        do j = 1, n
            do k = bptr(j), bptr(j + 1) - 1
                i = bcol(k)
                if (i <= j) cycle
                s = bptr(i)
                e = bptr(i + 1) - 1
                do while (s <= e)
                    t = (s + e)/2
                    if (bcol(t) == j) exit
                    if (bcol(t) < j) then
                        s = t + 1
                    else
                        e = t - 1
                    end if
                end do
                if (s > e) then
                    avg = 0
                else
                    avg = (bval(k) + bval(t))/2
                    bval(t) = avg
                end if
                bval(k) = avg
            end do
        end do
        matrix%n = n
        matrix%rowptr = bptr
        matrix%colidx = bcol(1:nz)
        matrix%val = bval(1:nz)
    contains
        subroutine push(c, v)
            integer, intent(in) :: c
            real(dp), intent(in) :: v
            if (mark(c) == 0) then
                nl = nl + 1
                list(nl) = c
                mark(c) = 1
                acc(c) = 0
            end if
            acc(c) = acc(c) + v
        end subroutine

        subroutine add_row(row, weight)
            integer, intent(in) :: row
            real(dp), intent(in) :: weight
            integer :: q
            do q = aptr(row), aptr(row + 1) - 1
                call push(acol(q), weight*aval(q))
            end do
        end subroutine

        subroutine flush_row(nz, col, val)
            integer, intent(inout) :: nz
            integer, allocatable, intent(inout) :: col(:)
            real(dp), allocatable, intent(inout) :: val(:)
            integer, allocatable :: ic(:)
            real(dp), allocatable :: rv(:)
            integer :: q, a, b, tmp
            ! Insertion sort keeps each row canonical.
            do a = 2, nl
                b = a
                do while (b > 1)
                    if (list(b) >= list(b - 1)) exit
                    tmp = list(b); list(b) = list(b - 1); list(b - 1) = tmp
                    b = b - 1
                end do
            end do
            if (nz + nl > size(col)) then
                allocate (ic(2*(nz + nl)), rv(2*(nz + nl)))
                ic(1:nz) = col(1:nz)
                rv(1:nz) = val(1:nz)
                call move_alloc(ic, col)
                call move_alloc(rv, val)
            end if
            do q = 1, nl
                nz = nz + 1
                col(nz) = list(q)
                val(nz) = acc(list(q))
                mark(list(q)) = 0
                acc(list(q)) = 0
            end do
        end subroutine
    end subroutine
end module
