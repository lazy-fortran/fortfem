module fortfem_elasticity_p1_csc_2d
    !! Vector P1 time-harmonic elasticity K - inertia*M on a triangle mesh, added
    !! into a CSC matrix through fortfem_block_csc, with JVP and VJP.
    !!
    !! Vertex v carries the dofs (2v - 1, 2v) = (u_x, u_y). Materials are given per
    !! triangle: complex Lame parameters lambda, mu and inertia = rho*omega**2.
    !! The pattern may hold further groups (boundary edges, interface blocks); the
    !! routines here only touch the cell group. Memory is O(triangles + nnz).
    !!
    !! Cotangents follow sparse_direct_solve_constrained_vjp: for a complex input z,
    !! z_bar = sum conj(dA/dz)*A_bar, and a real input x gets Re(sum conj(dA/dx)*
    !! A_bar), so Re<A_bar, A_dot> = <x_bar, x_dot> + Re<z_bar, z_dot>.
    use fortfem_block_csc, only: block_csc_add, block_csc_add_group, &
        block_csc_finalize, block_csc_gather, block_csc_init, block_csc_pattern_t
    use fortfem_element_to_csr, only: build_node_to_element
    use fortfem_generated_p1_elasticity_cell, only: generated_p1_elasticity_cell
    use fortfem_generated_p1_elasticity_cell_jvp, only: &
        generated_p1_elasticity_cell_jvp
    use fortfem_generated_p1_elasticity_cell_vjp, only: &
        generated_p1_elasticity_cell_vjp
    use fortfem_generated_triangle_measure_fp64, only: &
        generated_triangle_measure_fp64
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: elasticity_p1_cell_dofs
    public :: elasticity_p1_pattern
    public :: elasticity_p1_assemble
    public :: elasticity_p1_assemble_jvp
    public :: elasticity_p1_assemble_vjp

contains

    subroutine elasticity_p1_cell_dofs(triangles, dofs)
        !! Connectivity of the cell group: dofs(2k - 1:2k, e) of vertex k of e.
        integer, intent(in) :: triangles(:, :)
        integer, allocatable, intent(out) :: dofs(:, :)

        allocate(dofs(6, size(triangles, 2)))
        dofs(1::2, :) = 2*triangles - 1
        dofs(2::2, :) = 2*triangles
    end subroutine elasticity_p1_cell_dofs

    subroutine elasticity_p1_pattern(vertex_count, triangles, pattern, cell_group)
        !! Pattern of the elasticity operator alone. Callers coupling further blocks
        !! add elasticity_p1_cell_dofs as one group of their own pattern instead.
        integer, intent(in) :: vertex_count, triangles(:, :)
        type(block_csc_pattern_t), intent(out) :: pattern
        integer, intent(out) :: cell_group

        integer, allocatable :: dofs(:, :)

        call elasticity_p1_cell_dofs(triangles, dofs)
        call block_csc_init(pattern, 2*vertex_count)
        call block_csc_add_group(pattern, dofs, cell_group)
        call block_csc_finalize(pattern)
    end subroutine elasticity_p1_pattern

    subroutine elasticity_p1_assemble(pattern, cell_group, vertices, triangles, &
            lambda, mu, inertia, values, status)
        !! values += K - inertia*M (values is the CSC value array of the pattern).
        !! status = 1 for inconsistent shapes or a degenerate triangle.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: cell_group, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :)
        complex(dp), intent(in) :: lambda(:), mu(:), inertia(:)
        complex(dp), intent(inout) :: values(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: blocks(:, :, :)
        real(dp), allocatable :: geometry(:, :)
        integer :: e

        call prepare(pattern, vertices, triangles, lambda, mu, inertia, geometry, &
            status)
        if (status /= 0) return
        allocate(blocks(6, 6, size(triangles, 2)))
        !$omp parallel do default(none) &
        !$omp shared(geometry, lambda, mu, inertia, blocks) private(e) &
        !$omp schedule(static)
        do e = 1, size(geometry, 2)
            call generated_p1_elasticity_cell(geometry(1, e), geometry(2, e), &
                geometry(3, e), geometry(4, e), geometry(5, e), lambda(e), mu(e), &
                inertia(e), blocks(:, :, e))
        end do
        !$omp end parallel do
        call block_csc_add(pattern, cell_group, blocks, values)
    end subroutine elasticity_p1_assemble

    subroutine elasticity_p1_assemble_jvp(pattern, cell_group, vertices, &
            triangles, lambda, mu, inertia, vertices_dot, lambda_dot, mu_dot, &
            inertia_dot, values_dot, status)
        !! values_dot += directional derivative of K - inertia*M. Assembly is linear
        !! in the blocks, so this is the assembly of the cell JVPs.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: cell_group, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        complex(dp), intent(in) :: lambda(:), mu(:), inertia(:)
        complex(dp), intent(in) :: lambda_dot(:), mu_dot(:), inertia_dot(:)
        complex(dp), intent(inout) :: values_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: blocks(:, :, :)
        real(dp), allocatable :: geometry(:, :)
        real(dp) :: tangent(4)
        integer :: e

        call prepare(pattern, vertices, triangles, lambda, mu, inertia, geometry, &
            status)
        if (status /= 0) return
        status = 1
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (size(lambda_dot) /= size(lambda) .or. size(mu_dot) /= size(mu) .or. &
            size(inertia_dot) /= size(inertia)) return
        status = 0
        allocate(blocks(6, 6, size(triangles, 2)))
        !$omp parallel do default(none) &
        !$omp shared(geometry, triangles, vertices_dot, lambda, mu, inertia) &
        !$omp shared(lambda_dot, mu_dot, inertia_dot, blocks) &
        !$omp private(e, tangent) schedule(static)
        do e = 1, size(geometry, 2)
            tangent = affine_entries(vertices_dot(:, triangles(:, e)))
            call generated_p1_elasticity_cell_jvp(geometry(1, e), geometry(2, e), &
                geometry(3, e), geometry(4, e), geometry(5, e), lambda(e), mu(e), &
                inertia(e), tangent(1), tangent(2), tangent(3), tangent(4), &
                lambda_dot(e), mu_dot(e), inertia_dot(e), blocks(:, :, e))
        end do
        !$omp end parallel do
        call block_csc_add(pattern, cell_group, blocks, values_dot)
    end subroutine elasticity_p1_assemble_jvp

    subroutine elasticity_p1_assemble_vjp(pattern, cell_group, vertices, &
            triangles, lambda, mu, inertia, values_bar, vertices_bar, lambda_bar, &
            mu_bar, inertia_bar, status)
        !! Cotangents of the cell contribution for the value cotangent values_bar:
        !! the cell VJPs of the gathered blocks (gather is the adjoint of add).
        !! The outputs are overwritten.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: cell_group, triangles(:, :)
        real(dp), intent(in) :: vertices(:, :)
        complex(dp), intent(in) :: lambda(:), mu(:), inertia(:), values_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        complex(dp), intent(out) :: lambda_bar(:), mu_bar(:), inertia_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: blocks_bar(:, :, :)
        real(dp), allocatable :: geometry(:, :), local_bar(:, :, :)
        integer, allocatable :: element_ptr(:), element_idx(:)
        complex(dp) :: a_bar, b_bar, c_bar, d_bar, l_bar, m_bar, i_bar
        integer :: e, k, local, vertex

        vertices_bar = 0.0_dp
        lambda_bar = (0.0_dp, 0.0_dp)
        mu_bar = (0.0_dp, 0.0_dp)
        inertia_bar = (0.0_dp, 0.0_dp)
        call prepare(pattern, vertices, triangles, lambda, mu, inertia, geometry, &
            status)
        if (status /= 0) return
        status = 1
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(lambda_bar) /= size(lambda) .or. size(mu_bar) /= size(mu) .or. &
            size(inertia_bar) /= size(inertia)) return
        if (size(values_bar) /= pattern%nnz) return
        status = 0
        allocate(blocks_bar(6, 6, size(triangles, 2)))
        allocate(local_bar(2, 3, size(triangles, 2)))
        call block_csc_gather(pattern, cell_group, values_bar, blocks_bar)
        ! The generated VJP is the plain transpose: feed conj(A_bar) and conjugate.
        !$omp parallel do default(none) &
        !$omp shared(geometry, lambda, mu, inertia, blocks_bar, local_bar) &
        !$omp shared(lambda_bar, mu_bar, inertia_bar) &
        !$omp private(e, a_bar, b_bar, c_bar, d_bar, l_bar, m_bar, i_bar) &
        !$omp schedule(static)
        do e = 1, size(geometry, 2)
            call generated_p1_elasticity_cell_vjp(geometry(1, e), geometry(2, e), &
                geometry(3, e), geometry(4, e), geometry(5, e), lambda(e), mu(e), &
                inertia(e), conjg(blocks_bar(:, :, e)), a_bar, b_bar, c_bar, &
                d_bar, l_bar, m_bar, i_bar)
            lambda_bar(e) = conjg(l_bar)
            mu_bar(e) = conjg(m_bar)
            inertia_bar(e) = conjg(i_bar)
            ! a = x2 - x1, b = y2 - y1, c = x3 - x1, d = y3 - y1
            local_bar(:, 2, e) = [real(a_bar, dp), real(b_bar, dp)]
            local_bar(:, 3, e) = [real(c_bar, dp), real(d_bar, dp)]
            local_bar(:, 1, e) = -local_bar(:, 2, e) - local_bar(:, 3, e)
        end do
        !$omp end parallel do
        ! Vertex-parallel sums in ascending triangle order: thread-count invariant.
        call build_node_to_element(size(vertices, 2), triangles, element_ptr, &
            element_idx)
        !$omp parallel do default(none) &
        !$omp shared(vertices_bar, element_ptr, element_idx, triangles, local_bar) &
        !$omp private(vertex, k, e, local) schedule(static)
        do vertex = 1, size(vertices, 2)
            do k = element_ptr(vertex), element_ptr(vertex + 1) - 1
                e = element_idx(k)
                do local = 1, 3
                    if (triangles(local, e) /= vertex) cycle
                    vertices_bar(:, vertex) = vertices_bar(:, vertex) + &
                        local_bar(:, local, e)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine elasticity_p1_assemble_vjp

    subroutine prepare(pattern, vertices, triangles, lambda, mu, inertia, &
            geometry, status)
        !! Affine entries (a, b, c, d) and orientation sign(det) per triangle.
        type(block_csc_pattern_t), intent(in) :: pattern
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: lambda(:), mu(:), inertia(:)
        real(dp), allocatable, intent(out) :: geometry(:, :)
        integer, intent(out) :: status

        real(dp) :: det, scale
        integer :: e, nt
        logical :: degenerate

        status = 1
        nt = size(triangles, 2)
        if (size(vertices, 1) /= 2 .or. size(triangles, 1) /= 3) return
        if (2*size(vertices, 2) /= pattern%n) return
        if (any(triangles < 1) .or. any(triangles > size(vertices, 2))) return
        if (size(lambda) /= nt .or. size(mu) /= nt .or. size(inertia) /= nt) return
        allocate(geometry(5, nt))
        degenerate = .false.
        !$omp parallel do default(none) shared(geometry, vertices, triangles, nt) &
        !$omp private(e, det, scale) reduction(.or.: degenerate) schedule(static)
        do e = 1, nt
            geometry(1:4, e) = affine_entries(vertices(:, triangles(:, e)))
            call generated_triangle_measure_fp64(geometry(1, e), geometry(2, e), &
                geometry(3, e), geometry(4, e), det, scale)
            if (abs(det) <= epsilon(1.0_dp)*scale) degenerate = .true.
            geometry(5, e) = sign(1.0_dp, det)
        end do
        !$omp end parallel do
        if (degenerate) return
        status = 0
    end subroutine prepare

    pure function affine_entries(corners) result(entries)
        !! (a, b, c, d) = (x2 - x1, y2 - y1, x3 - x1, y3 - y1) of corners(2, 3).
        real(dp), intent(in) :: corners(:, :)
        real(dp) :: entries(4)

        entries = [corners(:, 2) - corners(:, 1), corners(:, 3) - corners(:, 1)]
    end function affine_entries

end module fortfem_elasticity_p1_csc_2d
