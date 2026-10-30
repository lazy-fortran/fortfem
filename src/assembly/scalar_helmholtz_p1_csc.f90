module fortfem_scalar_helmholtz_p1_csc
    !! P1 scalar Helmholtz cells sum_d g_d (d_d v, d_d u) - k**2 m (v, u) on
    !! intervals, triangles and tetrahedra, added into a CSC matrix through
    !! fortfem_block_csc, with JVP and VJP.
    !!
    !! The diagonal complex coefficients g(dim, cell) and m(cell) come from a
    !! complex-coordinate stretch (fortfem_cartesian_helmholtz_pml); g = m = 1 is
    !! the plain Helmholtz operator. The wave number k is real. Element matrices are
    !! the generated kernels of gen_scalar_helmholtz_products. Vertex v is dof v.
    !! The pattern may hold further groups; these routines only touch the given
    !! cell group. Memory is O(cells + nnz).
    !!
    !! A cell is rejected (status 1) when its length is not positive (intervals),
    !! its area is at most tiny (triangles, either orientation) or its signed
    !! volume is at most tiny (tetrahedra, positive orientation only).
    !!
    !! Cotangents follow sparse_direct_solve_constrained_vjp: a complex input z
    !! gets z_bar = sum conj(dA/dz) A_bar and a real input x gets
    !! Re(sum conj(dA/dx) A_bar), so Re<A_bar, A_dot> = <x_bar, x_dot> +
    !! Re<z_bar, z_dot>.
    use fortfem_block_csc, only: block_csc_add, block_csc_gather, &
        block_csc_pattern_t
    use fortfem_element_to_csr, only: build_node_to_element
    use fortfem_generated_p1_helmholtz_interval, only: &
        generated_p1_helmholtz_interval
    use fortfem_generated_p1_helmholtz_interval_jvp, only: &
        generated_p1_helmholtz_interval_jvp
    use fortfem_generated_p1_helmholtz_interval_vjp, only: &
        generated_p1_helmholtz_interval_vjp
    use fortfem_generated_p1_helmholtz_tetra, only: generated_p1_helmholtz_tetra
    use fortfem_generated_p1_helmholtz_tetra_jvp, only: &
        generated_p1_helmholtz_tetra_jvp
    use fortfem_generated_p1_helmholtz_tetra_vjp, only: &
        generated_p1_helmholtz_tetra_vjp
    use fortfem_generated_p1_helmholtz_triangle, only: &
        generated_p1_helmholtz_triangle
    use fortfem_generated_p1_helmholtz_triangle_jvp, only: &
        generated_p1_helmholtz_triangle_jvp
    use fortfem_generated_p1_helmholtz_triangle_vjp, only: &
        generated_p1_helmholtz_triangle_vjp
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: helmholtz_p1_assemble
    public :: helmholtz_p1_assemble_jvp
    public :: helmholtz_p1_assemble_vjp

contains

    subroutine helmholtz_p1_assemble(pattern, group, vertices, cells, gradient, &
            mass, wave_number, values, status)
        !! values += sum over cells of the element matrices. vertices(dim, nv) and
        !! cells(dim + 1, nc) with dim = 1, 2 or 3; gradient(dim, nc), mass(nc).
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, cells(:, :)
        real(dp), intent(in) :: vertices(:, :), wave_number
        complex(dp), intent(in) :: gradient(:, :), mass(:)
        complex(dp), intent(inout) :: values(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: blocks(:, :, :)
        real(dp), allocatable :: geometry(:, :)
        integer :: e, nloc

        call prepare(pattern, vertices, cells, gradient, mass, geometry, status)
        if (status /= 0) return
        nloc = size(cells, 1)
        allocate(blocks(nloc, nloc, size(cells, 2)))
        !$omp parallel do default(none) &
        !$omp shared(nloc, geometry, gradient, mass, wave_number, blocks) &
        !$omp private(e) schedule(static)
        do e = 1, size(geometry, 2)
            select case (nloc)
            case (2)
                call generated_p1_helmholtz_interval(geometry(1, e), &
                    gradient(1, e), mass(e), wave_number, blocks(:, :, e))
            case (3)
                call generated_p1_helmholtz_triangle(geometry(1, e), &
                    geometry(2, e), geometry(3, e), geometry(4, e), &
                    geometry(5, e), gradient(1, e), gradient(2, e), mass(e), &
                    wave_number, blocks(:, :, e))
            case default
                call generated_p1_helmholtz_tetra(geometry(1, e), &
                    geometry(2, e), geometry(3, e), geometry(4, e), &
                    geometry(5, e), geometry(6, e), geometry(7, e), &
                    geometry(8, e), geometry(9, e), geometry(10, e), &
                    gradient(1, e), gradient(2, e), gradient(3, e), mass(e), &
                    wave_number, blocks(:, :, e))
            end select
        end do
        !$omp end parallel do
        call block_csc_add(pattern, group, blocks, values)
    end subroutine helmholtz_p1_assemble

    subroutine helmholtz_p1_assemble_jvp(pattern, group, vertices, cells, &
            gradient, mass, wave_number, vertices_dot, gradient_dot, mass_dot, &
            wave_number_dot, values_dot, status)
        !! values_dot += directional derivative of the assembled cells; assembly
        !! is linear in the blocks, so this assembles the cell JVPs.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, cells(:, :)
        real(dp), intent(in) :: vertices(:, :), wave_number
        complex(dp), intent(in) :: gradient(:, :), mass(:)
        real(dp), intent(in) :: vertices_dot(:, :), wave_number_dot
        complex(dp), intent(in) :: gradient_dot(:, :), mass_dot(:)
        complex(dp), intent(inout) :: values_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: blocks(:, :, :)
        real(dp), allocatable :: geometry(:, :)
        real(dp) :: tangent(9)
        integer :: e, nloc

        call prepare(pattern, vertices, cells, gradient, mass, geometry, status)
        if (status /= 0) return
        status = 1
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (any(shape(gradient_dot) /= shape(gradient))) return
        if (size(mass_dot) /= size(mass)) return
        status = 0
        nloc = size(cells, 1)
        allocate(blocks(nloc, nloc, size(cells, 2)))
        !$omp parallel do default(none) &
        !$omp shared(nloc, geometry, cells, vertices_dot, gradient, mass) &
        !$omp shared(wave_number, gradient_dot, mass_dot, wave_number_dot, blocks) &
        !$omp private(e, tangent) schedule(static)
        do e = 1, size(geometry, 2)
            tangent = 0.0_dp
            call cell_entries(vertices_dot(:, cells(:, e)), tangent)
            select case (nloc)
            case (2)
                call generated_p1_helmholtz_interval_jvp(geometry(1, e), &
                    gradient(1, e), mass(e), wave_number, tangent(1), &
                    gradient_dot(1, e), mass_dot(e), wave_number_dot, &
                    blocks(:, :, e))
            case (3)
                call generated_p1_helmholtz_triangle_jvp(geometry(1, e), &
                    geometry(2, e), geometry(3, e), geometry(4, e), &
                    geometry(5, e), gradient(1, e), gradient(2, e), mass(e), &
                    wave_number, tangent(1), tangent(2), tangent(3), &
                    tangent(4), gradient_dot(1, e), gradient_dot(2, e), &
                    mass_dot(e), wave_number_dot, blocks(:, :, e))
            case default
                call generated_p1_helmholtz_tetra_jvp(geometry(1, e), &
                    geometry(2, e), geometry(3, e), geometry(4, e), &
                    geometry(5, e), geometry(6, e), geometry(7, e), &
                    geometry(8, e), geometry(9, e), geometry(10, e), &
                    gradient(1, e), gradient(2, e), gradient(3, e), mass(e), &
                    wave_number, tangent(1), tangent(2), tangent(3), &
                    tangent(4), tangent(5), tangent(6), tangent(7), &
                    tangent(8), tangent(9), gradient_dot(1, e), &
                    gradient_dot(2, e), gradient_dot(3, e), mass_dot(e), &
                    wave_number_dot, blocks(:, :, e))
            end select
        end do
        !$omp end parallel do
        call block_csc_add(pattern, group, blocks, values_dot)
    end subroutine helmholtz_p1_assemble_jvp

    subroutine helmholtz_p1_assemble_vjp(pattern, group, vertices, cells, &
            gradient, mass, wave_number, values_bar, vertices_bar, gradient_bar, &
            mass_bar, wave_number_bar, status)
        !! Cotangents of the cell contribution for the value cotangent values_bar
        !! (gather is the adjoint of add). The outputs are overwritten.
        type(block_csc_pattern_t), intent(in) :: pattern
        integer, intent(in) :: group, cells(:, :)
        real(dp), intent(in) :: vertices(:, :), wave_number
        complex(dp), intent(in) :: gradient(:, :), mass(:), values_bar(:)
        real(dp), intent(out) :: vertices_bar(:, :)
        complex(dp), intent(out) :: gradient_bar(:, :), mass_bar(:)
        real(dp), intent(out) :: wave_number_bar
        integer, intent(out) :: status

        complex(dp), allocatable :: blocks_bar(:, :, :)
        real(dp), allocatable :: geometry(:, :), local_bar(:, :, :), k_bar(:)
        complex(dp) :: bars(14)
        integer :: e, nloc, ngeo, dim

        vertices_bar = 0.0_dp
        gradient_bar = (0.0_dp, 0.0_dp)
        mass_bar = (0.0_dp, 0.0_dp)
        wave_number_bar = 0.0_dp
        call prepare(pattern, vertices, cells, gradient, mass, geometry, status)
        if (status /= 0) return
        status = 1
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (any(shape(gradient_bar) /= shape(gradient))) return
        if (size(mass_bar) /= size(mass)) return
        if (size(values_bar) /= pattern%nnz) return
        status = 0
        nloc = size(cells, 1)
        dim = nloc - 1
        ngeo = dim*dim
        allocate(blocks_bar(nloc, nloc, size(cells, 2)))
        allocate(local_bar(dim, nloc, size(cells, 2)), k_bar(size(cells, 2)))
        call block_csc_gather(pattern, group, values_bar, blocks_bar)
        ! The generated VJP is the plain transpose: feed conj(A_bar), conjugate.
        !$omp parallel do default(none) &
        !$omp shared(nloc, dim, ngeo, geometry, gradient, mass, wave_number) &
        !$omp shared(blocks_bar, local_bar, gradient_bar, mass_bar, k_bar) &
        !$omp private(e, bars) schedule(static)
        do e = 1, size(geometry, 2)
            select case (nloc)
            case (2)
                call generated_p1_helmholtz_interval_vjp(geometry(1, e), &
                    gradient(1, e), mass(e), wave_number, &
                    conjg(blocks_bar(:, :, e)), bars(1), bars(2), bars(3), &
                    bars(4))
            case (3)
                call generated_p1_helmholtz_triangle_vjp(geometry(1, e), &
                    geometry(2, e), geometry(3, e), geometry(4, e), &
                    geometry(5, e), gradient(1, e), gradient(2, e), mass(e), &
                    wave_number, conjg(blocks_bar(:, :, e)), bars(1), bars(2), &
                    bars(3), bars(4), bars(5), bars(6), bars(7), bars(8))
            case default
                call generated_p1_helmholtz_tetra_vjp(geometry(1, e), &
                    geometry(2, e), geometry(3, e), geometry(4, e), &
                    geometry(5, e), geometry(6, e), geometry(7, e), &
                    geometry(8, e), geometry(9, e), geometry(10, e), &
                    gradient(1, e), gradient(2, e), gradient(3, e), mass(e), &
                    wave_number, conjg(blocks_bar(:, :, e)), bars(1), bars(2), &
                    bars(3), bars(4), bars(5), bars(6), bars(7), bars(8), &
                    bars(9), bars(10), bars(11), bars(12), bars(13), bars(14))
            end select
            call cell_entries_transpose(real(bars(1:ngeo), dp), &
                local_bar(:, :, e))
            gradient_bar(:, e) = conjg(bars(ngeo + 1:ngeo + dim))
            mass_bar(e) = conjg(bars(ngeo + dim + 1))
            k_bar(e) = real(bars(ngeo + dim + 2), dp)
        end do
        !$omp end parallel do
        ! Ascending cell order: independent of the thread count.
        do e = 1, size(k_bar)
            wave_number_bar = wave_number_bar + k_bar(e)
        end do
        call scatter_to_vertices(cells, local_bar, vertices_bar)
    end subroutine helmholtz_p1_assemble_vjp

    subroutine prepare(pattern, vertices, cells, gradient, mass, geometry, status)
        !! Cell geometry: the entries of x_(k+1) - x_1 and, for triangles and
        !! tetrahedra, the orientation sign(det) as the last row.
        type(block_csc_pattern_t), intent(in) :: pattern
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        complex(dp), intent(in) :: gradient(:, :), mass(:)
        real(dp), allocatable, intent(out) :: geometry(:, :)
        integer, intent(out) :: status

        real(dp) :: det
        integer :: dim, e, nc, ngeo
        logical :: degenerate

        status = 1
        dim = size(vertices, 1)
        nc = size(cells, 2)
        if (dim < 1 .or. dim > 3) return
        if (size(cells, 1) /= dim + 1) return
        if (size(vertices, 2) /= pattern%n) return
        if (any(cells < 1) .or. any(cells > size(vertices, 2))) return
        if (size(gradient, 1) /= dim .or. size(gradient, 2) /= nc) return
        if (size(mass) /= nc) return
        ngeo = dim*dim
        allocate(geometry(merge(1, ngeo + 1, dim == 1), nc))
        degenerate = .false.
        !$omp parallel do default(none) &
        !$omp shared(geometry, vertices, cells, nc, dim, ngeo) &
        !$omp private(e, det) reduction(.or.: degenerate) schedule(static)
        do e = 1, nc
            call cell_entries(vertices(:, cells(:, e)), geometry(1:ngeo, e))
            select case (dim)
            case (1)
                if (geometry(1, e) <= 0.0_dp) degenerate = .true.
            case (2)
                det = geometry(1, e)*geometry(4, e) - geometry(2, e)*geometry(3, e)
                if (0.5_dp*abs(det) <= tiny(1.0_dp)) degenerate = .true.
                geometry(5, e) = sign(1.0_dp, det)
            case default
                det = determinant3(geometry(1:9, e))
                if (det <= tiny(1.0_dp)) degenerate = .true.
                geometry(10, e) = 1.0_dp
            end select
        end do
        !$omp end parallel do
        if (degenerate) return
        status = 0
    end subroutine prepare

    pure subroutine cell_entries(corners, entries)
        !! entries = the columns x_(k+1) - x_1, k = 1..dim, of corners(dim, dim+1).
        real(dp), intent(in) :: corners(:, :)
        real(dp), intent(inout) :: entries(:)
        integer :: dim, k

        dim = size(corners, 1)
        do k = 1, dim
            entries(dim*(k - 1) + 1:dim*k) = corners(:, k + 1) - corners(:, 1)
        end do
    end subroutine cell_entries

    pure subroutine cell_entries_transpose(entries_bar, corners_bar)
        !! Adjoint of cell_entries.
        real(dp), intent(in) :: entries_bar(:)
        real(dp), intent(out) :: corners_bar(:, :)
        integer :: dim, k

        dim = size(corners_bar, 1)
        corners_bar(:, 1) = 0.0_dp
        do k = 1, dim
            corners_bar(:, k + 1) = entries_bar(dim*(k - 1) + 1:dim*k)
            corners_bar(:, 1) = corners_bar(:, 1) - corners_bar(:, k + 1)
        end do
    end subroutine cell_entries_transpose

    pure real(dp) function determinant3(columns) result(det)
        real(dp), intent(in) :: columns(9)

        det = columns(1)*(columns(5)*columns(9) - columns(6)*columns(8)) - &
            columns(2)*(columns(4)*columns(9) - columns(6)*columns(7)) + &
            columns(3)*(columns(4)*columns(8) - columns(5)*columns(7))
    end function determinant3

    subroutine scatter_to_vertices(cells, local_bar, vertices_bar)
        !! vertices_bar(:, v) = sum of local_bar over the cell corners at v, in
        !! ascending cell order (vertex-parallel, thread-count invariant).
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: local_bar(:, :, :)
        real(dp), intent(inout) :: vertices_bar(:, :)

        integer, allocatable :: element_ptr(:), element_idx(:)
        integer :: e, k, local, vertex

        call build_node_to_element(size(vertices_bar, 2), cells, element_ptr, &
            element_idx)
        !$omp parallel do default(none) &
        !$omp shared(vertices_bar, element_ptr, element_idx, cells, local_bar) &
        !$omp private(vertex, k, e, local) schedule(static)
        do vertex = 1, size(vertices_bar, 2)
            do k = element_ptr(vertex), element_ptr(vertex + 1) - 1
                e = element_idx(k)
                do local = 1, size(cells, 1)
                    if (cells(local, e) /= vertex) cycle
                    vertices_bar(:, vertex) = vertices_bar(:, vertex) + &
                        local_bar(:, local, e)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine scatter_to_vertices

end module fortfem_scalar_helmholtz_p1_csc
