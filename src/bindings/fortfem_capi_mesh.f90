module fortfem_capi_mesh
    use iso_c_binding, only: c_double, c_double_complex, c_int
    use fortfem_kinds, only: dp
    use fortfem_assembly_nedelec_2d, only: &
        assemble_nedelec_axisymmetric_fourier_csc
    use fortfem_assembly_mixed_2d, only: assemble_nedelec_rt_mass_csc
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_mesh_extension_2d, only: extend_triangle_mesh
    use fortfem_rt_field_2d, only: evaluate_rt_field_2d, &
        reconstruct_axisymmetric_fourier_toroidal, rt_l2_norm
    use fortsparse, only: csc_t, fortsparse_status_t
    implicit none
    private

    integer, parameter :: max_meshes = 32
    type(mesh_2d_t), save :: meshes(max_meshes)
    logical, save :: mesh_occupied(max_meshes) = .false.

    public :: fortfem_triangle_edge_map
    public :: fortfem_triangle_mesh_create
    public :: fortfem_triangle_mesh_extend
    public :: fortfem_triangle_mesh_edges
    public :: fortfem_triangle_mesh_free
    public :: fortfem_rt0_evaluate
    public :: fortfem_rt0_l2_norm
    public :: fortfem_rt0_toroidal
    public :: fortfem_nedelec_axisymmetric_fourier_csc
    public :: fortfem_nedelec_rt0_mass_csc

contains

    subroutine fortfem_triangle_mesh_create( &
            n_vertices, vertices, n_triangles, triangles, &
            handle, n_edges, status) &
            bind(c, name="fortfem_triangle_mesh_create")
        integer(c_int), value :: n_vertices, n_triangles
        real(c_double), intent(in) :: vertices(*)
        integer(c_int), intent(in) :: triangles(*)
        integer(c_int), intent(out) :: handle, n_edges, status

        integer :: slot

        handle = 0_c_int
        n_edges = 0_c_int
        status = -1_c_int
        if (n_vertices < 3_c_int .or. n_triangles < 1_c_int) return

        slot = find_free_mesh_slot()
        if (slot == 0) then
            status = -3_c_int
            return
        end if
        call initialize_mesh( &
            meshes(slot), n_vertices, vertices, n_triangles, triangles, status)
        if (status /= 0_c_int) return

        mesh_occupied(slot) = .true.
        handle = int(slot, c_int)
        n_edges = int(meshes(slot)%n_edges, c_int)
    end subroutine fortfem_triangle_mesh_create

    subroutine fortfem_triangle_mesh_edges( &
            handle, edge_capacity, n_edges, edges, &
            triangle_edge_dofs, triangle_edge_signs, status) &
            bind(c, name="fortfem_triangle_mesh_edges")
        integer(c_int), value :: handle, edge_capacity
        integer(c_int), intent(out) :: n_edges, edges(*)
        integer(c_int), intent(out) :: triangle_edge_dofs(*)
        integer(c_int), intent(out) :: triangle_edge_signs(*)
        integer(c_int), intent(out) :: status

        integer :: slot

        n_edges = 0_c_int
        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return
        call export_edge_map(meshes(slot), edge_capacity, n_edges, edges, &
            triangle_edge_dofs, triangle_edge_signs, status)
    end subroutine fortfem_triangle_mesh_edges

    subroutine fortfem_triangle_mesh_extend( &
            core_handle, n_outer_vertices, outer_vertices, &
            handle, n_vertices, n_triangles, n_edges, &
            first_boundary_dof, status) &
            bind(c, name="fortfem_triangle_mesh_extend")
        integer(c_int), value :: core_handle, n_outer_vertices
        real(c_double), intent(in) :: outer_vertices(*)
        integer(c_int), intent(out) :: handle, n_vertices, n_triangles
        integer(c_int), intent(out) :: n_edges, first_boundary_dof, status

        real(dp), allocatable :: fortran_outer_vertices(:, :)
        integer :: component, core_slot, extension_status, slot, vertex

        handle = 0_c_int
        n_vertices = 0_c_int
        n_triangles = 0_c_int
        n_edges = 0_c_int
        first_boundary_dof = 0_c_int
        status = -3_c_int
        core_slot = int(core_handle)
        if (.not. valid_mesh_handle(core_slot)) return

        status = -1_c_int
        if (n_outer_vertices < 3_c_int) return
        slot = find_free_mesh_slot()
        if (slot == 0) then
            status = -3_c_int
            return
        end if
        allocate(fortran_outer_vertices(2, int(n_outer_vertices)))
        do vertex = 1, int(n_outer_vertices)
            do component = 1, 2
                fortran_outer_vertices(component, vertex) = real( &
                    outer_vertices(2 * (vertex - 1) + component), dp)
            end do
        end do

        call extend_triangle_mesh(meshes(core_slot), fortran_outer_vertices, &
            meshes(slot), extension_status)
        if (extension_status /= 0) then
            call meshes(slot)%destroy()
            status = int(extension_status, c_int)
            return
        end if

        mesh_occupied(slot) = .true.
        handle = int(slot, c_int)
        n_vertices = int(meshes(slot)%n_vertices, c_int)
        n_triangles = int(meshes(slot)%n_triangles, c_int)
        n_edges = int(meshes(slot)%n_edges, c_int)
        first_boundary_dof = int(meshes(slot)%n_interior_dofs, c_int)
        status = 0_c_int
    end subroutine fortfem_triangle_mesh_extend

    subroutine fortfem_triangle_mesh_free(handle, status) &
            bind(c, name="fortfem_triangle_mesh_free")
        integer(c_int), value :: handle
        integer(c_int), intent(out) :: status

        integer :: slot

        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return

        call meshes(slot)%destroy()
        mesh_occupied(slot) = .false.
        status = 0_c_int
    end subroutine fortfem_triangle_mesh_free

    subroutine fortfem_rt0_evaluate( &
            handle, triangle, x, y, n_dofs, dofs, &
            value, divergence, status) bind(c, name="fortfem_rt0_evaluate")
        integer(c_int), value :: handle, triangle, n_dofs
        real(c_double), value :: x, y
        complex(c_double_complex), intent(in) :: dofs(*)
        complex(c_double_complex), intent(out) :: value(*), divergence
        integer(c_int), intent(out) :: status

        complex(dp), allocatable :: fortran_dofs(:)
        complex(dp) :: fortran_value(2), fortran_divergence
        integer :: component, slot

        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return
        status = -1_c_int
        if (n_dofs /= int(meshes(slot)%n_edges, c_int)) return
        if (triangle < 0_c_int .or. &
            triangle >= int(meshes(slot)%n_triangles, c_int)) return

        call copy_complex_dofs(n_dofs, dofs, fortran_dofs)
        call evaluate_rt_field_2d(meshes(slot), int(triangle) + 1, &
            real(x, dp), real(y, dp), fortran_dofs, &
            fortran_value, fortran_divergence)
        do component = 1, 2
            value(component) = to_c_complex(fortran_value(component))
        end do
        divergence = to_c_complex(fortran_divergence)
        status = 0_c_int
    end subroutine fortfem_rt0_evaluate

    subroutine fortfem_rt0_l2_norm( &
            handle, n_dofs, dofs, norm, status) &
            bind(c, name="fortfem_rt0_l2_norm")
        integer(c_int), value :: handle, n_dofs
        complex(c_double_complex), intent(in) :: dofs(*)
        real(c_double), intent(out) :: norm
        integer(c_int), intent(out) :: status

        complex(dp), allocatable :: fortran_dofs(:)
        integer :: slot

        norm = 0.0_c_double
        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return
        status = -1_c_int
        if (n_dofs /= int(meshes(slot)%n_edges, c_int)) return

        call copy_complex_dofs(n_dofs, dofs, fortran_dofs)
        norm = real(rt_l2_norm(meshes(slot), fortran_dofs), c_double)
        status = 0_c_int
    end subroutine fortfem_rt0_l2_norm

    subroutine fortfem_rt0_toroidal( &
            handle, mode, n_dofs, dofs, toroidal, status) &
            bind(c, name="fortfem_rt0_toroidal")
        integer(c_int), value :: handle, mode, n_dofs
        complex(c_double_complex), intent(in) :: dofs(*)
        complex(c_double_complex), intent(out) :: toroidal(*)
        integer(c_int), intent(out) :: status

        complex(dp), allocatable :: fortran_dofs(:), fortran_toroidal(:)
        integer :: slot, triangle

        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return
        status = -1_c_int
        if (n_dofs /= int(meshes(slot)%n_edges, c_int)) return
        if (mode == 0_c_int) return

        call copy_complex_dofs(n_dofs, dofs, fortran_dofs)
        allocate(fortran_toroidal(meshes(slot)%n_triangles))
        call reconstruct_axisymmetric_fourier_toroidal( &
            meshes(slot), int(mode), fortran_dofs, fortran_toroidal)
        do triangle = 1, meshes(slot)%n_triangles
            toroidal(triangle) = to_c_complex(fortran_toroidal(triangle))
        end do
        status = 0_c_int
    end subroutine fortfem_rt0_toroidal

    subroutine fortfem_nedelec_axisymmetric_fourier_csc( &
            handle, mode, quadrature_degree, nnz_capacity, &
            n_dofs, nnz, col_ptr, row_ind, values, status) &
            bind(c, name="fortfem_nedelec_axisymmetric_fourier_csc")
        integer(c_int), value :: handle, mode, quadrature_degree
        integer(c_int), value :: nnz_capacity
        integer(c_int), intent(out) :: n_dofs, nnz
        integer(c_int), intent(out) :: col_ptr(*), row_ind(*)
        real(c_double), intent(out) :: values(*)
        integer(c_int), intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        integer :: slot

        n_dofs = 0_c_int
        nnz = 0_c_int
        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return

        n_dofs = int(meshes(slot)%n_edges, c_int)
        status = -1_c_int
        if (quadrature_degree < 1_c_int) return
        if (.not. valid_axisymmetric_mesh(meshes(slot))) return

        call assemble_nedelec_axisymmetric_fourier_csc( &
            meshes(slot), int(mode), int(quadrature_degree), &
            matrix, sparse_status)
        if (sparse_status%code /= 0) then
            status = int(sparse_status%code, c_int)
            return
        end if

        call export_real_csc(matrix, nnz_capacity, n_dofs, nnz, &
            col_ptr, row_ind, values, status)
    end subroutine fortfem_nedelec_axisymmetric_fourier_csc

    subroutine fortfem_nedelec_rt0_mass_csc( &
            handle, quadrature_degree, nnz_capacity, &
            n_dofs, nnz, col_ptr, row_ind, values, status) &
            bind(c, name="fortfem_nedelec_rt0_mass_csc")
        integer(c_int), value :: handle, quadrature_degree, nnz_capacity
        integer(c_int), intent(out) :: n_dofs, nnz
        integer(c_int), intent(out) :: col_ptr(*), row_ind(*)
        real(c_double), intent(out) :: values(*)
        integer(c_int), intent(out) :: status

        type(csc_t) :: matrix
        type(fortsparse_status_t) :: sparse_status
        integer :: slot

        n_dofs = 0_c_int
        nnz = 0_c_int
        status = -3_c_int
        slot = int(handle)
        if (.not. valid_mesh_handle(slot)) return

        n_dofs = int(meshes(slot)%n_edges, c_int)
        status = -1_c_int
        if (quadrature_degree < 1_c_int) return
        if (.not. counterclockwise_mesh(meshes(slot))) return

        call assemble_nedelec_rt_mass_csc( &
            meshes(slot), int(quadrature_degree), matrix, sparse_status)
        if (sparse_status%code /= 0) then
            status = int(sparse_status%code, c_int)
            return
        end if
        call export_real_csc(matrix, nnz_capacity, n_dofs, nnz, &
            col_ptr, row_ind, values, status)
    end subroutine fortfem_nedelec_rt0_mass_csc

    subroutine fortfem_triangle_edge_map( &
            n_vertices, vertices, n_triangles, triangles, &
            edge_capacity, n_edges, edges, triangle_edge_dofs, &
            triangle_edge_signs, status) &
            bind(c, name="fortfem_triangle_edge_map")
        integer(c_int), value :: n_vertices, n_triangles, edge_capacity
        real(c_double), intent(in) :: vertices(*)
        integer(c_int), intent(in) :: triangles(*)
        integer(c_int), intent(out) :: n_edges, edges(*)
        integer(c_int), intent(out) :: triangle_edge_dofs(*)
        integer(c_int), intent(out) :: triangle_edge_signs(*)
        integer(c_int), intent(out) :: status

        integer(c_int) :: handle, free_status

        if (edge_capacity < 1_c_int) then
            n_edges = 0_c_int
            status = -1_c_int
            return
        end if
        call fortfem_triangle_mesh_create( &
            n_vertices, vertices, n_triangles, triangles, &
            handle, n_edges, status)
        if (status /= 0_c_int) return
        call fortfem_triangle_mesh_edges(handle, edge_capacity, n_edges, &
            edges, triangle_edge_dofs, triangle_edge_signs, status)
        call fortfem_triangle_mesh_free(handle, free_status)
    end subroutine fortfem_triangle_edge_map

    subroutine initialize_mesh( &
            mesh, n_vertices, vertices, n_triangles, triangles, status)
        type(mesh_2d_t), intent(inout) :: mesh
        integer(c_int), intent(in) :: n_vertices, n_triangles
        real(c_double), intent(in) :: vertices(*)
        integer(c_int), intent(in) :: triangles(*)
        integer(c_int), intent(out) :: status

        integer :: triangle_id, vertex_id, component, local_vertex

        status = -1_c_int
        mesh%n_vertices = int(n_vertices)
        mesh%n_triangles = int(n_triangles)
        mesh%has_triangles = .true.
        allocate(mesh%vertices(2, mesh%n_vertices))
        allocate(mesh%triangles(3, mesh%n_triangles))

        do vertex_id = 1, mesh%n_vertices
            do component = 1, 2
                mesh%vertices(component, vertex_id) = real( &
                    vertices(2 * (vertex_id - 1) + component), dp)
            end do
        end do

        do triangle_id = 1, mesh%n_triangles
            do local_vertex = 1, 3
                vertex_id = int( &
                    triangles(3 * (triangle_id - 1) + local_vertex))
                if (vertex_id < 0 .or. vertex_id >= mesh%n_vertices) then
                    call mesh%destroy()
                    return
                end if
                mesh%triangles(local_vertex, triangle_id) = vertex_id + 1
            end do
        end do

        call mesh%build_edge_connectivity()
        call mesh%build_edge_dof_numbering()
        status = 0_c_int
    end subroutine initialize_mesh

    subroutine export_edge_map( &
            mesh, edge_capacity, n_edges, edges, &
            triangle_edge_dofs, triangle_edge_signs, status)
        type(mesh_2d_t), intent(in) :: mesh
        integer(c_int), intent(in) :: edge_capacity
        integer(c_int), intent(out) :: n_edges, edges(*)
        integer(c_int), intent(out) :: triangle_edge_dofs(*)
        integer(c_int), intent(out) :: triangle_edge_signs(*)
        integer(c_int), intent(out) :: status

        integer :: triangle_id, local_edge, edge_id, dof_id
        integer :: edge_dofs(3), edge_signs(3)

        n_edges = int(mesh%n_edges, c_int)
        if (mesh%n_edges > int(edge_capacity)) then
            status = -2_c_int
            return
        end if

        do dof_id = 0, mesh%n_edges - 1
            edge_id = mesh%dof_to_edge(dof_id + 1)
            edges(2 * dof_id + 1) = int( &
                mesh%edges(1, edge_id) - 1, c_int)
            edges(2 * dof_id + 2) = int( &
                mesh%edges(2, edge_id) - 1, c_int)
        end do

        do triangle_id = 1, mesh%n_triangles
            call mesh%get_triangle_edge_dofs( &
                triangle_id, edge_dofs, edge_signs)
            do local_edge = 1, 3
                triangle_edge_dofs(3 * (triangle_id - 1) + local_edge) = &
                    int(edge_dofs(local_edge), c_int)
                triangle_edge_signs(3 * (triangle_id - 1) + local_edge) = &
                    int(edge_signs(local_edge), c_int)
            end do
        end do
        status = 0_c_int
    end subroutine export_edge_map

    integer function find_free_mesh_slot() result(slot)
        integer :: candidate

        slot = 0
        do candidate = 1, max_meshes
            if (.not. mesh_occupied(candidate)) then
                slot = candidate
                return
            end if
        end do
    end function find_free_mesh_slot

    logical function valid_mesh_handle(slot) result(valid)
        integer, intent(in) :: slot

        valid = slot >= 1 .and. slot <= max_meshes
        if (valid) valid = mesh_occupied(slot)
    end function valid_mesh_handle

    logical function valid_axisymmetric_mesh(mesh) result(valid)
        type(mesh_2d_t), intent(in) :: mesh

        valid = all(mesh%vertices(1, :) > 0.0_dp) .and. &
            counterclockwise_mesh(mesh)
    end function valid_axisymmetric_mesh

    logical function counterclockwise_mesh(mesh) result(valid)
        type(mesh_2d_t), intent(in) :: mesh

        real(dp) :: determinant
        integer :: triangle

        valid = .true.
        do triangle = 1, mesh%n_triangles
            determinant = &
                (mesh%vertices(1, mesh%triangles(2, triangle)) - &
                mesh%vertices(1, mesh%triangles(1, triangle))) * &
                (mesh%vertices(2, mesh%triangles(3, triangle)) - &
                mesh%vertices(2, mesh%triangles(1, triangle))) - &
                (mesh%vertices(1, mesh%triangles(3, triangle)) - &
                mesh%vertices(1, mesh%triangles(1, triangle))) * &
                (mesh%vertices(2, mesh%triangles(2, triangle)) - &
                mesh%vertices(2, mesh%triangles(1, triangle)))
            if (determinant <= 0.0_dp) then
                valid = .false.
                return
            end if
        end do
    end function counterclockwise_mesh

    subroutine export_real_csc( &
            matrix, nnz_capacity, n_dofs, nnz, &
            col_ptr, row_ind, values, status)
        type(csc_t), intent(in) :: matrix
        integer(c_int), intent(in) :: nnz_capacity
        integer(c_int), intent(out) :: n_dofs, nnz
        integer(c_int), intent(out) :: col_ptr(*), row_ind(*)
        real(c_double), intent(out) :: values(*)
        integer(c_int), intent(out) :: status

        integer :: entry

        n_dofs = int(matrix%ncol, c_int)
        nnz = int(matrix%nnz, c_int)
        if (nnz_capacity < nnz) then
            status = -2_c_int
            return
        end if
        do entry = 1, matrix%ncol + 1
            col_ptr(entry) = int(matrix%col_ptr(entry) - 1, c_int)
        end do
        do entry = 1, matrix%nnz
            row_ind(entry) = int(matrix%row_idx(entry) - 1, c_int)
            values(entry) = real(matrix%val(entry), c_double)
        end do
        status = 0_c_int
    end subroutine export_real_csc

    subroutine copy_complex_dofs(n_dofs, c_dofs, fortran_dofs)
        integer(c_int), intent(in) :: n_dofs
        complex(c_double_complex), intent(in) :: c_dofs(*)
        complex(dp), allocatable, intent(out) :: fortran_dofs(:)

        integer :: dof

        allocate(fortran_dofs(int(n_dofs)))
        do dof = 1, int(n_dofs)
            fortran_dofs(dof) = cmplx( &
                real(c_dofs(dof), dp), aimag(c_dofs(dof)), dp)
        end do
    end subroutine copy_complex_dofs

    pure complex(c_double_complex) function to_c_complex(value) result(c_value)
        complex(dp), intent(in) :: value

        c_value = cmplx( &
            real(value, c_double), aimag(value), c_double_complex)
    end function to_c_complex

end module fortfem_capi_mesh
