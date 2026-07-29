program test_capi_nedelec_field
    use check, only: check_condition, check_summary
    use iso_c_binding, only: c_double, c_double_complex, c_int
    implicit none

    real(c_double), parameter :: vertices(8) = [ &
        0.0_c_double, 0.0_c_double, 1.0_c_double, 0.0_c_double, &
        1.0_c_double, 1.0_c_double, 0.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: triangles(6) = [ &
        0_c_int, 1_c_int, 2_c_int, 0_c_int, 2_c_int, 3_c_int]
    complex(c_double_complex), parameter :: scale = &
        cmplx(1.0, 2.0, c_double_complex)
    integer(c_int) :: edges(10), triangle_edge_dofs(6)
    integer(c_int) :: triangle_edge_signs(6)
    complex(c_double_complex) :: curl, dofs(5), value(2)
    integer(c_int) :: edge, handle, n_edges, status
    real(c_double) :: edge_R, edge_Z, midpoint_R, midpoint_Z
    integer :: vertex_a, vertex_b
    logical :: all_passed

    interface
        subroutine fortfem_triangle_mesh_create( &
                n_vertices, vertices, n_triangles, triangles, &
                handle, n_edges, status) &
                bind(c, name="fortfem_triangle_mesh_create")
            import c_double, c_int
            integer(c_int), value :: n_vertices, n_triangles
            real(c_double), intent(in) :: vertices(*)
            integer(c_int), intent(in) :: triangles(*)
            integer(c_int), intent(out) :: handle, n_edges, status
        end subroutine fortfem_triangle_mesh_create

        subroutine fortfem_triangle_mesh_edges( &
                handle, edge_capacity, n_edges, edges, &
                triangle_edge_dofs, triangle_edge_signs, status) &
                bind(c, name="fortfem_triangle_mesh_edges")
            import c_int
            integer(c_int), value :: handle, edge_capacity
            integer(c_int), intent(out) :: n_edges, edges(*)
            integer(c_int), intent(out) :: triangle_edge_dofs(*)
            integer(c_int), intent(out) :: triangle_edge_signs(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_triangle_mesh_edges

        subroutine fortfem_nedelec_evaluate( &
                handle, triangle, x, y, n_dofs, dofs, &
                value, curl, status) bind(c, name="fortfem_nedelec_evaluate")
            import c_double, c_double_complex, c_int
            integer(c_int), value :: handle, triangle, n_dofs
            real(c_double), value :: x, y
            complex(c_double_complex), intent(in) :: dofs(*)
            complex(c_double_complex), intent(out) :: value(*), curl
            integer(c_int), intent(out) :: status
        end subroutine fortfem_nedelec_evaluate

        subroutine fortfem_triangle_mesh_free(handle, status) &
                bind(c, name="fortfem_triangle_mesh_free")
            import c_int
            integer(c_int), value :: handle
            integer(c_int), intent(out) :: status
        end subroutine fortfem_triangle_mesh_free
    end interface

    all_passed = .true.
    call fortfem_triangle_mesh_create( &
        4_c_int, vertices, 2_c_int, triangles, handle, n_edges, status)
    call fortfem_triangle_mesh_edges( &
        handle, 5_c_int, n_edges, edges, triangle_edge_dofs, &
        triangle_edge_signs, status)
    do edge = 0_c_int, n_edges - 1_c_int
        vertex_a = edges(2 * edge + 1)
        vertex_b = edges(2 * edge + 2)
        edge_R = vertices(2 * vertex_b + 1) - vertices(2 * vertex_a + 1)
        edge_Z = vertices(2 * vertex_b + 2) - vertices(2 * vertex_a + 2)
        midpoint_R = 0.5_c_double * ( &
            vertices(2 * vertex_a + 1) + vertices(2 * vertex_b + 1))
        midpoint_Z = 0.5_c_double * ( &
            vertices(2 * vertex_a + 2) + vertices(2 * vertex_b + 2))
        dofs(edge + 1) = scale * ( &
            (1.0_c_double + 2.0_c_double * midpoint_Z) * edge_R + &
            (3.0_c_double - 2.0_c_double * midpoint_R) * edge_Z)
    end do

    call fortfem_nedelec_evaluate( &
        handle, 0_c_int, 0.75_c_double, 0.25_c_double, n_edges, &
        dofs, value, curl, status)
    call record_condition(status == 0_c_int .and. &
        maxval(abs(value - scale * [1.5_c_double, 1.5_c_double])) < &
        1.0e-13_c_double, &
        "C Nedelec evaluation reproduces an affine complex field")
    call record_condition(abs(curl + 4.0_c_double * scale) < &
        1.0e-13_c_double, "C Nedelec evaluation returns exact curl")

    call fortfem_triangle_mesh_free(handle, status)
    call check_summary("C Nedelec field API")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_nedelec_field
