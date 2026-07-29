program test_capi_rt0
    use iso_c_binding, only: c_double, c_double_complex, c_int
    use check, only: check_condition, check_summary
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
    complex(c_double_complex) :: dofs(5), value(2), divergence
    complex(c_double_complex) :: toroidal(2)
    integer(c_int) :: handle, n_edges, status, edge_id
    integer :: vertex_a, vertex_b
    real(c_double) :: point_a(2), point_b(2), edge_vector(2), norm
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

        subroutine fortfem_rt0_evaluate( &
                handle, triangle, x, y, n_dofs, dofs, &
                value, divergence, status) bind(c, name="fortfem_rt0_evaluate")
            import c_double, c_double_complex, c_int
            integer(c_int), value :: handle, triangle, n_dofs
            real(c_double), value :: x, y
            complex(c_double_complex), intent(in) :: dofs(*)
            complex(c_double_complex), intent(out) :: value(*), divergence
            integer(c_int), intent(out) :: status
        end subroutine fortfem_rt0_evaluate

        subroutine fortfem_rt0_l2_norm( &
                handle, n_dofs, dofs, norm, status) &
                bind(c, name="fortfem_rt0_l2_norm")
            import c_double, c_double_complex, c_int
            integer(c_int), value :: handle, n_dofs
            complex(c_double_complex), intent(in) :: dofs(*)
            real(c_double), intent(out) :: norm
            integer(c_int), intent(out) :: status
        end subroutine fortfem_rt0_l2_norm

        subroutine fortfem_rt0_toroidal( &
                handle, mode, n_dofs, dofs, toroidal, status) &
                bind(c, name="fortfem_rt0_toroidal")
            import c_double_complex, c_int
            integer(c_int), value :: handle, mode, n_dofs
            complex(c_double_complex), intent(in) :: dofs(*)
            complex(c_double_complex), intent(out) :: toroidal(*)
            integer(c_int), intent(out) :: status
        end subroutine fortfem_rt0_toroidal

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

    do edge_id = 0_c_int, n_edges - 1_c_int
        vertex_a = edges(2 * edge_id + 1)
        vertex_b = edges(2 * edge_id + 2)
        point_a = vertices(2 * vertex_a + 1:2 * vertex_a + 2)
        point_b = vertices(2 * vertex_b + 1:2 * vertex_b + 2)
        edge_vector = point_b - point_a
        dofs(edge_id + 1) = scale * ( &
            edge_vector(2) * (1.0_c_double + point_a(1) + point_b(1)) - &
            edge_vector(1) * (3.0_c_double + point_a(2) + point_b(2)))
    end do

    call fortfem_rt0_evaluate(handle, 0_c_int, 0.75_c_double, &
        0.25_c_double, n_edges, dofs, value, divergence, status)
    call record_condition(status == 0_c_int .and. maxval(abs(value - &
        scale * [2.5_c_double, 3.5_c_double])) < 1.0e-13_c_double, &
        "C RT0 evaluation reproduces an affine complex field")
    call record_condition(abs(divergence - 4.0_c_double * scale) < &
        1.0e-13_c_double, "C RT0 evaluation returns exact divergence")

    call fortfem_rt0_l2_norm(handle, n_edges, dofs, norm, status)
    call record_condition(status == 0_c_int .and. abs(norm - &
        sqrt(5.0_c_double * 62.0_c_double / 3.0_c_double)) < &
        1.0e-13_c_double, "C RT0 L2 norm matches an analytical integral")

    call fortfem_rt0_toroidal( &
        handle, 2_c_int, n_edges, dofs, toroidal, status)
    call record_condition(status == 0_c_int .and. maxval(abs(toroidal - &
        cmplx(-4.0, 2.0, c_double_complex))) < 1.0e-13_c_double, &
        "C RT0 toroidal coefficients enforce zero divergence")

    call fortfem_triangle_mesh_free(handle, status)
    call check_summary("C RT0 field API")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_capi_rt0
