program test_freefem_maxwell_reference
    use check, only: check_condition, check_summary
    use fortfem_assembly_mixed_2d, only: assemble_nedelec_rt_mass_csc
    use fortfem_assembly_nedelec_2d, only: &
        assemble_nedelec_axisymmetric_fourier
    use fortfem_edge_interpolation_2d, only: interpolate_rt_edge_dofs
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_nedelec_field_2d, only: evaluate_nedelec_field_2d
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    use iso_c_binding, only: c_double, c_double_complex, c_int
    implicit none

    real(c_double), parameter :: c_vertices(18) = [ &
        2.0_c_double, 0.0_c_double, 2.5_c_double, 0.0_c_double, &
        3.0_c_double, 0.0_c_double, 2.0_c_double, 0.5_c_double, &
        2.5_c_double, 0.5_c_double, 3.0_c_double, 0.5_c_double, &
        2.0_c_double, 1.0_c_double, 2.5_c_double, 1.0_c_double, &
        3.0_c_double, 1.0_c_double]
    integer(c_int), parameter :: c_triangles(24) = [ &
        0_c_int, 1_c_int, 4_c_int, 0_c_int, 4_c_int, 3_c_int, &
        1_c_int, 2_c_int, 5_c_int, 1_c_int, 5_c_int, 4_c_int, &
        3_c_int, 4_c_int, 7_c_int, 3_c_int, 7_c_int, 6_c_int, &
        4_c_int, 5_c_int, 8_c_int, 4_c_int, 8_c_int, 7_c_int]
    real(dp), parameter :: reference_matrix(8, 8) = reshape([ &
        41.069535155375213_dp, -18.660516146451499_dp, 0.0_dp, &
        -21.338036046717779_dp, 21.083731882452010_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, -18.660516146451499_dp, 36.595048348973911_dp, &
        -17.323145412045168_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, -17.323145412045168_dp, 37.189111924768532_dp, 0.0_dp, &
        0.0_dp, 18.379328829614302_dp, -18.676360537521308_dp, 0.0_dp, &
        -21.338036046717779_dp, 0.0_dp, 0.0_dp, 44.486190410572519_dp, &
        -21.326558275342396_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        21.083731882452010_dp, 0.0_dp, 0.0_dp, -21.326558275342396_dp, &
        44.971843148581499_dp, 0.0_dp, 0.0_dp, -22.673172281430439_dp, &
        0.0_dp, 0.0_dp, 18.379328829614302_dp, 0.0_dp, 0.0_dp, &
        41.069535155375213_dp, -18.660516146451499_dp, &
        -21.338036046717779_dp, 0.0_dp, 0.0_dp, -18.676360537521308_dp, &
        0.0_dp, 0.0_dp, -18.660516146451499_dp, 36.595048348973911_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -22.673172281430439_dp, &
        -21.338036046717779_dp, 0.0_dp, 44.486190410572519_dp], [8, 8])
    real(dp), parameter :: reference_rhs(8) = [ &
        4.083333333333333_dp, 2.833333333333333_dp, &
        -1.0000000000000002_dp, 3.2500000000000013_dp, &
        -1.1666666666666681_dp, 4.7500000000000009_dp, &
        3.4166666666666670_dp, 3.8333333333333344_dp]
    real(dp), parameter :: reference_solution(8) = [ &
        0.11007462701844933_dp, 0.075735449754011447_dp, &
        -0.12213949558924840_dp, 0.28756367201823718_dp, &
        0.33731818404732899_dp, 0.61357589369381260_dp, &
        0.34390413048519836_dp, 0.55239414623134686_dp]
    real(dp), parameter :: reference_nodal(2, 9) = reshape([ &
        0.15147089950802289_dp, 0.0_dp, &
        0.35497808999957570_dp, 0.22014925403689867_dp, &
        0.0_dp, 0.0_dp, 0.68780826097039671_dp, 0.0_dp, &
        -0.12236349492493148_dp, 1.2271517873876252_dp, &
        0.67463636809465799_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.2271517873876252_dp, 0.0_dp, &
        1.1047882924626937_dp], [2, 9])

    type(csc_t) :: mixed_matrix
    type(fortsparse_status_t) :: sparse_status
    type(mesh_2d_t) :: mesh
    complex(dp) :: field_dofs(16), nodal_value(2), residual_curl
    real(dp) :: matrix(16, 16), current_dofs(16), load(16)
    real(dp) :: nodal_error, residual(8)
    integer(c_int) :: handle, n_edges, status
    integer :: point, triangle
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

        subroutine fortfem_nedelec_evaluate_vertex( &
                handle, vertex, n_dofs, dofs, value, curl, status) &
                bind(c, name="fortfem_nedelec_evaluate_vertex")
            import c_double_complex, c_int
            integer(c_int), value :: handle, vertex, n_dofs
            complex(c_double_complex), intent(in) :: dofs(*)
            complex(c_double_complex), intent(out) :: value(2), curl
            integer(c_int), intent(out) :: status
        end subroutine fortfem_nedelec_evaluate_vertex

        subroutine fortfem_triangle_mesh_free(handle, status) &
                bind(c, name="fortfem_triangle_mesh_free")
            import c_int
            integer(c_int), value :: handle
            integer(c_int), intent(out) :: status
        end subroutine fortfem_triangle_mesh_free
    end interface

    all_passed = .true.
    call mesh%create_rectangular(3, 3, 2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()

    call assemble_nedelec_axisymmetric_fourier(mesh, 2, 5, matrix)
    call record_condition( &
        maxval(abs(matrix(:8, :8) - reference_matrix)) < 5.0e-12_dp, &
        "FortFEM interior operator matches FreeFem 4.15")

    call interpolate_rt_edge_dofs(mesh, current_field, 2, current_dofs)
    call assemble_nedelec_rt_mass_csc(mesh, 2, mixed_matrix, sparse_status)
    load = csc_matvec(mixed_matrix, current_dofs)
    call record_condition(sparse_status%code == 0, &
        "Mixed Nedelec-RT0 reference operator assembles")
    call record_condition(maxval(abs(load(:8) - reference_rhs)) < 5.0e-13_dp, &
        "FortFEM mixed load matches FreeFem 4.15")

    residual = matmul(matrix(:8, :8), reference_solution) - reference_rhs
    call record_condition(maxval(abs(residual)) < 5.0e-12_dp, &
        "FreeFem solution satisfies the FortFEM interior system")

    field_dofs = cmplx(0.0_dp, 0.0_dp, dp)
    field_dofs(:8) = cmplx(reference_solution, 0.0_dp, dp)
    nodal_error = 0.0_dp
    do point = 1, mesh%n_vertices
        triangle = last_incident_triangle(mesh, point)
        call evaluate_nedelec_field_2d(mesh, triangle, &
            mesh%vertices(1, point), mesh%vertices(2, point), &
            field_dofs, nodal_value)
        nodal_error = max(nodal_error, &
            maxval(abs(real(nodal_value, dp) - reference_nodal(:, point))))
    end do
    call record_condition(nodal_error < 5.0e-12_dp, &
        "FortFEM nodal evaluation matches FreeFem interpolation")

    call fortfem_triangle_mesh_create( &
        9_c_int, c_vertices, 8_c_int, c_triangles, handle, n_edges, status)
    call record_condition(status == 0_c_int .and. n_edges == 16_c_int, &
        "FreeFem reference mesh is retained through the C API")
    nodal_error = 0.0_dp
    do point = 1, 9
        call fortfem_nedelec_evaluate_vertex( &
            handle, int(point - 1, c_int), 16_c_int, field_dofs, &
            nodal_value, residual_curl, status)
        nodal_error = max(nodal_error, &
            maxval(abs(real(nodal_value, dp) - reference_nodal(:, point))))
    end do
    call record_condition(status == 0_c_int .and. nodal_error < 5.0e-12_dp, &
        "Retained-mesh vertex evaluation matches FreeFem interpolation")
    call fortfem_triangle_mesh_free(handle, status)

    call check_summary("FreeFem Maxwell reference")
    if (.not. all_passed) error stop 1

contains

    pure subroutine current_field(x, y, value)
        real(dp), intent(in) :: x, y
        real(dp), intent(out) :: value(2)

        value = [1.0_dp + x + 2.0_dp * y, 3.0_dp + 4.0_dp * x + 5.0_dp * y]
    end subroutine current_field

    integer function last_incident_triangle(local_mesh, point) result(triangle)
        type(mesh_2d_t), intent(in) :: local_mesh
        integer, intent(in) :: point

        integer :: candidate

        triangle = 0
        do candidate = 1, local_mesh%n_triangles
            if (any(local_mesh%triangles(:, candidate) == point)) then
                triangle = candidate
            end if
        end do
        if (triangle == 0) error stop "Reference vertex has no incident triangle"
    end function last_incident_triangle

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_freefem_maxwell_reference
