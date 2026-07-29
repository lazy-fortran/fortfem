program test_nedelec_lowest_order
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_basis_edge_2d, only: evaluate_edge_basis_2d, &
        evaluate_edge_basis_curl_2d, evaluate_edge_basis_2d_piola
    use check, only: check_condition, check_summary
    implicit none

    logical :: all_passed = .true.

    call test_reference_edge_moments()
    call test_reference_curls()
    call test_covariant_piola_edge_moments()
    call check_summary("Lowest-order Nedelec element")
    if (.not. all_passed) error stop 1

contains

    subroutine test_reference_edge_moments()
        real(dp), parameter :: midpoint(2, 3) = reshape([ &
            0.5_dp, 0.0_dp, &
            0.5_dp, 0.5_dp, &
            0.0_dp, 0.5_dp], [2, 3])
        real(dp), parameter :: edge_vector(2, 3) = reshape([ &
            1.0_dp, 0.0_dp, &
            -1.0_dp, 1.0_dp, &
            0.0_dp, -1.0_dp], [2, 3])
        real(dp) :: values(2, 3), moment
        integer :: basis_id, edge_id

        do edge_id = 1, 3
            call evaluate_edge_basis_2d(midpoint(1, edge_id), &
                midpoint(2, edge_id), 0.5_dp, values)
            do basis_id = 1, 3
                moment = dot_product(values(:, basis_id), &
                    edge_vector(:, edge_id))
                call record_condition(abs(moment - kronecker_delta( &
                    basis_id, edge_id)) < 1.0e-13_dp, &
                    "Nedelec reference edge moments are unisolvent")
            end do
        end do
    end subroutine test_reference_edge_moments

    subroutine test_reference_curls()
        real(dp) :: curls(3)

        call evaluate_edge_basis_curl_2d(0.23_dp, 0.31_dp, 0.5_dp, curls)
        call record_condition(maxval(abs(curls - 2.0_dp)) < 1.0e-13_dp, &
            "Nedelec reference basis has constant curl two")
    end subroutine test_reference_curls

    subroutine test_covariant_piola_edge_moments()
        real(dp), parameter :: midpoint(2, 3) = reshape([ &
            0.5_dp, 0.0_dp, &
            0.5_dp, 0.5_dp, &
            0.0_dp, 0.5_dp], [2, 3])
        type(mesh_2d_t) :: mesh
        real(dp) :: values(2, 3), physical_edge(2), moment
        integer :: basis_id, edge_id, start_vertex, end_vertex

        mesh%n_vertices = 3
        mesh%n_triangles = 1
        allocate(mesh%vertices(2, 3), mesh%triangles(3, 1))
        mesh%vertices(:, 1) = [0.2_dp, -0.3_dp]
        mesh%vertices(:, 2) = [2.1_dp, 0.4_dp]
        mesh%vertices(:, 3) = [-0.4_dp, 1.7_dp]
        mesh%triangles(:, 1) = [1, 2, 3]

        do edge_id = 1, 3
            start_vertex = edge_id
            end_vertex = mod(edge_id, 3) + 1
            physical_edge = mesh%vertices(:, end_vertex) - &
                mesh%vertices(:, start_vertex)
            call evaluate_edge_basis_2d_piola(mesh, 1, &
                midpoint(1, edge_id), midpoint(2, edge_id), values)
            do basis_id = 1, 3
                moment = dot_product(values(:, basis_id), physical_edge)
                call record_condition(abs(moment - kronecker_delta( &
                    basis_id, edge_id)) < 1.0e-13_dp, &
                    "Covariant Piola map preserves tangential moments")
            end do
        end do
    end subroutine test_covariant_piola_edge_moments

    pure real(dp) function kronecker_delta(i, j) result(value)
        integer, intent(in) :: i, j

        value = 0.0_dp
        if (i == j) value = 1.0_dp
    end function kronecker_delta

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_lowest_order
