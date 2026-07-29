program test_nedelec_local_assembly
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_assembly_nedelec_2d, only: &
        assemble_nedelec_curl_mass_element
    use check, only: check_condition, check_summary
    implicit none

    type(mesh_2d_t) :: mesh
    real(dp) :: element_matrix(3, 3)
    real(dp), parameter :: expected(3, 3) = reshape([ &
        7.0_dp / 3.0_dp, 2.0_dp, 11.0_dp / 6.0_dp, &
        2.0_dp, 13.0_dp / 6.0_dp, 2.0_dp, &
        11.0_dp / 6.0_dp, 2.0_dp, 7.0_dp / 3.0_dp], [3, 3])
    logical :: all_passed

    all_passed = .true.
    mesh%n_vertices = 3
    mesh%n_triangles = 1
    allocate(mesh%vertices(2, 3), mesh%triangles(3, 1))
    mesh%vertices(:, 1) = [0.0_dp, 0.0_dp]
    mesh%vertices(:, 2) = [1.0_dp, 0.0_dp]
    mesh%vertices(:, 3) = [0.0_dp, 1.0_dp]
    mesh%triangles(:, 1) = [1, 2, 3]

    call assemble_nedelec_curl_mass_element(mesh, 1, element_matrix)

    call record_condition( &
        maxval(abs(element_matrix - expected)) < 1.0e-13_dp, &
        "Reference curl-curl plus mass matrix matches exact integration")
    call record_condition( &
        maxval(abs(element_matrix - transpose(element_matrix))) < &
        1.0e-13_dp, &
        "Nedelec element operator is symmetric")

    call check_summary("Nedelec local assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_nedelec_local_assembly
