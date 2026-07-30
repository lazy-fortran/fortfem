program test_tetra_nedelec_global_dof_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_face_moment_transforms, only: &
        transform_tetra_face_moments
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: order = 5
    integer, parameter :: dof_count = 140
    integer, parameter :: face_dof_count = 20
    integer, allocatable :: edge_orientations(:, :), edges(:, :)
    integer, allocatable :: face_permutations(:, :, :), faces(:, :)
    integer, allocatable :: global_dofs(:, :)
    real(dp) :: canonical(face_dof_count), recovered(face_dof_count)
    real(dp) :: local(face_dof_count), transform(dof_count, dof_count)
    integer :: face_start, moment, status, tetrahedra(4, 2)
    logical :: all_passed

    all_passed = .true.
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [3, 1, 2, 5]
    call build_tetra_nedelec_dof_map( &
        order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)

    call record_condition(status == 0, &
        "Arbitrary-order tetrahedral topology builds")
    call record_condition(size(edges, 2) == 9 .and. size(faces, 2) == 7, &
        "Shared face reuses three edges and one canonical face")
    call record_condition(size(global_dofs, 1) == dof_count .and. &
        maxval(global_dofs) == 245, &
        "Order-five topology has the analytic global dimension")
    face_start = 6 * order + 1
    call record_condition(all( &
        global_dofs(face_start:face_start + face_dof_count - 1, 1) == &
        global_dofs(face_start:face_start + face_dof_count - 1, 2)), &
        "Every shared-face moment has one global degree of freedom")
    call record_condition(all(face_permutations(:, 1, 1) == [1, 2, 3]) .and. &
        all(face_permutations(:, 1, 2) == [3, 1, 2]), &
        "Face permutations encode local vertex order against sorted vertices")

    do moment = 1, face_dof_count
        canonical(moment) = real(2 * moment - 7, dp) / real(moment + 2, dp)
    end do
    call build_tetra_nedelec_basis_transform( &
        order, edge_orientations(:, 2), face_permutations(:, :, 2), &
        transform, status)
    call record_condition( &
        all([transform(1, 1), transform(2, 2), &
        transform(3, 3), transform(4, 4)] == &
        [-1.0_dp, 1.0_dp, -1.0_dp, 1.0_dp]), &
        "Reversed edge uses the analytic alternating moment signs")
    local = matmul( &
        transform( &
        face_start:face_start + face_dof_count - 1, &
        face_start:face_start + face_dof_count - 1), canonical)
    call transform_tetra_face_moments( &
        order, face_permutations(:, 1, 2), local, recovered, status)
    call record_condition(status == 0 .and. &
        maxval(abs(recovered - canonical)) < 2.0e-11_dp, &
        "Element transform gives shared faces canonical moment coefficients")

    call build_tetra_nedelec_dof_map( &
        0, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)
    call record_condition(status /= 0, "Topology rejects order zero")
    tetrahedra(:, 2) = [1, 1, 2, 5]
    call build_tetra_nedelec_dof_map( &
        order, tetrahedra, edges, faces, global_dofs, edge_orientations, &
        face_permutations, status)
    call record_condition(status /= 0, &
        "Topology rejects a tetrahedron with repeated vertices")

    call check_summary("Arbitrary-order tetrahedral Nedelec topology")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_global_dof_map
