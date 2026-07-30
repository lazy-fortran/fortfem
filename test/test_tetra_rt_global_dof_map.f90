program test_tetra_rt_global_dof_map
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        build_tetra_rt_basis_transform, build_tetra_rt_dof_map
    use fortfem_tetra_face_moment_transforms, only: &
        transform_tetra_rt_face_moments
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 6
    integer, parameter :: dof_count = 280
    integer, parameter :: face_dof_count = 28
    integer, allocatable :: face_orientations(:, :)
    integer, allocatable :: face_permutations(:, :, :), faces(:, :)
    integer, allocatable :: global_dofs(:, :)
    real(dp) :: canonical(face_dof_count), local(face_dof_count)
    real(dp) :: recovered(face_dof_count), transform(dof_count, dof_count)
    integer :: face, face_start, moment, status, tetrahedra(4, 2)
    logical :: all_passed

    all_passed = .true.
    tetrahedra(:, 1) = [1, 2, 3, 4]
    tetrahedra(:, 2) = [3, 2, 1, 5]
    call build_tetra_rt_dof_map( &
        degree, tetrahedra, faces, global_dofs, face_orientations, &
        face_permutations, status)

    call record_condition(status == 0, &
        "Arbitrary-order tetrahedral RT topology builds")
    call record_condition(size(faces, 2) == 7, &
        "Two tetrahedra share one canonical RT face")
    call record_condition(size(global_dofs, 1) == dof_count .and. &
        maxval(global_dofs) == 532, &
        "Degree-six RT topology has the analytic global dimension")
    call record_condition(all( &
        global_dofs(2*face_dof_count + 1:3*face_dof_count, 1) == &
        global_dofs(2*face_dof_count + 1:3*face_dof_count, 2)), &
        "Every shared-face RT moment has one global degree of freedom")
    call record_condition( &
        face_orientations(3, 1) == -face_orientations(3, 2), &
        "Adjacent tetrahedra use opposite outward normal orientations")

    do moment = 1, face_dof_count
        canonical(moment) = real(3*moment - 8, dp)/real(moment + 2, dp)
    end do
    call build_tetra_rt_basis_transform( &
        degree, face_orientations(:, 2), face_permutations(:, :, 2), &
        transform, status)
    call record_condition(status == 0, &
        "Degree-six RT basis orientation transform builds")
    do face = 1, 4
        face_start = (face - 1)*face_dof_count
        local = matmul( &
            transform( &
            face_start + 1:face_start + face_dof_count, &
            face_start + 1:face_start + face_dof_count), canonical)
        call transform_tetra_rt_face_moments( &
            degree, face_permutations(:, face, 2), local, recovered, status)
        recovered = face_orientations(face, 2)*recovered
        if (maxval(abs(recovered - canonical)) >= 1.0e-8_dp) then
            write (*, '(A,I0,A,ES12.4)') "RT face ", face, &
                " orientation error ", maxval(abs(recovered - canonical))
        end if
        call record_condition(status == 0 .and. &
            maxval(abs(recovered - canonical)) < 1.0e-8_dp, &
            "RT element transform gives canonical oriented face moments")
    end do

    call build_tetra_rt_dof_map( &
        -1, tetrahedra, faces, global_dofs, face_orientations, &
        face_permutations, status)
    call record_condition(status /= 0, "RT topology rejects negative degree")
    tetrahedra(:, 2) = [1, 1, 2, 5]
    call build_tetra_rt_dof_map( &
        degree, tetrahedra, faces, global_dofs, face_orientations, &
        face_permutations, status)
    call record_condition(status /= 0, &
        "RT topology rejects a tetrahedron with repeated vertices")

    call check_summary("Arbitrary-order tetrahedral RT topology")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_rt_global_dof_map
