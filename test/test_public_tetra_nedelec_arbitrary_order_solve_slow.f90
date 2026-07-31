program test_public_tetra_nedelec_arbitrary_order_solve_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_tetra_nedelec_curl_mass
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, tetra_nedelec_first_kind_t
    use fortfem_tetra_nedelec_global_dof_map, only: &
        build_tetra_nedelec_basis_transform, build_tetra_nedelec_dof_map
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: tetrahedra(4, 2) = reshape( &
        [1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: vertices(3, 5) = reshape( &
        [0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:)
    real(dp) :: maximum_error
    integer :: order

    do order = 1, 6
        call solve_tetra_nedelec_curl_mass( &
            vertices, tetrahedra, order, constant_source, &
            1.0_dp, 1.0_dp, solution, status)
        if (status%code /= 0) error stop "public tetra Nedelec solve failed"
        call measure_constant_error(order, solution, maximum_error)
        call check_condition(maximum_error < solve_tolerance(order), &
            "Public tetrahedral Hcurl solve reproduces a constant field")
        deallocate(solution)
    end do
    call check_summary("Public tetrahedral arbitrary-order Nedelec solve")

contains

    pure real(dp) function solve_tolerance(order) result(tolerance)
        integer, intent(in) :: order

        tolerance = merge(3.0e-10_dp, 2.0e-8_dp, order <= 4)
    end function solve_tolerance

    pure subroutine constant_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)

        associate(unused => x + y + z)
        end associate
        value = [1.25_dp, -0.75_dp, 0.5_dp]
    end subroutine constant_source

    subroutine measure_constant_error(order, coefficients, error)
        integer, intent(in) :: order
        real(dp), intent(in) :: coefficients(:)
        real(dp), intent(out) :: error

        type(tetra_nedelec_first_kind_t) :: basis
        integer, allocatable :: edge_orientations(:, :), edges(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), local_dofs(:)
        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: jacobian(3, 3), value(3)
        integer :: dof_count, local_status, tetrahedron

        call build_tetra_nedelec_dof_map( &
            order, tetrahedra, edges, faces, global_dofs, &
            edge_orientations, face_permutations, local_status)
        if (local_status /= 0) error stop "tetra Hcurl map failed"
        call initialize_tetra_nedelec_first_kind(order, basis, local_status)
        if (local_status /= 0) error stop "tetra Hcurl basis failed"
        dof_count = size(global_dofs, 1)
        allocate(basis_transform(dof_count, dof_count))
        allocate(local_dofs(dof_count))
        allocate(reference_values(3, dof_count), reference_curls(3, dof_count))
        allocate(physical_values(3, dof_count), physical_curls(3, dof_count))
        error = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            call build_tetra_nedelec_basis_transform( &
                order, edge_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) error stop "tetra Hcurl transform failed"
            local_dofs = matmul( &
                basis_transform, coefficients(global_dofs(:, tetrahedron)))
            jacobian(:, 1) = vertices(:, tetrahedra(2, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 2) = vertices(:, tetrahedra(3, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 3) = vertices(:, tetrahedra(4, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            call evaluate_tetra_nedelec_first_kind( &
                basis, [0.21_dp, 0.24_dp, 0.19_dp], reference_values, &
                reference_curls, local_status)
            if (local_status /= 0) error stop "tetra Hcurl evaluation failed"
            call map_tetra_nedelec_covariant( &
                jacobian, reference_values, reference_curls, physical_values, &
                physical_curls, local_status)
            if (local_status /= 0) error stop "tetra Hcurl Piola map failed"
            value = matmul(physical_values, local_dofs)
            error = max(error, maxval(abs( &
                value - [1.25_dp, -0.75_dp, 0.5_dp])))
            error = max(error, maxval(abs(matmul(physical_curls, local_dofs))))
        end do
    end subroutine measure_constant_error

end program test_public_tetra_nedelec_arbitrary_order_solve_slow
