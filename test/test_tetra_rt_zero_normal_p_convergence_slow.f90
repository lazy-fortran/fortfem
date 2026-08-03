program test_tetra_rt_zero_normal_p_convergence_slow
    use check, only: check_condition, check_summary
    use fortfem_feec, only: solve_tetra_rt_div_mass
    use fortfem_kinds, only: dp
    use fortfem_tetra_piola_maps, only: map_tetra_rt_contravariant
    use fortfem_tetra_rt_arbitrary_order, only: &
        evaluate_tetra_rt, initialize_tetra_rt, tetra_rt_t
    use fortfem_tetra_rt_global_dof_map, only: &
        build_tetra_rt_basis_transform, build_tetra_rt_dof_map
    use fortnum_linalg, only: det3
    use fortsparse, only: fortsparse_status_t
    implicit none

    type(fortsparse_status_t) :: status
    real(dp), allocatable :: solution(:)
    real(dp) :: errors(2, 6), vertices(3, 8)
    integer :: degree, tetrahedra(4, 6)

    call build_cube_mesh(vertices, tetrahedra)
    do degree = 0, 5
        call solve_tetra_rt_div_mass( &
            vertices, tetrahedra, degree, vector_source, &
            1.0_dp, 1.0_dp, solution, status, .true.)
        if (status%code /= 0) error stop "zero-normal tetra Hdiv solve failed"
        call measure_error( &
            vertices, tetrahedra, degree, solution, errors(:, degree + 1))
        write(*, '(a,i0,a,2(es12.4,1x))') &
            "zero-normal RT degree ", degree, " field/div errors ", &
            errors(:, degree + 1)
        deallocate(solution)
    end do

    call check_condition(all(errors(1, 2:6) < errors(1, 1:5)), &
        "zero-normal field error decreases with degree")
    call check_condition(errors(2, 6) < 0.01_dp*errors(2, 1), &
        "zero-normal divergence converges by two orders of magnitude")
    call check_condition(errors(1, 6) < 1.0e-2_dp, &
        "degree-five solve reaches the analytical Hdiv field")
    call check_condition(errors(2, 6) < 5.0e-2_dp, &
        "degree-five solve reaches the analytical divergence")
    call check_summary("Tetrahedral zero-normal Hdiv p-convergence")

contains

    pure subroutine vector_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = (pi*pi + 1.0_dp)*[ &
            sin(pi*x), sin(pi*y), sin(pi*z)]
    end subroutine vector_source

    pure function exact_field(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = sin(pi*point)
    end function exact_field

    pure function exact_divergence(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = pi*sum(cos(pi*point))
    end function exact_divergence

    subroutine build_cube_mesh(mesh_vertices, cells)
        real(dp), intent(out) :: mesh_vertices(3, 8)
        integer, intent(out) :: cells(4, 6)

        real(dp) :: jacobian(3, 3)
        integer :: cell, temporary

        mesh_vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        mesh_vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        mesh_vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        mesh_vertices(:, 4) = [1.0_dp, 1.0_dp, 0.0_dp]
        mesh_vertices(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        mesh_vertices(:, 6) = [1.0_dp, 0.0_dp, 1.0_dp]
        mesh_vertices(:, 7) = [0.0_dp, 1.0_dp, 1.0_dp]
        mesh_vertices(:, 8) = [1.0_dp, 1.0_dp, 1.0_dp]
        cells(:, 1) = [1, 2, 4, 8]
        cells(:, 2) = [1, 2, 6, 8]
        cells(:, 3) = [1, 3, 4, 8]
        cells(:, 4) = [1, 3, 7, 8]
        cells(:, 5) = [1, 5, 6, 8]
        cells(:, 6) = [1, 5, 7, 8]
        do cell = 1, size(cells, 2)
            jacobian(:, 1) = &
                mesh_vertices(:, cells(2, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 2) = &
                mesh_vertices(:, cells(3, cell)) - &
                mesh_vertices(:, cells(1, cell))
            jacobian(:, 3) = &
                mesh_vertices(:, cells(4, cell)) - &
                mesh_vertices(:, cells(1, cell))
            if (det3(jacobian) < 0.0_dp) then
                temporary = cells(3, cell)
                cells(3, cell) = cells(4, cell)
                cells(4, cell) = temporary
            end if
        end do
    end subroutine build_cube_mesh

    subroutine measure_error( &
            mesh_vertices, cells, degree, coefficients, error)
        real(dp), intent(in) :: mesh_vertices(:, :), coefficients(:)
        integer, intent(in) :: cells(:, :), degree
        real(dp), intent(out) :: error(2)

        type(tetra_rt_t) :: basis
        integer, allocatable :: face_orientations(:, :)
        integer, allocatable :: face_permutations(:, :, :), faces(:, :)
        integer, allocatable :: global_dofs(:, :)
        real(dp), allocatable :: basis_transform(:, :), local_dofs(:)
        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: divergence, jacobian(3, 3), point(3)
        real(dp) :: reference_point(3), value(3)
        integer :: dof_count, local_status, tetrahedron

        call build_tetra_rt_dof_map( &
            degree, cells, faces, global_dofs, face_orientations, &
            face_permutations, local_status)
        if (local_status /= 0) error stop "zero-normal Hdiv map failed"
        call initialize_tetra_rt(degree, basis, local_status)
        if (local_status /= 0) error stop "zero-normal Hdiv basis failed"
        dof_count = size(global_dofs, 1)
        allocate(basis_transform(dof_count, dof_count))
        allocate(local_dofs(dof_count))
        allocate(reference_values(3, dof_count))
        allocate(reference_divergences(dof_count))
        allocate(physical_values(3, dof_count))
        allocate(physical_divergences(dof_count))
        reference_point = [0.21_dp, 0.24_dp, 0.19_dp]
        error = 0.0_dp
        do tetrahedron = 1, size(cells, 2)
            call build_tetra_rt_basis_transform( &
                degree, face_orientations(:, tetrahedron), &
                face_permutations(:, :, tetrahedron), basis_transform, &
                local_status)
            if (local_status /= 0) error stop "zero-normal transform failed"
            local_dofs = matmul( &
                basis_transform, coefficients(global_dofs(:, tetrahedron)))
            jacobian(:, 1) = mesh_vertices(:, cells(2, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 2) = mesh_vertices(:, cells(3, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            jacobian(:, 3) = mesh_vertices(:, cells(4, tetrahedron)) - &
                mesh_vertices(:, cells(1, tetrahedron))
            call evaluate_tetra_rt( &
                basis, reference_point, reference_values, &
                reference_divergences, local_status)
            if (local_status /= 0) error stop "zero-normal evaluation failed"
            call map_tetra_rt_contravariant( &
                jacobian, reference_values, reference_divergences, &
                physical_values, physical_divergences, local_status)
            if (local_status /= 0) error stop "zero-normal Piola map failed"
            point = mesh_vertices(:, cells(1, tetrahedron)) + &
                matmul(jacobian, reference_point)
            value = matmul(physical_values, local_dofs)
            divergence = dot_product(physical_divergences, local_dofs)
            error(1) = max( &
                error(1), maxval(abs(value - exact_field(point))))
            error(2) = max( &
                error(2), abs(divergence - exact_divergence(point)))
        end do
    end subroutine measure_error

end program test_tetra_rt_zero_normal_p_convergence_slow
