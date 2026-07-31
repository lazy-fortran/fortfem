module fortfem_tetra_vector_evaluation
    !! Differentiable fixed-reference observations of tetrahedral vector fields.
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        evaluate_tetra_nedelec_first_kind, tetra_nedelec_first_kind_t
    use fortfem_tetra_piola_maps, only: &
        map_tetra_nedelec_covariant, map_tetra_nedelec_covariant_jvp, &
        map_tetra_nedelec_covariant_vjp, map_tetra_rt_contravariant, &
        map_tetra_rt_contravariant_jvp, map_tetra_rt_contravariant_vjp
    use fortfem_tetra_rt_arbitrary_order, only: evaluate_tetra_rt, tetra_rt_t
    implicit none

    private

    public :: evaluate_tetra_nedelec_interpolant
    public :: evaluate_tetra_nedelec_interpolant_jvp
    public :: evaluate_tetra_nedelec_interpolant_vjp
    public :: evaluate_tetra_rt_interpolant
    public :: evaluate_tetra_rt_interpolant_jvp
    public :: evaluate_tetra_rt_interpolant_vjp

contains

    subroutine evaluate_tetra_nedelec_interpolant( &
            vertices, basis, dofs, reference_point, value, curl, status)
        real(dp), intent(in) :: vertices(3, 4), dofs(:), reference_point(3)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(out) :: value(3), curl(3)
        integer, intent(out) :: status

        real(dp), allocatable :: physical_curls(:, :), physical_values(:, :)
        real(dp), allocatable :: reference_curls(:, :), reference_values(:, :)
        real(dp) :: jacobian(3, 3)

        value = 0.0_dp
        curl = 0.0_dp
        status = 1
        if (size(dofs) < 1) return
        allocate(reference_values(3, size(dofs)))
        allocate(reference_curls(3, size(dofs)))
        allocate(physical_values(3, size(dofs)))
        allocate(physical_curls(3, size(dofs)))
        call evaluate_tetra_nedelec_first_kind( &
            basis, reference_point, reference_values, reference_curls, status)
        if (status /= 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        call map_tetra_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        if (status /= 0) return
        value = matmul(physical_values, dofs)
        curl = matmul(physical_curls, dofs)
    end subroutine evaluate_tetra_nedelec_interpolant

    subroutine evaluate_tetra_nedelec_interpolant_jvp( &
            vertices, basis, dofs, reference_point, vertices_dot, dofs_dot, &
            value_dot, curl_dot, status)
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: dofs(:), dofs_dot(:), reference_point(3)
        real(dp), intent(out) :: value_dot(3), curl_dot(3)
        integer, intent(out) :: status

        real(dp), allocatable :: physical_curls(:, :)
        real(dp), allocatable :: physical_curls_dot(:, :)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: physical_values_dot(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_curls_dot(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: reference_values_dot(:, :)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3)

        value_dot = 0.0_dp
        curl_dot = 0.0_dp
        status = 1
        if (size(dofs) < 1) return
        if (size(dofs_dot) /= size(dofs)) return
        allocate(reference_values(3, size(dofs)))
        allocate(reference_values_dot(3, size(dofs)), source=0.0_dp)
        allocate(reference_curls(3, size(dofs)))
        allocate(reference_curls_dot(3, size(dofs)), source=0.0_dp)
        allocate(physical_values(3, size(dofs)))
        allocate(physical_values_dot(3, size(dofs)))
        allocate(physical_curls(3, size(dofs)))
        allocate(physical_curls_dot(3, size(dofs)))
        call evaluate_tetra_nedelec_first_kind( &
            basis, reference_point, reference_values, reference_curls, status)
        if (status /= 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        call tetrahedron_jacobian(vertices_dot, jacobian_dot)
        call map_tetra_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        if (status /= 0) return
        call map_tetra_nedelec_covariant_jvp( &
            jacobian, reference_values, reference_curls, jacobian_dot, &
            reference_values_dot, reference_curls_dot, physical_values_dot, &
            physical_curls_dot, status)
        if (status /= 0) return
        value_dot = matmul(physical_values_dot, dofs) + &
            matmul(physical_values, dofs_dot)
        curl_dot = matmul(physical_curls_dot, dofs) + &
            matmul(physical_curls, dofs_dot)
    end subroutine evaluate_tetra_nedelec_interpolant_jvp

    subroutine evaluate_tetra_nedelec_interpolant_vjp( &
            vertices, basis, dofs, reference_point, value_bar, curl_bar, &
            vertices_bar, dofs_bar, status)
        real(dp), intent(in) :: vertices(3, 4), dofs(:), reference_point(3)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: value_bar(3), curl_bar(3)
        real(dp), intent(out) :: vertices_bar(3, 4), dofs_bar(:)
        integer, intent(out) :: status

        real(dp), allocatable :: physical_curls(:, :)
        real(dp), allocatable :: physical_curls_bar(:, :)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: physical_values_bar(:, :)
        real(dp), allocatable :: reference_curls(:, :)
        real(dp), allocatable :: reference_curls_bar(:, :)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: reference_values_bar(:, :)
        real(dp) :: jacobian(3, 3), jacobian_bar(3, 3)
        integer :: dof

        vertices_bar = 0.0_dp
        dofs_bar = 0.0_dp
        status = 1
        if (size(dofs) < 1) return
        if (size(dofs_bar) /= size(dofs)) return
        allocate(reference_values(3, size(dofs)))
        allocate(reference_values_bar(3, size(dofs)))
        allocate(reference_curls(3, size(dofs)))
        allocate(reference_curls_bar(3, size(dofs)))
        allocate(physical_values(3, size(dofs)))
        allocate(physical_values_bar(3, size(dofs)))
        allocate(physical_curls(3, size(dofs)))
        allocate(physical_curls_bar(3, size(dofs)))
        call evaluate_tetra_nedelec_first_kind( &
            basis, reference_point, reference_values, reference_curls, status)
        if (status /= 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        call map_tetra_nedelec_covariant( &
            jacobian, reference_values, reference_curls, physical_values, &
            physical_curls, status)
        if (status /= 0) return
        dofs_bar = matmul(transpose(physical_values), value_bar) + &
            matmul(transpose(physical_curls), curl_bar)
        do dof = 1, size(dofs)
            physical_values_bar(:, dof) = value_bar*dofs(dof)
            physical_curls_bar(:, dof) = curl_bar*dofs(dof)
        end do
        call map_tetra_nedelec_covariant_vjp( &
            jacobian, reference_values, reference_curls, physical_values_bar, &
            physical_curls_bar, jacobian_bar, reference_values_bar, &
            reference_curls_bar, status)
        if (status /= 0) return
        call scatter_jacobian_bar(jacobian_bar, vertices_bar)
    end subroutine evaluate_tetra_nedelec_interpolant_vjp

    subroutine evaluate_tetra_rt_interpolant( &
            vertices, basis, dofs, reference_point, value, divergence, status)
        real(dp), intent(in) :: vertices(3, 4), dofs(:), reference_point(3)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(out) :: value(3), divergence
        integer, intent(out) :: status

        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp) :: jacobian(3, 3)

        value = 0.0_dp
        divergence = 0.0_dp
        status = 1
        if (size(dofs) < 1) return
        allocate(reference_values(3, size(dofs)))
        allocate(reference_divergences(size(dofs)))
        allocate(physical_values(3, size(dofs)))
        allocate(physical_divergences(size(dofs)))
        call evaluate_tetra_rt( &
            basis, reference_point, reference_values, reference_divergences, &
            status)
        if (status /= 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        call map_tetra_rt_contravariant( &
            jacobian, reference_values, reference_divergences, &
            physical_values, physical_divergences, status)
        if (status /= 0) return
        value = matmul(physical_values, dofs)
        divergence = dot_product(physical_divergences, dofs)
    end subroutine evaluate_tetra_rt_interpolant

    subroutine evaluate_tetra_rt_interpolant_jvp( &
            vertices, basis, dofs, reference_point, vertices_dot, dofs_dot, &
            value_dot, divergence_dot, status)
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: dofs(:), dofs_dot(:), reference_point(3)
        real(dp), intent(out) :: value_dot(3), divergence_dot
        integer, intent(out) :: status

        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_divergences_dot(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: physical_values_dot(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_divergences_dot(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: reference_values_dot(:, :)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3)

        value_dot = 0.0_dp
        divergence_dot = 0.0_dp
        status = 1
        if (size(dofs) < 1) return
        if (size(dofs_dot) /= size(dofs)) return
        allocate(reference_values(3, size(dofs)))
        allocate(reference_values_dot(3, size(dofs)), source=0.0_dp)
        allocate(reference_divergences(size(dofs)))
        allocate(reference_divergences_dot(size(dofs)), source=0.0_dp)
        allocate(physical_values(3, size(dofs)))
        allocate(physical_values_dot(3, size(dofs)))
        allocate(physical_divergences(size(dofs)))
        allocate(physical_divergences_dot(size(dofs)))
        call evaluate_tetra_rt( &
            basis, reference_point, reference_values, reference_divergences, &
            status)
        if (status /= 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        call tetrahedron_jacobian(vertices_dot, jacobian_dot)
        call map_tetra_rt_contravariant( &
            jacobian, reference_values, reference_divergences, &
            physical_values, physical_divergences, status)
        if (status /= 0) return
        call map_tetra_rt_contravariant_jvp( &
            jacobian, reference_values, reference_divergences, jacobian_dot, &
            reference_values_dot, reference_divergences_dot, &
            physical_values_dot, physical_divergences_dot, status)
        if (status /= 0) return
        value_dot = matmul(physical_values_dot, dofs) + &
            matmul(physical_values, dofs_dot)
        divergence_dot = dot_product(physical_divergences_dot, dofs) + &
            dot_product(physical_divergences, dofs_dot)
    end subroutine evaluate_tetra_rt_interpolant_jvp

    subroutine evaluate_tetra_rt_interpolant_vjp( &
            vertices, basis, dofs, reference_point, value_bar, divergence_bar, &
            vertices_bar, dofs_bar, status)
        real(dp), intent(in) :: vertices(3, 4), dofs(:), reference_point(3)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: value_bar(3), divergence_bar
        real(dp), intent(out) :: vertices_bar(3, 4), dofs_bar(:)
        integer, intent(out) :: status

        real(dp), allocatable :: physical_divergences(:)
        real(dp), allocatable :: physical_divergences_bar(:)
        real(dp), allocatable :: physical_values(:, :)
        real(dp), allocatable :: physical_values_bar(:, :)
        real(dp), allocatable :: reference_divergences(:)
        real(dp), allocatable :: reference_divergences_bar(:)
        real(dp), allocatable :: reference_values(:, :)
        real(dp), allocatable :: reference_values_bar(:, :)
        real(dp) :: jacobian(3, 3), jacobian_bar(3, 3)
        integer :: dof

        vertices_bar = 0.0_dp
        dofs_bar = 0.0_dp
        status = 1
        if (size(dofs) < 1) return
        if (size(dofs_bar) /= size(dofs)) return
        allocate(reference_values(3, size(dofs)))
        allocate(reference_values_bar(3, size(dofs)))
        allocate(reference_divergences(size(dofs)))
        allocate(reference_divergences_bar(size(dofs)))
        allocate(physical_values(3, size(dofs)))
        allocate(physical_values_bar(3, size(dofs)))
        allocate(physical_divergences(size(dofs)))
        allocate(physical_divergences_bar(size(dofs)))
        call evaluate_tetra_rt( &
            basis, reference_point, reference_values, reference_divergences, &
            status)
        if (status /= 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        call map_tetra_rt_contravariant( &
            jacobian, reference_values, reference_divergences, &
            physical_values, physical_divergences, status)
        if (status /= 0) return
        dofs_bar = matmul(transpose(physical_values), value_bar) + &
            divergence_bar*physical_divergences
        do dof = 1, size(dofs)
            physical_values_bar(:, dof) = value_bar*dofs(dof)
            physical_divergences_bar(dof) = divergence_bar*dofs(dof)
        end do
        call map_tetra_rt_contravariant_vjp( &
            jacobian, reference_values, reference_divergences, &
            physical_values_bar, physical_divergences_bar, jacobian_bar, &
            reference_values_bar, reference_divergences_bar, status)
        if (status /= 0) return
        call scatter_jacobian_bar(jacobian_bar, vertices_bar)
    end subroutine evaluate_tetra_rt_interpolant_vjp

    pure subroutine tetrahedron_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetrahedron_jacobian

    pure subroutine scatter_jacobian_bar(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(3, 3)
        real(dp), intent(out) :: vertices_bar(3, 4)

        vertices_bar(:, 1) = &
            -jacobian_bar(:, 1) - jacobian_bar(:, 2) - jacobian_bar(:, 3)
        vertices_bar(:, 2) = jacobian_bar(:, 1)
        vertices_bar(:, 3) = jacobian_bar(:, 2)
        vertices_bar(:, 4) = jacobian_bar(:, 3)
    end subroutine scatter_jacobian_bar

end module fortfem_tetra_vector_evaluation
