module fortfem_tetra_discontinuous_projection
    use fortfem_kinds, only: dp
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, tetra_discontinuous_dof_count, &
        tetra_discontinuous_t
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortnum_linalg, only: dense_solve, det3
    implicit none

    private

    public :: project_physical_tetra_discontinuous
    public :: project_sampled_physical_tetra_discontinuous
    public :: project_sampled_physical_tetra_discontinuous_jvp
    public :: project_sampled_physical_tetra_discontinuous_vjp

    abstract interface
        pure real(dp) function physical_scalar_field_3d(x, y, z) result(value)
            import :: dp
            real(dp), intent(in) :: x, y, z
        end function physical_scalar_field_3d
    end interface

contains

    subroutine project_physical_tetra_discontinuous( &
            basis, vertices, field, coefficients, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        procedure(physical_scalar_field_3d) :: field
        real(dp), allocatable, intent(out) :: coefficients(:)
        integer, intent(out) :: status

        real(dp), allocatable :: mass(:, :), right_hand_side(:), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        integer :: column, degree, dof_count, node, row

        status = 1
        degree = degree_from_dof_count(tetra_discontinuous_dof_count(basis))
        if (degree < 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call tetra_duffy_quadrature( &
            2*degree + 8, x, y, z, weights, status)
        if (status /= 0) return
        dof_count = tetra_discontinuous_dof_count(basis)
        allocate ( &
            mass(dof_count, dof_count), right_hand_side(dof_count), &
            values(dof_count), coefficients(dof_count))
        mass = 0.0_dp
        right_hand_side = 0.0_dp
        do node = 1, size(weights)
            call evaluate_tetra_discontinuous( &
                basis, x(node), y(node), z(node), values, status)
            if (status /= 0) return
            physical_point = vertices(:, 1) + &
                matmul(jacobian, [x(node), y(node), z(node)])
            do column = 1, dof_count
                right_hand_side(column) = right_hand_side(column) + &
                    weights(node)*values(column)*field( &
                    physical_point(1), physical_point(2), physical_point(3))
                do row = 1, dof_count
                    mass(row, column) = mass(row, column) + &
                        weights(node)*values(row)*values(column)
                end do
            end do
        end do
        call dense_solve(mass, right_hand_side, coefficients, status)
    end subroutine project_physical_tetra_discontinuous

    subroutine project_sampled_physical_tetra_discontinuous( &
            basis, vertices, field_samples, coefficients, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4), field_samples(:)
        real(dp), allocatable, intent(out) :: coefficients(:)
        integer, intent(out) :: status

        real(dp), allocatable :: mass(:, :), reference_points(:, :)
        real(dp), allocatable :: right_hand_side(:), sampling(:, :)

        call initialize_sampled_projection( &
            basis, vertices, mass, sampling, reference_points, status)
        if (status /= 0) return
        if (size(field_samples) /= size(sampling, 2)) then
            status = 1
            return
        end if
        right_hand_side = matmul(sampling, field_samples)
        allocate(coefficients(size(right_hand_side)))
        call dense_solve(mass, right_hand_side, coefficients, status)
    end subroutine project_sampled_physical_tetra_discontinuous

    subroutine project_sampled_physical_tetra_discontinuous_jvp( &
            basis, vertices, field_samples, field_spatial_gradients, &
            vertices_dot, field_parameter_dot, coefficients_dot, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        real(dp), intent(in) :: field_samples(:)
        real(dp), intent(in) :: field_spatial_gradients(:, :)
        real(dp), intent(in) :: field_parameter_dot(:)
        real(dp), intent(out) :: coefficients_dot(:)
        integer, intent(out) :: status

        real(dp), allocatable :: mass(:, :), reference_points(:, :)
        real(dp), allocatable :: right_hand_side_dot(:), sample_dot(:)
        real(dp), allocatable :: sampling(:, :)
        real(dp) :: jacobian_dot(3, 3), point_dot(3)
        integer :: node

        coefficients_dot = 0.0_dp
        call initialize_sampled_projection( &
            basis, vertices, mass, sampling, reference_points, status)
        if (status /= 0) return
        if (.not. valid_sampled_shapes( &
            field_samples, field_spatial_gradients, field_parameter_dot, &
            size(sampling, 2))) then
            status = 1
            return
        end if
        if (size(coefficients_dot) /= size(sampling, 1)) then
            status = 1
            return
        end if
        call tetrahedron_jacobian(vertices_dot, jacobian_dot)
        allocate(sample_dot(size(field_samples)))
        do node = 1, size(sample_dot)
            point_dot = vertices_dot(:, 1) + &
                matmul(jacobian_dot, reference_points(:, node))
            sample_dot(node) = field_parameter_dot(node) + &
                dot_product(field_spatial_gradients(:, node), point_dot)
        end do
        right_hand_side_dot = matmul(sampling, sample_dot)
        call dense_solve( &
            mass, right_hand_side_dot, coefficients_dot, status)
    end subroutine project_sampled_physical_tetra_discontinuous_jvp

    subroutine project_sampled_physical_tetra_discontinuous_vjp( &
            basis, vertices, field_samples, field_spatial_gradients, &
            coefficients_bar, vertices_bar, field_samples_bar, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4), field_samples(:)
        real(dp), intent(in) :: field_spatial_gradients(:, :)
        real(dp), intent(in) :: coefficients_bar(:)
        real(dp), intent(out) :: vertices_bar(3, 4), field_samples_bar(:)
        integer, intent(out) :: status

        real(dp), allocatable :: dual(:), mass(:, :)
        real(dp), allocatable :: reference_points(:, :), sampling(:, :)
        real(dp) :: barycentric(4), point_bar(3)
        integer :: node

        vertices_bar = 0.0_dp
        field_samples_bar = 0.0_dp
        call initialize_sampled_projection( &
            basis, vertices, mass, sampling, reference_points, status)
        if (status /= 0) return
        if (size(field_samples) /= size(sampling, 2)) then
            status = 1
            return
        end if
        if (size(field_spatial_gradients, 1) /= 3 .or. &
            size(field_spatial_gradients, 2) /= size(sampling, 2)) then
            status = 1
            return
        end if
        if (size(coefficients_bar) /= size(sampling, 1)) then
            status = 1
            return
        end if
        if (size(field_samples_bar) /= size(sampling, 2)) then
            status = 1
            return
        end if
        allocate(dual(size(coefficients_bar)))
        call dense_solve( &
            transpose(mass), coefficients_bar, dual, status)
        if (status /= 0) return
        field_samples_bar = matmul(transpose(sampling), dual)
        do node = 1, size(field_samples_bar)
            point_bar = field_samples_bar(node)* &
                field_spatial_gradients(:, node)
            barycentric = [ &
                1.0_dp - sum(reference_points(:, node)), &
                reference_points(:, node)]
            vertices_bar = vertices_bar + &
                spread(point_bar, 2, 4)*spread(barycentric, 1, 3)
        end do
    end subroutine project_sampled_physical_tetra_discontinuous_vjp

    subroutine initialize_sampled_projection( &
            basis, vertices, mass, sampling, reference_points, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), allocatable, intent(out) :: mass(:, :), sampling(:, :)
        real(dp), allocatable, intent(out) :: reference_points(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: values(:), weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3)
        integer :: degree, dof_count, node

        status = 1
        degree = degree_from_dof_count(tetra_discontinuous_dof_count(basis))
        if (degree < 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        determinant = det3(jacobian)
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call tetra_duffy_quadrature( &
            2*degree + 8, x, y, z, weights, status)
        if (status /= 0) return
        dof_count = tetra_discontinuous_dof_count(basis)
        allocate(mass(dof_count, dof_count), sampling(dof_count, size(weights)))
        allocate(reference_points(3, size(weights)), values(dof_count))
        mass = 0.0_dp
        do node = 1, size(weights)
            reference_points(:, node) = [x(node), y(node), z(node)]
            call evaluate_tetra_discontinuous( &
                basis, x(node), y(node), z(node), values, status)
            if (status /= 0) return
            sampling(:, node) = weights(node)*values
            mass = mass + weights(node)* &
                spread(values, 2, dof_count)*spread(values, 1, dof_count)
        end do
        status = 0
    end subroutine initialize_sampled_projection

    pure logical function valid_sampled_shapes( &
            samples, gradients, parameter_dot, count) result(valid)
        real(dp), intent(in) :: samples(:), gradients(:, :), parameter_dot(:)
        integer, intent(in) :: count

        valid = size(samples) == count
        if (.not. valid) return
        valid = size(parameter_dot) == count
        if (.not. valid) return
        valid = size(gradients, 1) == 3
        if (.not. valid) return
        valid = size(gradients, 2) == count
    end function valid_sampled_shapes

    pure subroutine tetrahedron_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetrahedron_jacobian

    pure integer function degree_from_dof_count(dof_count) result(degree)
        integer, intent(in) :: dof_count

        do degree = 0, 1000
            if (dof_count == (degree + 1)*(degree + 2)*(degree + 3)/6) return
            if ((degree + 1)*(degree + 2)*(degree + 3)/6 > dof_count) exit
        end do
        degree = -1
    end function degree_from_dof_count

end module fortfem_tetra_discontinuous_projection
