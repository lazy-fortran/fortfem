module fortfem_tetra_rt_interpolation
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_rt_arbitrary_order, only: &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_tetra_vector_samples, only: &
        tetra_vector_sample_gradients_t, tetra_vector_samples_t, &
        zero_tetra_vector_samples_like
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: det3, det3_jvp, det3_vjp, inv3, inv3_jvp, inv3_vjp
    use fortnum_special_jacobi, only: tetrahedron_koornwinder, &
        triangle_dubiner
    implicit none

    private

    public :: interpolate_reference_tetra_rt
    public :: interpolate_physical_tetra_rt
    public :: interpolate_sampled_physical_tetra_rt
    public :: interpolate_sampled_physical_tetra_rt_jvp
    public :: interpolate_sampled_physical_tetra_rt_vjp
    public :: tetra_rt_interpolation_points

    abstract interface
        pure subroutine physical_vector_field_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine physical_vector_field_3d

        subroutine reference_vector_field_3d(point, value)
            import :: dp
            real(dp), intent(in) :: point(3)
            real(dp), intent(out) :: value(3)
        end subroutine reference_vector_field_3d
    end interface

contains

    subroutine interpolate_physical_tetra_rt( &
            basis, vertices, field, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        procedure(physical_vector_field_3d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        integer :: inverse_status

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        status = 1
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return
        call interpolate_reference_tetra_rt( &
            basis, pulled_back_field, dofs, status)

    contains

        subroutine pulled_back_field(reference_point, value)
            real(dp), intent(in) :: reference_point(3)
            real(dp), intent(out) :: value(3)

            real(dp) :: physical_point(3), physical_value(3)

            physical_point = vertices(:, 1) + &
                matmul(jacobian, reference_point)
            call field( &
                physical_point(1), physical_point(2), physical_point(3), &
                physical_value)
            value = determinant*matmul(inverse_jacobian, physical_value)
        end subroutine pulled_back_field

    end subroutine interpolate_physical_tetra_rt

    subroutine interpolate_reference_tetra_rt(basis, field, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        procedure(reference_vector_field_3d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: area_normal(3), field_value(3), point(3)
        integer :: component, degree, face, moment, node, total
        integer :: x_degree, y_degree, z_degree

        status = 1
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        if (degree < 0) return
        allocate(dofs(tetra_rt_dof_count(basis)))
        dofs = 0.0_dp
        moment = 0
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do total = 0, degree
                do x_degree = 0, total
                    y_degree = total - x_degree
                    moment = moment + 1
                    do node = 1, size(weights)
                        call reference_face( &
                            face, x(node), y(node), point, area_normal)
                        call field(point, field_value)
                        dofs(moment) = dofs(moment) + weights(node)* &
                            dot_product(field_value, area_normal)* &
                            face_moment_value( &
                            degree, x_degree, y_degree, &
                            x(node), y(node))
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total = 0, degree - 1
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            call field(point, field_value)
                            dofs(moment) = dofs(moment) + weights(node)* &
                                field_value(component)*volume_moment_value( &
                                degree, x_degree, y_degree, z_degree, &
                                point)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) return
        status = 0
    end subroutine interpolate_reference_tetra_rt

    subroutine tetra_rt_interpolation_points(basis, vertices, points, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        type(tetra_vector_samples_t), intent(out) :: points
        integer, intent(out) :: status

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: area_normal(3), jacobian(3, 3), reference_point(3)
        integer :: degree, face, node

        status = 1
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        if (degree < 0) return
        call tetrahedron_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) return
        allocate(points%edge(3, 0, 0))
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        allocate(points%face(3, size(weights), 4))
        do face = 1, 4
            do node = 1, size(weights)
                call reference_face( &
                    face, x(node), y(node), reference_point, area_normal)
                points%face(:, node, face) = vertices(:, 1) + &
                    matmul(jacobian, reference_point)
            end do
        end do
        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        allocate(points%volume(3, size(weights)))
        do node = 1, size(weights)
            reference_point = [x(node), y(node), z(node)]
            points%volume(:, node) = vertices(:, 1) + &
                matmul(jacobian, reference_point)
        end do
        status = 0
    end subroutine tetra_rt_interpolation_points

    subroutine interpolate_sampled_physical_tetra_rt( &
            basis, vertices, samples, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        type(tetra_vector_samples_t) :: points, pulled
        real(dp) :: determinant, inverse(3, 3), jacobian(3, 3)
        integer :: face, node

        call tetra_rt_interpolation_points(basis, vertices, points, status)
        if (status /= 0) return
        if (.not. same_sample_shape(samples, points)) then
            status = 1
            return
        end if
        call tetrahedron_jacobian(vertices, jacobian)
        determinant = det3(jacobian)
        call inv3(jacobian, inverse, status)
        if (status /= 0) return
        call zero_tetra_vector_samples_like(samples, pulled)
        do face = 1, 4
            do node = 1, size(samples%face, 2)
                pulled%face(:, node, face) = determinant* &
                    matmul(inverse, samples%face(:, node, face))
            end do
        end do
        do node = 1, size(samples%volume, 2)
            pulled%volume(:, node) = determinant* &
                matmul(inverse, samples%volume(:, node))
        end do
        call interpolate_reference_tetra_rt_samples( &
            basis, pulled, dofs, status)
    end subroutine interpolate_sampled_physical_tetra_rt

    subroutine interpolate_sampled_physical_tetra_rt_jvp( &
            basis, vertices, samples, gradients, vertices_dot, parameter_dot, &
            dofs_dot, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4), vertices_dot(3, 4)
        type(tetra_vector_samples_t), intent(in) :: samples, parameter_dot
        type(tetra_vector_sample_gradients_t), intent(in) :: gradients
        real(dp), intent(out) :: dofs_dot(:)
        integer, intent(out) :: status

        type(tetra_vector_samples_t) :: points, pulled_dot
        real(dp) :: determinant, determinant_dot, field_dot(3)
        real(dp) :: inverse(3, 3), inverse_dot(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3)
        real(dp) :: point_dot(3), reference_point(3), area_normal(3)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        integer :: degree, face, node

        dofs_dot = 0.0_dp
        call tetra_rt_interpolation_points(basis, vertices, points, status)
        if (status /= 0) return
        if (.not. valid_product_shapes( &
            samples, gradients, parameter_dot, points)) then
            status = 1
            return
        end if
        if (size(dofs_dot) /= tetra_rt_dof_count(basis)) then
            status = 1
            return
        end if
        call tetrahedron_jacobian(vertices, jacobian)
        call tetrahedron_jacobian(vertices_dot, jacobian_dot)
        call inv3_jvp(jacobian, jacobian_dot, inverse, inverse_dot, status)
        if (status /= 0) return
        call det3_jvp(jacobian, jacobian_dot, determinant_dot)
        determinant = det3(jacobian)
        call zero_tetra_vector_samples_like(samples, pulled_dot)
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do node = 1, size(weights)
                call reference_face( &
                    face, x(node), y(node), reference_point, area_normal)
                point_dot = vertices_dot(:, 1) + &
                    matmul(jacobian_dot, reference_point)
                field_dot = parameter_dot%face(:, node, face) + &
                    matmul(gradients%face(:, :, node, face), point_dot)
                pulled_dot%face(:, node, face) = &
                    determinant_dot*matmul( &
                    inverse, samples%face(:, node, face)) + &
                    determinant*matmul( &
                    inverse_dot, samples%face(:, node, face)) + &
                    determinant*matmul(inverse, field_dot)
            end do
        end do
        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do node = 1, size(weights)
            reference_point = [x(node), y(node), z(node)]
            point_dot = vertices_dot(:, 1) + &
                matmul(jacobian_dot, reference_point)
            field_dot = parameter_dot%volume(:, node) + &
                matmul(gradients%volume(:, :, node), point_dot)
            pulled_dot%volume(:, node) = determinant_dot* &
                matmul(inverse, samples%volume(:, node)) + determinant* &
                matmul(inverse_dot, samples%volume(:, node)) + determinant* &
                matmul(inverse, field_dot)
        end do
        call interpolate_reference_tetra_rt_samples_into( &
            basis, pulled_dot, dofs_dot, status)
    end subroutine interpolate_sampled_physical_tetra_rt_jvp

    subroutine interpolate_sampled_physical_tetra_rt_vjp( &
            basis, vertices, samples, gradients, dofs_bar, vertices_bar, &
            samples_bar, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        type(tetra_vector_samples_t), intent(in) :: samples
        type(tetra_vector_sample_gradients_t), intent(in) :: gradients
        real(dp), intent(in) :: dofs_bar(:)
        real(dp), intent(out) :: vertices_bar(3, 4)
        type(tetra_vector_samples_t), intent(out) :: samples_bar
        integer, intent(out) :: status

        type(tetra_vector_samples_t) :: parameter_zero, points, pulled_bar
        real(dp) :: determinant, determinant_bar, field_bar(3)
        real(dp) :: inverse(3, 3), inverse_bar(3, 3)
        real(dp) :: determinant_jacobian_bar(3, 3)
        real(dp) :: inverse_jacobian_bar(3, 3), jacobian(3, 3)
        real(dp) :: jacobian_bar(3, 3), point_bar(3), reference_point(3)
        real(dp) :: area_normal(3)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        integer :: degree, face, node

        vertices_bar = 0.0_dp
        call tetra_rt_interpolation_points(basis, vertices, points, status)
        if (status /= 0) return
        call zero_tetra_vector_samples_like(samples, parameter_zero)
        if (.not. valid_product_shapes( &
            samples, gradients, parameter_zero, points)) then
            status = 1
            return
        end if
        if (size(dofs_bar) /= tetra_rt_dof_count(basis)) then
            status = 1
            return
        end if
        call interpolate_reference_tetra_rt_samples_vjp( &
            basis, samples, dofs_bar, pulled_bar, status)
        if (status /= 0) return
        call zero_tetra_vector_samples_like(samples, samples_bar)
        call tetrahedron_jacobian(vertices, jacobian)
        determinant = det3(jacobian)
        call inv3(jacobian, inverse, status)
        if (status /= 0) return
        determinant_bar = 0.0_dp
        inverse_bar = 0.0_dp
        jacobian_bar = 0.0_dp
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do node = 1, size(weights)
                call reference_face( &
                    face, x(node), y(node), reference_point, area_normal)
                call reverse_rt_sample( &
                    determinant, inverse, reference_point, &
                    samples%face(:, node, face), &
                    gradients%face(:, :, node, face), &
                    pulled_bar%face(:, node, face), field_bar, point_bar, &
                    determinant_bar, inverse_bar, jacobian_bar)
                samples_bar%face(:, node, face) = field_bar
                vertices_bar(:, 1) = vertices_bar(:, 1) + point_bar
            end do
        end do
        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do node = 1, size(weights)
            reference_point = [x(node), y(node), z(node)]
            call reverse_rt_sample( &
                determinant, inverse, reference_point, &
                samples%volume(:, node), gradients%volume(:, :, node), &
                pulled_bar%volume(:, node), field_bar, point_bar, &
                determinant_bar, inverse_bar, jacobian_bar)
            samples_bar%volume(:, node) = field_bar
            vertices_bar(:, 1) = vertices_bar(:, 1) + point_bar
        end do
        call inv3_vjp( &
            jacobian, inverse_bar, inverse, inverse_jacobian_bar, status)
        if (status /= 0) return
        call det3_vjp( &
            jacobian, determinant_bar, determinant_jacobian_bar)
        jacobian_bar = jacobian_bar + inverse_jacobian_bar + &
            determinant_jacobian_bar
        call scatter_jacobian_bar(jacobian_bar, vertices_bar)
        status = 0
    end subroutine interpolate_sampled_physical_tetra_rt_vjp

    subroutine interpolate_reference_tetra_rt_samples( &
            basis, samples, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        allocate(dofs(tetra_rt_dof_count(basis)))
        call interpolate_reference_tetra_rt_samples_into( &
            basis, samples, dofs, status)
    end subroutine interpolate_reference_tetra_rt_samples

    subroutine interpolate_reference_tetra_rt_samples_into( &
            basis, samples, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: area_normal(3), point(3)
        integer :: component, degree, face, moment, node, total
        integer :: x_degree, y_degree, z_degree

        dofs = 0.0_dp
        status = 1
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        if (degree < 0) return
        if (size(dofs) /= tetra_rt_dof_count(basis)) return
        moment = 0
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do total = 0, degree
                do x_degree = 0, total
                    y_degree = total - x_degree
                    moment = moment + 1
                    do node = 1, size(weights)
                        call reference_face( &
                            face, x(node), y(node), point, area_normal)
                        dofs(moment) = dofs(moment) + weights(node)* &
                            dot_product( &
                            samples%face(:, node, face), area_normal)* &
                            face_moment_value( &
                            degree, x_degree, y_degree, x(node), y(node))
                    end do
                end do
            end do
        end do
        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total = 0, degree - 1
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            dofs(moment) = dofs(moment) + weights(node)* &
                                samples%volume(component, node)* &
                                volume_moment_value( &
                                degree, x_degree, y_degree, z_degree, point)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) return
        status = 0
    end subroutine interpolate_reference_tetra_rt_samples_into

    subroutine interpolate_reference_tetra_rt_samples_vjp( &
            basis, samples, dofs_bar, samples_bar, status)
        type(tetra_rt_t), intent(in) :: basis
        type(tetra_vector_samples_t), intent(in) :: samples
        real(dp), intent(in) :: dofs_bar(:)
        type(tetra_vector_samples_t), intent(out) :: samples_bar
        integer, intent(out) :: status

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: area_normal(3), point(3), seed
        integer :: component, degree, face, moment, node, total
        integer :: x_degree, y_degree, z_degree

        call zero_tetra_vector_samples_like(samples, samples_bar)
        status = 1
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        if (degree < 0) return
        if (size(dofs_bar) /= tetra_rt_dof_count(basis)) return
        moment = 0
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do total = 0, degree
                do x_degree = 0, total
                    y_degree = total - x_degree
                    moment = moment + 1
                    do node = 1, size(weights)
                        call reference_face( &
                            face, x(node), y(node), point, area_normal)
                        seed = dofs_bar(moment)*weights(node)* &
                            face_moment_value( &
                            degree, x_degree, y_degree, x(node), y(node))
                        samples_bar%face(:, node, face) = &
                            samples_bar%face(:, node, face) + seed*area_normal
                    end do
                end do
            end do
        end do
        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total = 0, degree - 1
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            seed = dofs_bar(moment)*weights(node)* &
                                volume_moment_value( &
                                degree, x_degree, y_degree, z_degree, point)
                            samples_bar%volume(component, node) = &
                                samples_bar%volume(component, node) + seed
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs_bar)) return
        status = 0
    end subroutine interpolate_reference_tetra_rt_samples_vjp

    pure subroutine reverse_rt_sample( &
            determinant, inverse, reference_point, sample, gradient, &
            pulled_bar, sample_bar, point_bar, determinant_bar, inverse_bar, &
            jacobian_bar)
        real(dp), intent(in) :: determinant, inverse(3, 3), reference_point(3)
        real(dp), intent(in) :: sample(3), gradient(3, 3), pulled_bar(3)
        real(dp), intent(out) :: sample_bar(3), point_bar(3)
        real(dp), intent(inout) :: determinant_bar, inverse_bar(3, 3)
        real(dp), intent(inout) :: jacobian_bar(3, 3)

        sample_bar = determinant*matmul(transpose(inverse), pulled_bar)
        point_bar = matmul(transpose(gradient), sample_bar)
        determinant_bar = determinant_bar + &
            dot_product(pulled_bar, matmul(inverse, sample))
        inverse_bar = inverse_bar + determinant* &
            spread(pulled_bar, 2, 3)*spread(sample, 1, 3)
        jacobian_bar = jacobian_bar + &
            spread(point_bar, 2, 3)*spread(reference_point, 1, 3)
    end subroutine reverse_rt_sample

    pure subroutine tetrahedron_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetrahedron_jacobian

    pure subroutine scatter_jacobian_bar(jacobian_bar, vertices_bar)
        real(dp), intent(in) :: jacobian_bar(3, 3)
        real(dp), intent(inout) :: vertices_bar(3, 4)

        vertices_bar(:, 1) = vertices_bar(:, 1) - jacobian_bar(:, 1) - &
            jacobian_bar(:, 2) - jacobian_bar(:, 3)
        vertices_bar(:, 2) = vertices_bar(:, 2) + jacobian_bar(:, 1)
        vertices_bar(:, 3) = vertices_bar(:, 3) + jacobian_bar(:, 2)
        vertices_bar(:, 4) = vertices_bar(:, 4) + jacobian_bar(:, 3)
    end subroutine scatter_jacobian_bar

    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp) :: determinant

        determinant = det3(jacobian)
        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
    end function valid_jacobian

    pure logical function same_sample_shape(first, second) result(valid)
        type(tetra_vector_samples_t), intent(in) :: first, second

        valid = allocated(first%edge)
        if (.not. valid) return
        valid = allocated(first%face)
        if (.not. valid) return
        valid = allocated(first%volume)
        if (.not. valid) return
        valid = allocated(second%edge)
        if (.not. valid) return
        valid = allocated(second%face)
        if (.not. valid) return
        valid = allocated(second%volume)
        if (.not. valid) return
        valid = all(shape(first%edge) == shape(second%edge))
        if (.not. valid) return
        valid = all(shape(first%face) == shape(second%face))
        if (.not. valid) return
        valid = all(shape(first%volume) == shape(second%volume))
    end function same_sample_shape

    pure logical function valid_product_shapes( &
            samples, gradients, parameter_dot, points) result(valid)
        type(tetra_vector_samples_t), intent(in) :: samples, parameter_dot
        type(tetra_vector_samples_t), intent(in) :: points
        type(tetra_vector_sample_gradients_t), intent(in) :: gradients

        valid = same_sample_shape(samples, points)
        if (.not. valid) return
        valid = same_sample_shape(parameter_dot, points)
        if (.not. valid) return
        valid = allocated(gradients%edge)
        if (.not. valid) return
        valid = allocated(gradients%face)
        if (.not. valid) return
        valid = allocated(gradients%volume)
        if (.not. valid) return
        valid = all(shape(gradients%edge) == [3, 3, 0, 0])
        if (.not. valid) return
        valid = all(shape(gradients%face) == [ &
            3, 3, size(points%face, 2), size(points%face, 3)])
        if (.not. valid) return
        valid = all(shape(gradients%volume) == [3, 3, size(points%volume, 2)])
    end function valid_product_shapes

    pure function degree_from_dof_count(dof_count) result(degree)
        integer, intent(in) :: dof_count
        integer :: degree

        do degree = 0, 1000
            if (dof_count == &
                (degree + 1)*(degree + 2)*(degree + 4)/2) return
            if ((degree + 1)*(degree + 2)*(degree + 4)/2 > dof_count) exit
        end do
        degree = -1
    end function degree_from_dof_count

    pure real(dp) function face_moment_value( &
            degree, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: degree, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (degree <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function face_moment_value

    pure real(dp) function volume_moment_value( &
            degree, first_degree, second_degree, third_degree, point) &
            result(value)
        integer, intent(in) :: degree
        integer, intent(in) :: first_degree, second_degree, third_degree
        real(dp), intent(in) :: point(3)

        if (degree <= 5) then
            value = point(1)**first_degree*point(2)**second_degree* &
                point(3)**third_degree
        else
            value = tetrahedron_koornwinder( &
                first_degree, second_degree, third_degree, &
                point(1), point(2), point(3))
        end if
    end function volume_moment_value

    pure subroutine reference_face(face, u, v, point, area_normal)
        integer, intent(in) :: face
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: point(3), area_normal(3)

        select case (face)
        case (1)
            point = [0.0_dp, u, v]
            area_normal = [-1.0_dp, 0.0_dp, 0.0_dp]
        case (2)
            point = [u, 0.0_dp, v]
            area_normal = [0.0_dp, -1.0_dp, 0.0_dp]
        case (3)
            point = [v, u, 0.0_dp]
            area_normal = [0.0_dp, 0.0_dp, -1.0_dp]
        case default
            point = [1.0_dp - u - v, u, v]
            area_normal = [1.0_dp, 1.0_dp, 1.0_dp]
        end select
    end subroutine reference_face

end module fortfem_tetra_rt_interpolation
