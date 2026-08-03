program test_tetra_rt_interpolation_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: initialize_tetra_rt, interpolate_physical_tetra_rt, &
        interpolate_sampled_physical_tetra_rt, &
        interpolate_sampled_physical_tetra_rt_jvp, &
        interpolate_sampled_physical_tetra_rt_vjp, &
        tetra_rt_interpolation_points, tetra_rt_t, &
        tetra_vector_sample_gradients_t, tetra_vector_samples_t, &
        zero_tetra_vector_samples_like
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(tetra_rt_t) :: basis
    type(tetra_vector_samples_t) :: parameter_dot, points, points_minus
    type(tetra_vector_samples_t) :: points_plus, samples, samples_bar
    type(tetra_vector_samples_t) :: samples_minus, samples_plus
    type(tetra_vector_sample_gradients_t) :: gradients
    real(dp), allocatable :: dofs(:), dofs_bar(:), dofs_callback(:)
    real(dp), allocatable :: dofs_dot(:), dofs_minus(:), dofs_plus(:)
    real(dp) :: lhs, rhs, vertices(3, 4), vertices_bar(3, 4)
    real(dp) :: vertices_dot(3, 4)
    integer :: failures, index, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, 1.3_dp, 0.1_dp, -0.1_dp, &
        -0.15_dp, 1.2_dp, 0.2_dp, 0.2_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.025_dp, 0.01_dp, &
        0.02_dp, 0.015_dp, -0.025_dp, -0.01_dp, 0.03_dp, 0.02_dp], [3, 4])
    call initialize_tetra_rt(2, basis, status)
    call tetra_rt_interpolation_points(basis, vertices, points, status)
    call tetra_rt_interpolation_points( &
        basis, vertices + step*vertices_dot, points_plus, status)
    call tetra_rt_interpolation_points( &
        basis, vertices - step*vertices_dot, points_minus, status)
    call initialize_sample_data(points, samples, gradients, parameter_dot)
    call evaluate_samples(points_plus, parameter_dot, step, samples_plus)
    call evaluate_samples(points_minus, parameter_dot, -step, samples_minus)

    call interpolate_sampled_physical_tetra_rt( &
        basis, vertices, samples, dofs, status)
    call interpolate_physical_tetra_rt( &
        basis, vertices, field_xyz, dofs_callback, status)
    call check(maxval(abs(dofs - dofs_callback)) < 2.0e-11_dp, &
        "sampled RT interpolation matches callback")
    allocate(dofs_dot(size(dofs)), dofs_plus(size(dofs)), dofs_minus(size(dofs)))
    allocate(dofs_bar(size(dofs)))
    do index = 1, size(dofs)
        dofs_bar(index) = 0.01_dp*index - 0.08_dp
    end do
    call interpolate_sampled_physical_tetra_rt_jvp( &
        basis, vertices, samples, gradients, vertices_dot, parameter_dot, &
        dofs_dot, status)
    call interpolate_sampled_physical_tetra_rt( &
        basis, vertices + step*vertices_dot, samples_plus, dofs_plus, status)
    call interpolate_sampled_physical_tetra_rt( &
        basis, vertices - step*vertices_dot, samples_minus, dofs_minus, status)
    call check(maxval(abs(dofs_dot - &
        (dofs_plus - dofs_minus)/(2.0_dp*step))) < 1.0e-7_dp, &
        "sampled RT interpolation JVP")
    call interpolate_sampled_physical_tetra_rt_vjp( &
        basis, vertices, samples, gradients, dofs_bar, vertices_bar, &
        samples_bar, status)
    lhs = dot_product(dofs_bar, dofs_dot)
    rhs = sum(vertices_bar*vertices_dot) + sample_dot(samples_bar, parameter_dot)
    call check(abs(lhs - rhs) < 7.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled RT interpolation adjoint identity")

    if (failures /= 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine initialize_sample_data( &
            sample_points, field_samples, field_gradients, tangent_samples)
        type(tetra_vector_samples_t), intent(in) :: sample_points
        type(tetra_vector_samples_t), intent(out) :: field_samples
        type(tetra_vector_sample_gradients_t), intent(out) :: field_gradients
        type(tetra_vector_samples_t), intent(out) :: tangent_samples

        call zero_tetra_vector_samples_like(sample_points, field_samples)
        call zero_tetra_vector_samples_like(sample_points, tangent_samples)
        allocate(field_gradients%edge(3, 3, 0, 0))
        allocate(field_gradients%face(3, 3, size(sample_points%face, 2), 4))
        allocate(field_gradients%volume(3, 3, size(sample_points%volume, 2)))
        call fill_values(sample_points%face, field_samples%face, &
            field_gradients%face, tangent_samples%face)
        call fill_values_volume(sample_points%volume, field_samples%volume, &
            field_gradients%volume, tangent_samples%volume)
    end subroutine initialize_sample_data

    subroutine evaluate_samples( &
            sample_points, tangent_samples, scale, field_samples)
        type(tetra_vector_samples_t), intent(in) :: sample_points
        type(tetra_vector_samples_t), intent(in) :: tangent_samples
        real(dp), intent(in) :: scale
        type(tetra_vector_samples_t), intent(out) :: field_samples
        type(tetra_vector_sample_gradients_t) :: unused_gradients
        type(tetra_vector_samples_t) :: unused_tangent

        call initialize_sample_data( &
            sample_points, field_samples, unused_gradients, unused_tangent)
        field_samples%face = field_samples%face + scale*tangent_samples%face
        field_samples%volume = &
            field_samples%volume + scale*tangent_samples%volume
    end subroutine evaluate_samples

    subroutine fill_values(points_array, values_array, gradient_array, tangent)
        real(dp), intent(in) :: points_array(:, :, :)
        real(dp), intent(out) :: values_array(:, :, :)
        real(dp), intent(out) :: gradient_array(:, :, :, :)
        real(dp), intent(out) :: tangent(:, :, :)
        integer :: node, topology

        do topology = 1, size(points_array, 3)
            do node = 1, size(points_array, 2)
                call field_at_point(points_array(:, node, topology), &
                    values_array(:, node, topology))
                call field_gradient(gradient_array(:, :, node, topology))
                tangent(:, node, topology) = &
                    [0.003_dp*node, -0.002_dp*topology, 0.001_dp*(node + topology)]
            end do
        end do
    end subroutine fill_values

    subroutine fill_values_volume(points_array, values_array, gradient_array, &
            tangent)
        real(dp), intent(in) :: points_array(:, :)
        real(dp), intent(out) :: values_array(:, :)
        real(dp), intent(out) :: gradient_array(:, :, :)
        real(dp), intent(out) :: tangent(:, :)
        integer :: node

        do node = 1, size(points_array, 2)
            call field_at_point(points_array(:, node), values_array(:, node))
            call field_gradient(gradient_array(:, :, node))
            tangent(:, node) = &
                [0.003_dp*node, -0.002_dp*node, 0.001_dp*node]
        end do
    end subroutine fill_values_volume

    pure subroutine field_at_point(point, value)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)

        value = [point(1) + point(2), point(2) + point(3), &
            point(3) + point(1)]
    end subroutine field_at_point

    pure subroutine field_xyz(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)

        call field_at_point([x, y, z], value)
    end subroutine field_xyz

    pure subroutine field_gradient(gradient)
        real(dp), intent(out) :: gradient(3, 3)

        gradient = reshape([ &
            1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 1.0_dp], [3, 3])
    end subroutine field_gradient

    pure real(dp) function sample_dot(first, second) result(value)
        type(tetra_vector_samples_t), intent(in) :: first, second

        value = sum(first%face*second%face) + &
            sum(first%volume*second%volume)
    end function sample_dot

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_rt_interpolation_ad
