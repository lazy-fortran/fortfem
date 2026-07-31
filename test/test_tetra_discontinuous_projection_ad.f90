program test_tetra_discontinuous_projection_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: &
        initialize_tetra_discontinuous, &
        project_physical_tetra_discontinuous, &
        project_sampled_physical_tetra_discontinuous, &
        project_sampled_physical_tetra_discontinuous_jvp, &
        project_sampled_physical_tetra_discontinuous_vjp, &
        tetra_discontinuous_t
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    implicit none

    integer, parameter :: degree = 2
    real(dp), parameter :: step = 2.0e-6_dp
    type(tetra_discontinuous_t) :: basis
    real(dp), allocatable :: coefficients(:), coefficients_bar(:)
    real(dp), allocatable :: coefficients_callback(:)
    real(dp), allocatable :: coefficients_dot(:), coefficients_minus(:)
    real(dp), allocatable :: coefficients_plus(:), gradients(:, :)
    real(dp), allocatable :: samples(:), samples_bar(:), samples_dot(:)
    real(dp), allocatable :: samples_minus(:), samples_plus(:), weights(:)
    real(dp), allocatable :: x(:), y(:), z(:)
    real(dp) :: jacobian(3, 3), jacobian_dot(3, 3), lhs, point(3), point_dot(3)
    real(dp) :: rhs, vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    integer :: failures, node, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, 1.3_dp, 0.1_dp, -0.1_dp, &
        -0.15_dp, 1.2_dp, 0.2_dp, 0.2_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, -0.03_dp, 0.025_dp, 0.01_dp, &
        0.02_dp, 0.015_dp, -0.025_dp, -0.01_dp, 0.03_dp, 0.02_dp], [3, 4])
    call initialize_tetra_discontinuous(degree, basis, status)
    call tetra_duffy_quadrature(2*degree + 8, x, y, z, weights, status)
    allocate(samples(size(weights)), samples_dot(size(weights)))
    allocate(samples_plus(size(weights)), samples_minus(size(weights)))
    allocate(samples_bar(size(weights)), gradients(3, size(weights)))
    call tetra_jacobian(vertices, jacobian)
    call tetra_jacobian(vertices_dot, jacobian_dot)
    do node = 1, size(weights)
        point = vertices(:, 1) + matmul(jacobian, [x(node), y(node), z(node)])
        point_dot = vertices_dot(:, 1) + &
            matmul(jacobian_dot, [x(node), y(node), z(node)])
        samples(node) = field(point)
        gradients(:, node) = [2.0_dp*point(1), 2.0_dp, -0.5_dp]
        samples_dot(node) = 0.004_dp*node - 0.03_dp
        samples_plus(node) = field(point + step*point_dot) + &
            step*samples_dot(node)
        samples_minus(node) = field(point - step*point_dot) - &
            step*samples_dot(node)
    end do

    call project_sampled_physical_tetra_discontinuous( &
        basis, vertices, samples, coefficients, status)
    call project_physical_tetra_discontinuous( &
        basis, vertices, field_xyz, coefficients_callback, status)
    call check(maxval(abs(coefficients - coefficients_callback)) < 2.0e-11_dp, &
        "sampled projection matches callback projection")
    allocate(coefficients_dot(size(coefficients)))
    allocate(coefficients_plus(size(coefficients)), source=0.0_dp)
    allocate(coefficients_minus(size(coefficients)), source=0.0_dp)
    allocate(coefficients_bar(size(coefficients)))
    coefficients_bar = [(0.02_dp*node - 0.07_dp, node = 1, size(coefficients))]
    call project_sampled_physical_tetra_discontinuous_jvp( &
        basis, vertices, samples, gradients, vertices_dot, samples_dot, &
        coefficients_dot, status)
    call project_sampled_physical_tetra_discontinuous( &
        basis, vertices + step*vertices_dot, samples_plus, coefficients_plus, &
        status)
    call project_sampled_physical_tetra_discontinuous( &
        basis, vertices - step*vertices_dot, samples_minus, &
        coefficients_minus, status)
    call check_close(maxval(abs(coefficients_dot - &
        (coefficients_plus - coefficients_minus)/(2.0_dp*step))), 2.0e-7_dp, &
        "sampled discontinuous projection JVP")
    call project_sampled_physical_tetra_discontinuous_vjp( &
        basis, vertices, samples, gradients, coefficients_bar, vertices_bar, &
        samples_bar, status)
    lhs = dot_product(coefficients_bar, coefficients_dot)
    rhs = sum(vertices_bar*vertices_dot) + dot_product(samples_bar, samples_dot)
    call check(abs(lhs - rhs) < 3.0e-11_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "sampled discontinuous projection adjoint identity")

    if (failures /= 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    pure real(dp) function field(point) result(value)
        real(dp), intent(in) :: point(3)

        value = point(1)**2 + 2.0_dp*point(2) - 0.5_dp*point(3)
    end function field

    pure real(dp) function field_xyz(x_value, y_value, z_value) result(value)
        real(dp), intent(in) :: x_value, y_value, z_value

        value = x_value**2 + 2.0_dp*y_value - 0.5_dp*z_value
    end function field_xyz

    pure subroutine tetra_jacobian(local_vertices, local_jacobian)
        real(dp), intent(in) :: local_vertices(3, 4)
        real(dp), intent(out) :: local_jacobian(3, 3)

        local_jacobian(:, 1) = local_vertices(:, 2) - local_vertices(:, 1)
        local_jacobian(:, 2) = local_vertices(:, 3) - local_vertices(:, 1)
        local_jacobian(:, 3) = local_vertices(:, 4) - local_vertices(:, 1)
    end subroutine tetra_jacobian

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

    subroutine check_close(error, tolerance, label)
        real(dp), intent(in) :: error, tolerance
        character(*), intent(in) :: label

        if (error < tolerance) return
        failures = failures + 1
        write(error_unit, "(a,es12.4)") "FAIL: "//label//" error=", error
    end subroutine check_close

end program test_tetra_discontinuous_projection_ad
