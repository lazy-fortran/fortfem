program test_helmholtz_bem_assembly_ad
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_helmholtz_single_layer_p0_3d, &
        assemble_helmholtz_single_layer_p0_3d_jvp, &
        assemble_helmholtz_single_layer_p0_3d_vjp, &
        assemble_laplace_single_layer_p0_3d
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: quadrature_degree = 4
    real(dp), parameter :: step = 2.0e-7_dp
    integer, parameter :: triangles(3, 2) = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    real(dp) :: wave_number, wave_number_bar, wave_number_dot
    real(dp), allocatable :: laplace_matrix(:, :)
    complex(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    complex(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    complex(dp) :: matrix_bar(2, 2)
    real(dp) :: area_vector(3), areas(2), lhs, rhs
    integer :: first, status, status_minus, status_plus

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, &
        1.2_dp, 0.15_dp, -0.1_dp, &
        -0.25_dp, 1.1_dp, 0.2_dp, &
        0.3_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, &
        -0.015_dp, 0.025_dp, -0.02_dp, &
        0.01_dp, 0.035_dp, -0.005_dp, &
        -0.01_dp, 0.02_dp, 0.015_dp], [3, 4])
    wave_number = 2.3_dp
    wave_number_dot = -0.17_dp
    matrix_bar = reshape([ &
        cmplx(0.4_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.25_dp, dp), &
        cmplx(0.2_dp, 0.15_dp, dp), cmplx(0.6_dp, -0.35_dp, dp)], [2, 2])

    call assemble_helmholtz_single_layer_p0_3d( &
        vertices, triangles, wave_number, quadrature_degree, matrix, status)
    call assemble_helmholtz_single_layer_p0_3d_jvp( &
        vertices, triangles, wave_number, quadrature_degree, vertices_dot, &
        wave_number_dot, matrix_dot, status)
    call assemble_helmholtz_single_layer_p0_3d( &
        vertices + step*vertices_dot, triangles, &
        wave_number + step*wave_number_dot, quadrature_degree, matrix_plus, &
        status_plus)
    call assemble_helmholtz_single_layer_p0_3d( &
        vertices - step*vertices_dot, triangles, &
        wave_number - step*wave_number_dot, quadrature_degree, matrix_minus, &
        status_minus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Helmholtz BEM assembly JVP succeeds")
    call check_condition(maxval(abs( &
        matrix_dot - (matrix_plus - matrix_minus)/(2.0_dp*step))) < &
        6.0e-9_dp, "Helmholtz BEM assembly JVP matches central differences")
    call check_condition( &
        maxval(abs(matrix_dot - transpose(matrix_dot))) < 2.0e-15_dp, &
        "Helmholtz BEM assembly JVP preserves matrix symmetry")

    call assemble_helmholtz_single_layer_p0_3d_vjp( &
        vertices, triangles, wave_number, quadrature_degree, matrix_bar, &
        vertices_bar, wave_number_bar, status)
    lhs = real(sum(conjg(matrix_bar)*matrix_dot), dp)
    rhs = sum(vertices_bar*vertices_dot) + wave_number_bar*wave_number_dot
    call check_condition(status == 0, "Helmholtz BEM assembly VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 6.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Helmholtz BEM assembly products obey the adjoint identity")

    call assemble_laplace_single_layer_p0_3d( &
        vertices, triangles, quadrature_degree, laplace_matrix, status)
    call assemble_helmholtz_single_layer_p0_3d( &
        vertices, triangles, 0.0_dp, quadrature_degree, matrix, status_plus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. &
        maxval(abs(matrix - cmplx(laplace_matrix, 0.0_dp, dp))) < 5.0e-15_dp, &
        "zero-wave-number Helmholtz assembly reduces to Laplace")
    do first = 1, 2
        area_vector = cross_product( &
            vertices(:, triangles(2, first)) - &
            vertices(:, triangles(1, first)), &
            vertices(:, triangles(3, first)) - &
            vertices(:, triangles(1, first)))
        areas(first) = 0.5_dp*norm2(area_vector)
    end do
    call assemble_helmholtz_single_layer_p0_3d_jvp( &
        vertices, triangles, 0.0_dp, quadrature_degree, 0.0_dp*vertices_dot, &
        1.0_dp, matrix_dot, status)
    call check_condition( &
        status == 0 .and. maxval(abs(matrix_dot - cmplx(0.0_dp, &
        spread(areas, 2, 2)*spread(areas, 1, 2)/(4.0_dp*acos(-1.0_dp)), &
        dp))) < 2.0e-14_dp, &
        "zero-wave-number derivative matches the constant kernel integral")

    call check_summary("Helmholtz BEM assembly derivatives")

contains

    pure function cross_product(first_vector, second_vector) result(product)
        real(dp), intent(in) :: first_vector(3), second_vector(3)
        real(dp) :: product(3)

        product = [ &
            first_vector(2)*second_vector(3) - &
            first_vector(3)*second_vector(2), &
            first_vector(3)*second_vector(1) - &
            first_vector(1)*second_vector(3), &
            first_vector(1)*second_vector(2) - &
            first_vector(2)*second_vector(1)]
    end function cross_product
end program test_helmholtz_bem_assembly_ad
