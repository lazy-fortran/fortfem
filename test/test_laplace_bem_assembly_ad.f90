program test_laplace_bem_assembly_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_laplace_single_layer_p0_3d, &
        assemble_laplace_single_layer_p0_3d_jvp, &
        assemble_laplace_single_layer_p0_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: quadrature_degree = 4
    real(dp), parameter :: step = 2.0e-7_dp
    integer, parameter :: triangles(3, 2) = reshape([1, 2, 3, 1, 3, 4], [3, 2])
    real(dp) :: vertices(3, 4), vertices_bar(3, 4), vertices_dot(3, 4)
    real(dp), allocatable :: matrix(:, :), matrix_dot(:, :)
    real(dp), allocatable :: matrix_minus(:, :), matrix_plus(:, :)
    real(dp) :: matrix_bar(2, 2), lhs, rhs
    integer :: status, status_minus, status_plus

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
    matrix_bar = reshape([0.4_dp, -0.3_dp, 0.2_dp, 0.6_dp], [2, 2])

    call assemble_laplace_single_layer_p0_3d( &
        vertices, triangles, quadrature_degree, matrix, status)
    call assemble_laplace_single_layer_p0_3d_jvp( &
        vertices, triangles, quadrature_degree, vertices_dot, matrix_dot, &
        status)
    call assemble_laplace_single_layer_p0_3d( &
        vertices + step*vertices_dot, triangles, quadrature_degree, &
        matrix_plus, status_plus)
    call assemble_laplace_single_layer_p0_3d( &
        vertices - step*vertices_dot, triangles, quadrature_degree, &
        matrix_minus, status_minus)
    call check_condition( &
        status == 0 .and. status_plus == 0 .and. status_minus == 0, &
        "Laplace BEM assembly JVP succeeds")
    call check_condition(maxval(abs( &
        matrix_dot - (matrix_plus - matrix_minus)/(2.0_dp*step))) < &
        4.0e-9_dp, "Laplace BEM assembly JVP matches central differences")
    call check_condition( &
        maxval(abs(matrix_dot - transpose(matrix_dot))) < 2.0e-15_dp, &
        "Laplace BEM assembly JVP preserves matrix symmetry")

    call assemble_laplace_single_layer_p0_3d_vjp( &
        vertices, triangles, quadrature_degree, matrix_bar, vertices_bar, &
        status)
    lhs = sum(matrix_bar*matrix_dot)
    rhs = sum(vertices_bar*vertices_dot)
    call check_condition(status == 0, "Laplace BEM assembly VJP succeeds")
    call check_condition( &
        abs(lhs - rhs) < 4.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "Laplace BEM assembly products obey the adjoint identity")

    call check_summary("Laplace BEM assembly derivatives")
end program test_laplace_bem_assembly_ad
