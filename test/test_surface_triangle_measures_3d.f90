program test_surface_triangle_measures_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_surface_triangle_measures_3d, &
        assemble_surface_triangle_measures_3d_jvp, &
        assemble_surface_triangle_measures_3d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 4), vertices_dot(3, 4), vertices_bar(3, 4)
    integer, parameter :: triangles(3, 2) = reshape([1, 2, 3, 1, 4, 2], [3, 2])
    real(dp), allocatable :: areas(:), areas_dot(:), areas_minus(:), areas_plus(:)
    real(dp), allocatable :: normals(:, :), normals_dot(:, :)
    real(dp), allocatable :: normals_minus(:, :), normals_plus(:, :)
    real(dp) :: area_bar(2), normal_bar(3, 2), lhs, rhs
    integer :: status, status_minus, status_plus

    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    vertices_dot = reshape([ &
        0.02_dp, -0.01_dp, 0.03_dp, &
        -0.03_dp, 0.04_dp, 0.01_dp, &
        0.01_dp, 0.02_dp, -0.02_dp, &
        -0.015_dp, 0.01_dp, 0.025_dp], [3, 4])
    area_bar = [0.3_dp, -0.2_dp]
    normal_bar = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, &
        -0.4_dp, 0.5_dp, -0.6_dp], [3, 2])

    call assemble_surface_triangle_measures_3d( &
        vertices, triangles, areas, normals, status)
    call check_condition(status == 0, &
        "oriented surface measures accept a valid triangle mesh")
    call check_condition(maxval(abs(areas - [0.5_dp, 0.5_dp])) < 1.0e-14_dp, &
        "surface measures return exact triangle areas")
    call check_condition(maxval(abs(normals - reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], [3, 2]))) < &
        1.0e-14_dp, "surface measures preserve oriented unit normals")

    call assemble_surface_triangle_measures_3d_jvp( &
        vertices, triangles, vertices_dot, areas, areas_dot, normals, &
        normals_dot, status)
    call assemble_surface_triangle_measures_3d( &
        vertices - step*vertices_dot, triangles, areas_minus, normals_minus, &
        status_minus)
    call assemble_surface_triangle_measures_3d( &
        vertices + step*vertices_dot, triangles, areas_plus, normals_plus, &
        status_plus)
    call check_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "surface measure JVP accepts a valid direction")
    call check_condition(maxval(abs(areas_dot - &
        (areas_plus - areas_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "surface area JVP matches central differences")
    call check_condition(maxval(abs(normals_dot - &
        (normals_plus - normals_minus)/(2.0_dp*step))) < 3.0e-9_dp, &
        "surface normal JVP matches central differences")

    call assemble_surface_triangle_measures_3d_vjp( &
        vertices, triangles, area_bar, normal_bar, vertices_bar, status)
    lhs = dot_product(area_bar, areas_dot) + sum(normal_bar*normals_dot)
    rhs = sum(vertices_bar*vertices_dot)
    call check_condition(status == 0, "surface measure VJP accepts valid bars")
    call check_condition(abs(lhs - rhs) < &
        3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "surface measure VJP satisfies the real dot-product identity")

    vertices(:, 4) = vertices(:, 1)
    call assemble_surface_triangle_measures_3d( &
        vertices, triangles, areas, normals, status)
    call check_condition(status /= 0, "surface measures reject a degenerate face")
    call check_summary("Oriented surface triangle measures")
end program test_surface_triangle_measures_3d
