program test_triangle_affine_map_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_core, only: invert_triangle_affine_map, &
        invert_triangle_affine_map_jvp, invert_triangle_affine_map_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: lhs, point(2), point_bar(2), point_dot(2)
    real(dp) :: reference(2), reference_bar(2), reference_dot(2)
    real(dp) :: reference_minus(2), reference_plus(2), rhs
    real(dp) :: vertices(2, 3), vertices_bar(2, 3), vertices_dot(2, 3)
    integer :: failures, status

    failures = 0
    vertices = reshape([ &
        0.1_dp, -0.2_dp, 1.4_dp, 0.15_dp, -0.3_dp, 1.2_dp], [2, 3])
    point = [0.37_dp, 0.28_dp]
    vertices_dot = reshape([ &
        0.03_dp, -0.01_dp, -0.02_dp, 0.04_dp, 0.01_dp, -0.025_dp], [2, 3])
    point_dot = [-0.015_dp, 0.022_dp]
    reference_bar = [0.8_dp, -0.45_dp]

    call invert_triangle_affine_map( &
        vertices, point, reference, status)
    call check(status == 0, "affine inverse succeeds")
    call check(maxval(abs( &
        vertices(:, 1) + matmul( &
        reshape([ &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)], [2, 2]), reference) - point)) < &
        2.0e-15_dp, "affine inverse reconstructs physical point")
    call invert_triangle_affine_map_jvp( &
        vertices, point, vertices_dot, point_dot, reference_dot, status)
    call invert_triangle_affine_map( &
        vertices + step*vertices_dot, point + step*point_dot, reference_plus, &
        status)
    call invert_triangle_affine_map( &
        vertices - step*vertices_dot, point - step*point_dot, reference_minus, &
        status)
    call check(maxval(abs( &
        reference_dot - (reference_plus - reference_minus)/(2.0_dp*step))) < &
        2.0e-9_dp, "affine inverse JVP matches full re-inversion")
    call invert_triangle_affine_map_vjp( &
        vertices, point, reference_bar, vertices_bar, point_bar, status)
    lhs = dot_product(reference_bar, reference_dot)
    rhs = sum(vertices_bar*vertices_dot) + dot_product(point_bar, point_dot)
    call check(abs(lhs - rhs) < 2.0e-13_dp, &
        "affine inverse products obey adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_triangle_affine_map_ad
