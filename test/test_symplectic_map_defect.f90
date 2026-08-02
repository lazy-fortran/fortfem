program test_symplectic_map_defect
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_symplectic_map_defect, &
        assemble_symplectic_map_defect_jvp, assemble_symplectic_map_defect_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: step = 1.0e-7_dp
    real(dp), parameter :: map(4, 4) = reshape([ &
        1.0_dp, 0.0_dp, 0.2_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.3_dp, &
        0.1_dp, 0.0_dp, 1.1_dp, 0.0_dp, &
        0.0_dp, -0.2_dp, 0.0_dp, 0.9_dp], [4, 4])
    real(dp), parameter :: symplectic_form(4, 4) = reshape([ &
        0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp], [4, 4])
    real(dp), parameter :: map_dot(4, 4) = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, -0.6_dp, &
        0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 1.1_dp, -1.2_dp, &
        1.3_dp, -1.4_dp, 1.5_dp, -1.6_dp], [4, 4])
    real(dp), parameter :: form_dot(4, 4) = reshape([ &
        0.0_dp, 0.1_dp, -0.2_dp, 0.3_dp, -0.4_dp, 0.5_dp, &
        -0.6_dp, 0.7_dp, 0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, &
        1.2_dp, -1.3_dp, 1.4_dp, -1.5_dp], [4, 4])
    real(dp) :: defect(4, 4), defect_dot(4, 4), defect_expected(4, 4)
    real(dp) :: defect_plus(4, 4), defect_minus(4, 4)
    real(dp) :: map_bar(4, 4), form_bar(4, 4), defect_bar(4, 4)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_symplectic_map_defect( &
        map, symplectic_form, defect, status)
    defect_expected = matmul(transpose(map), matmul(symplectic_form, map)) - &
        symplectic_form
    call check_condition(status%code == 0 .and. &
        maxval(abs(defect - defect_expected)) < 1.0e-14_dp, &
        "symplectic-map defect matches the independent matrix oracle")

    call assemble_symplectic_map_defect_jvp( &
        map, symplectic_form, map_dot, form_dot, defect_dot, status)
    defect_expected = matmul(transpose(map_dot), matmul(symplectic_form, map)) + &
        matmul(transpose(map), matmul(form_dot, map)) + &
        matmul(transpose(map), matmul(symplectic_form, map_dot)) - form_dot
    call check_condition(status%code == 0 .and. &
        maxval(abs(defect_dot - defect_expected)) < 1.0e-14_dp, &
        "symplectic-map JVP applies the full map and form product rule")

    call assemble_symplectic_map_defect( &
        map + step*map_dot, symplectic_form + step*form_dot, defect_plus, status)
    call assemble_symplectic_map_defect( &
        map - step*map_dot, symplectic_form - step*form_dot, defect_minus, status)
    call check_condition(maxval(abs(defect_dot - &
        (defect_plus - defect_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "symplectic-map JVP matches an independent central difference")

    defect_bar = reshape([0.2_dp, -0.3_dp, 0.4_dp, -0.5_dp, 0.6_dp, -0.7_dp, &
        0.8_dp, -0.9_dp, 1.0_dp, -1.1_dp, 1.2_dp, -1.3_dp, 1.4_dp, -1.5_dp, &
        1.6_dp, -1.7_dp], [4, 4])
    call assemble_symplectic_map_defect_vjp( &
        map, symplectic_form, defect_bar, map_bar, form_bar, status)
    lhs = sum(defect_bar*defect_dot)
    rhs = sum(map_bar*map_dot) + sum(form_bar*form_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-13_dp, &
        "symplectic-map VJP satisfies the real dot-product identity")

    call check_summary("symplectic map defect")
end program test_symplectic_map_defect
