program test_feec_commuting_projection
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_feec_commuting_projection, &
        assemble_feec_commuting_projection_jvp, &
        assemble_feec_commuting_projection_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp) :: gd(3, 2), cd(2, 3), dd(1, 2)
    real(dp) :: gc(3, 2), cc(2, 3), dc(1, 2)
    real(dp) :: p0(2, 2), p1(3, 3), p2(2, 2), p3(1, 1)
    real(dp) :: gd_dot(3, 2), cd_dot(2, 3), dd_dot(1, 2)
    real(dp) :: gc_dot(3, 2), cc_dot(2, 3), dc_dot(1, 2)
    real(dp) :: p0_dot(2, 2), p1_dot(3, 3), p2_dot(2, 2), p3_dot(1, 1)
    real(dp) :: gd_bar(3, 2), cd_bar(2, 3), dd_bar(1, 2)
    real(dp) :: gc_bar(3, 2), cc_bar(2, 3), dc_bar(1, 2)
    real(dp) :: p0_bar(2, 2), p1_bar(3, 3), p2_bar(2, 2), p3_bar(1, 1)
    real(dp) :: defect_g(3, 2), defect_c(2, 3), defect_d(1, 2)
    real(dp) :: defect_g_dot(3, 2), defect_c_dot(2, 3), defect_d_dot(1, 2)
    real(dp) :: defect_g_bar(3, 2), defect_c_bar(2, 3), defect_d_bar(1, 2)
    real(dp) :: expected_g(3, 2), expected_c(2, 3), expected_d(1, 2)
    real(dp) :: plus_g(3, 2), plus_c(2, 3), plus_d(1, 2)
    real(dp) :: minus_g(3, 2), minus_c(2, 3), minus_d(1, 2)
    real(dp) :: lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.

    gd = reshape([1.0_dp, -2.0_dp, 0.5_dp, 0.25_dp, 1.5_dp, -0.75_dp], [3, 2])
    cd = reshape([0.4_dp, -1.0_dp, 2.0_dp, 0.8_dp, 1.25_dp, -0.5_dp], [2, 3])
    dd = reshape([1.1_dp, -0.7_dp], [1, 2])
    gc = reshape([0.3_dp, 1.2_dp, -0.4_dp, 0.9_dp, 1.1_dp, -1.3_dp], [3, 2])
    cc = reshape([-0.8_dp, 0.6_dp, 1.4_dp, 0.2_dp, -1.1_dp, 0.5_dp], [2, 3])
    dc = reshape([0.75_dp, -1.4_dp], [1, 2])
    p0 = reshape([1.0_dp, 0.2_dp, -0.1_dp, 0.9_dp], [2, 2])
    p1 = reshape([1.1_dp, 0.0_dp, 0.2_dp, -0.1_dp, 0.8_dp, 0.3_dp, &
        0.05_dp, -0.2_dp, 1.2_dp], [3, 3])
    p2 = reshape([0.9_dp, -0.15_dp, 0.4_dp, 1.05_dp], [2, 2])
    p3 = reshape([1.25_dp], [1, 1])

    gd_dot = reshape([-0.2_dp, 0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, 0.6_dp], [3, 2])
    cd_dot = reshape([0.15_dp, -0.25_dp, 0.35_dp, 0.45_dp, -0.55_dp, 0.65_dp], [2, 3])
    dd_dot = reshape([-0.3_dp, 0.7_dp], [1, 2])
    gc_dot = reshape([0.45_dp, -0.35_dp, 0.25_dp, 0.15_dp, -0.05_dp, 0.55_dp], [3, 2])
    cc_dot = reshape([-0.65_dp, 0.75_dp, -0.85_dp, 0.95_dp, 0.1_dp, -0.2_dp], [2, 3])
    dc_dot = reshape([0.8_dp, -0.9_dp], [1, 2])
    p0_dot = reshape([0.12_dp, -0.22_dp, 0.32_dp, -0.42_dp], [2, 2])
    p1_dot = reshape([0.11_dp, -0.21_dp, 0.31_dp, -0.41_dp, 0.51_dp, -0.61_dp, &
        0.71_dp, -0.81_dp, 0.91_dp], [3, 3])
    p2_dot = reshape([-0.13_dp, 0.23_dp, -0.33_dp, 0.43_dp], [2, 2])
    p3_dot = reshape([0.53_dp], [1, 1])

    call assemble_feec_commuting_projection( &
        gd, cd, dd, gc, cc, dc, p0, p1, p2, p3, &
        defect_g, defect_c, defect_d, status)
    expected_g = matmul(gd, p0) - matmul(p1, gc)
    expected_c = matmul(cd, p1) - matmul(p2, cc)
    expected_d = matmul(dd, p2) - matmul(p3, dc)
    call record_condition(status%code == 0 .and. &
        maxval(abs(defect_g - expected_g)) < 1.0e-14_dp .and. &
        maxval(abs(defect_c - expected_c)) < 1.0e-14_dp .and. &
        maxval(abs(defect_d - expected_d)) < 1.0e-14_dp, &
        "commuting defects match independent matrix-product oracle")

    call assemble_feec_commuting_projection_jvp( &
        gd, cd, dd, gc, cc, dc, p0, p1, p2, p3, gd_dot, cd_dot, dd_dot, &
        gc_dot, cc_dot, dc_dot, p0_dot, p1_dot, p2_dot, p3_dot, &
        defect_g_dot, defect_c_dot, defect_d_dot, status)
    expected_g = matmul(gd_dot, p0) + matmul(gd, p0_dot) - &
        matmul(p1_dot, gc) - matmul(p1, gc_dot)
    expected_c = matmul(cd_dot, p1) + matmul(cd, p1_dot) - &
        matmul(p2_dot, cc) - matmul(p2, cc_dot)
    expected_d = matmul(dd_dot, p2) + matmul(dd, p2_dot) - &
        matmul(p3_dot, dc) - matmul(p3, dc_dot)
    call record_condition(status%code == 0 .and. &
        maxval(abs(defect_g_dot - expected_g)) < 1.0e-14_dp .and. &
        maxval(abs(defect_c_dot - expected_c)) < 1.0e-14_dp .and. &
        maxval(abs(defect_d_dot - expected_d)) < 1.0e-14_dp, &
        "commuting defect JVP applies every product-rule term")

    call assemble_feec_commuting_projection( &
        gd + eps*gd_dot, cd + eps*cd_dot, dd + eps*dd_dot, &
        gc + eps*gc_dot, cc + eps*cc_dot, dc + eps*dc_dot, &
        p0 + eps*p0_dot, p1 + eps*p1_dot, p2 + eps*p2_dot, p3 + eps*p3_dot, &
        plus_g, plus_c, plus_d, status)
    call assemble_feec_commuting_projection( &
        gd - eps*gd_dot, cd - eps*cd_dot, dd - eps*dd_dot, &
        gc - eps*gc_dot, cc - eps*cc_dot, dc - eps*dc_dot, &
        p0 - eps*p0_dot, p1 - eps*p1_dot, p2 - eps*p2_dot, p3 - eps*p3_dot, &
        minus_g, minus_c, minus_d, status)
    call record_condition( &
        maxval(abs((plus_g - minus_g)/(2.0_dp*eps) - defect_g_dot)) < 2.0e-8_dp .and. &
        maxval(abs((plus_c - minus_c)/(2.0_dp*eps) - defect_c_dot)) < 2.0e-8_dp .and. &
        maxval(abs((plus_d - minus_d)/(2.0_dp*eps) - defect_d_dot)) < 2.0e-8_dp, &
        "commuting defect JVP agrees with a complete central difference")

    defect_g_bar = reshape([0.4_dp, -0.3_dp, 0.2_dp, -0.1_dp, 0.6_dp, -0.5_dp], [3, 2])
    defect_c_bar = reshape([0.7_dp, -0.8_dp, 0.9_dp, -1.0_dp, 0.3_dp, 0.2_dp], [2, 3])
    defect_d_bar = reshape([-0.6_dp, 0.5_dp], [1, 2])
    call assemble_feec_commuting_projection_vjp( &
        gd, cd, dd, gc, cc, dc, p0, p1, p2, p3, defect_g_bar, defect_c_bar, &
        defect_d_bar, gd_bar, cd_bar, dd_bar, gc_bar, cc_bar, dc_bar, &
        p0_bar, p1_bar, p2_bar, p3_bar, status)
    lhs = sum(defect_g_bar*defect_g_dot) + sum(defect_c_bar*defect_c_dot) + &
        sum(defect_d_bar*defect_d_dot)
    rhs = sum(gd_bar*gd_dot) + sum(cd_bar*cd_dot) + sum(dd_bar*dd_dot) + &
        sum(gc_bar*gc_dot) + sum(cc_bar*cc_dot) + sum(dc_bar*dc_dot) + &
        sum(p0_bar*p0_dot) + sum(p1_bar*p1_dot) + sum(p2_bar*p2_dot) + &
        sum(p3_bar*p3_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "commuting defect VJP satisfies the real dot-product identity")

    gd = reshape([-1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [3, 2])
    cd = reshape([1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp], [2, 3])
    dd = reshape([2.0_dp, -1.0_dp], [1, 2])
    gc = gd
    cc = cd
    dc = dd
    p0 = 0.0_dp
    p1 = 0.0_dp
    p2 = 0.0_dp
    p3 = 0.0_dp
    p0(1, 1) = 1.0_dp
    p0(2, 2) = 1.0_dp
    p1(1, 1) = 1.0_dp
    p1(2, 2) = 1.0_dp
    p1(3, 3) = 1.0_dp
    p2(1, 1) = 1.0_dp
    p2(2, 2) = 1.0_dp
    p3(1, 1) = 1.0_dp
    call assemble_feec_commuting_projection( &
        gd, cd, dd, gc, cc, dc, p0, p1, p2, p3, &
        defect_g, defect_c, defect_d, status)
    call record_condition(status%code == 0 .and. &
        maxval(abs(defect_g)) < 1.0e-14_dp .and. &
        maxval(abs(defect_c)) < 1.0e-14_dp .and. &
        maxval(abs(defect_d)) < 1.0e-14_dp, &
        "identity projections preserve an exact discrete/continuous diagram")

    call assemble_feec_commuting_projection( &
        gd, cd, dd, gc, cc, dc, p0, p1(:, 1:2), p2, p3, &
        defect_g, defect_c, defect_d, status)
    call record_condition(status%code /= 0, &
        "commuting projection rejects incompatible projection dimensions")
    call check_summary("FEEC commuting projection")

    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_feec_commuting_projection
