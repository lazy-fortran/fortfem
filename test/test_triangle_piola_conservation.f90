program test_triangle_piola_conservation
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_triangle_piola_maps, only: map_triangle_nedelec_covariant, &
        map_triangle_nedelec_covariant_jvp, map_triangle_rt_contravariant, &
        map_triangle_rt_contravariant_jvp
    implicit none
    real(dp) :: j(2, 2), jd(2, 2), v(2, 2), vd(2, 2), s(2), sd(2)
    real(dp) :: w(2, 2), wd(2, 2), t(2), td(2), det, detd
    real(dp) :: tangent(2), normal(2), physical_tangent(2), physical_normal(2)
    integer :: sample, status

    v = reshape([2.0_dp, -3.0_dp, -5.0_dp, 7.0_dp], [2, 2])
    vd = reshape([0.3_dp, -0.7_dp, 0.2_dp, 0.9_dp], [2, 2])
    s = [11.0_dp, -13.0_dp]
    sd = [-0.2_dp, 0.4_dp]
    jd = reshape([0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
    tangent = [0.6_dp, -0.8_dp]
    normal = [-tangent(2), tangent(1)]
    do sample = 1, 9
        j = reshape([1.0_dp + sample/10.0_dp, 0.2_dp, &
            -0.3_dp, 1.2_dp + sample/20.0_dp], [2, 2])
        det = j(1, 1)*j(2, 2) - j(1, 2)*j(2, 1)
        detd = jd(1, 1)*j(2, 2) + j(1, 1)*jd(2, 2) - &
            jd(1, 2)*j(2, 1) - j(1, 2)*jd(2, 1)
        physical_tangent = matmul(j, tangent)
        physical_normal = [-physical_tangent(2), physical_tangent(1)]

        call map_triangle_nedelec_covariant(j, v, s, w, t, status)
        call check_condition(status == 0, "Covariant map accepts positive area")
        call check_condition(maxval(abs(matmul(transpose(w), physical_tangent) &
            - matmul(transpose(v), tangent))) < 1.0e-13_dp, &
            "Covariant map preserves oriented tangential line moments")
        call check_condition(maxval(abs(det*t - s)) < 1.0e-13_dp, &
            "Covariant map preserves oriented curl area integrals")
        call map_triangle_nedelec_covariant_jvp(j, v, s, jd, vd, sd, wd, td, status)
        call check_condition(status == 0, "Covariant tangent succeeds")
        call check_condition(maxval(abs(matmul(transpose(j), wd) + &
            matmul(transpose(jd), w) - vd)) < 1.0e-13_dp, &
            "Covariant tangent satisfies differentiated line conservation")
        call check_condition(maxval(abs(det*td + detd*t - sd)) < 1.0e-13_dp, &
            "Covariant tangent satisfies differentiated curl conservation")

        call map_triangle_rt_contravariant(j, v, s, w, t, status)
        call check_condition(status == 0, "Contravariant map accepts positive area")
        call check_condition(maxval(abs(matmul(transpose(w), physical_normal) &
            - matmul(transpose(v), normal))) < 1.0e-13_dp, &
            "Contravariant map preserves oriented normal flux moments")
        call check_condition(maxval(abs(det*t - s)) < 1.0e-13_dp, &
            "Contravariant map preserves divergence area integrals")
        call map_triangle_rt_contravariant_jvp(j, v, s, jd, vd, sd, wd, td, status)
        call check_condition(status == 0, "Contravariant tangent succeeds")
        call check_condition(maxval(abs(det*wd + detd*w - &
            matmul(jd, v) - matmul(j, vd))) < 1.0e-13_dp, &
            "Contravariant tangent satisfies differentiated flux conservation")
        call check_condition(maxval(abs(det*td + detd*t - sd)) < 1.0e-13_dp, &
            "Contravariant tangent satisfies differentiated divergence conservation")
    end do
    call check_summary("Triangle Piola independent conservation")
end program test_triangle_piola_conservation
