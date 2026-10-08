program test_tetra_piola_conservation
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_tetra_piola_maps, only: map_tetra_nedelec_covariant, &
        map_tetra_nedelec_covariant_jvp, map_tetra_rt_contravariant, &
        map_tetra_rt_contravariant_jvp
    implicit none
    real(dp) :: j(3,3), jd(3,3), v(3,2), vd(3,2), c(3,2), cd(3,2)
    real(dp) :: w(3,2), wd(3,2), k(3,2), kd(3,2), s(2), sd(2), t(2), td(2)
    real(dp) :: a(3), b(3), n(3), pn(3), pnd(3), det, detd
    integer :: sample, status
    v = reshape([2.0_dp,-3.0_dp,4.0_dp,-5.0_dp,7.0_dp,-2.0_dp],[3,2])
    vd = reshape([0.3_dp,-0.7_dp,0.1_dp,0.2_dp,0.9_dp,-0.5_dp],[3,2])
    c = reshape([11.0_dp,-13.0_dp,2.0_dp,-3.0_dp,4.0_dp,5.0_dp],[3,2])
    cd = reshape([-0.2_dp,0.4_dp,0.7_dp,0.1_dp,-0.3_dp,0.5_dp],[3,2])
    s = [11.0_dp,-13.0_dp]; sd = [-0.2_dp,0.4_dp]
    jd = reshape([0.1_dp,-0.2_dp,0.05_dp,0.3_dp,0.4_dp,-0.1_dp, &
        -0.15_dp,0.25_dp,0.2_dp],[3,3])
    a = [0.6_dp,-0.8_dp,0.3_dp]; b = [-0.2_dp,0.5_dp,0.9_dp]
    n = cross(a,b)
    do sample = 1,9
        j = reshape([1.0_dp+sample/10.0_dp,0.2_dp,-0.1_dp, &
            -0.3_dp,1.2_dp+sample/20.0_dp,0.15_dp, &
            0.12_dp,-0.22_dp,1.1_dp+sample/30.0_dp],[3,3])
        det = dot_product(j(:,1),cross(j(:,2),j(:,3)))
        detd = dot_product(jd(:,1),cross(j(:,2),j(:,3))) + &
            dot_product(j(:,1),cross(jd(:,2),j(:,3))+cross(j(:,2),jd(:,3)))
        pn = cross(matmul(j,a),matmul(j,b))
        pnd = cross(matmul(jd,a),matmul(j,b))+cross(matmul(j,a),matmul(jd,b))
        call map_tetra_nedelec_covariant(j,v,c,w,k,status)
        call check_condition(status==0,'Covariant positive volume accepted')
        call check_condition(maxval(abs(matmul(transpose(w),matmul(j,a)) - &
            matmul(transpose(v),a)))<1.0e-12_dp,'Oriented line moments preserved')
        call check_condition(maxval(abs(matmul(transpose(k),pn) - &
            matmul(transpose(c),n)))<1.0e-12_dp,'Oriented curl flux preserved')
        call map_tetra_nedelec_covariant_jvp(j,v,c,jd,vd,cd,wd,kd,status)
        call check_condition(status==0,'Covariant tangent accepted')
        call check_condition(maxval(abs(matmul(transpose(wd),matmul(j,a)) + &
            matmul(transpose(w),matmul(jd,a))-matmul(transpose(vd),a))) &
            <1.0e-12_dp,'Differentiated line conservation')
        call check_condition(maxval(abs(matmul(transpose(kd),pn) + &
            matmul(transpose(k),pnd)-matmul(transpose(cd),n))) &
            <1.0e-12_dp,'Differentiated curl flux conservation')
        call map_tetra_rt_contravariant(j,v,s,w,t,status)
        call check_condition(status==0,'Contravariant positive volume accepted')
        call check_condition(maxval(abs(matmul(transpose(w),pn) - &
            matmul(transpose(v),n)))<1.0e-12_dp,'Oriented normal flux preserved')
        call check_condition(maxval(abs(det*t-s))<1.0e-12_dp, &
            'Divergence volume integral preserved')
        call map_tetra_rt_contravariant_jvp(j,v,s,jd,vd,sd,wd,td,status)
        call check_condition(status==0,'Contravariant tangent accepted')
        call check_condition(maxval(abs(matmul(transpose(wd),pn) + &
            matmul(transpose(w),pnd)-matmul(transpose(vd),n))) &
            <1.0e-12_dp,'Differentiated normal flux conservation')
        call check_condition(maxval(abs(det*td+detd*t-sd))<1.0e-12_dp, &
            'Differentiated divergence volume conservation')
    end do
    call check_summary('Tetrahedron Piola independent conservation')
contains
    pure function cross(x,y) result(z)
        real(dp),intent(in) :: x(3),y(3)
        real(dp) :: z(3)
        z=[x(2)*y(3)-x(3)*y(2),x(3)*y(1)-x(1)*y(3),x(1)*y(2)-x(2)*y(1)]
    end function cross
end program test_tetra_piola_conservation
