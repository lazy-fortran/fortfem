program test_legacy_affine_basis_geometry
    use fortfem_kinds,only:dp
    use basis_p1_2d_module,only:basis_p1_2d_t
    use basis_p2_2d_module,only:basis_p2_2d_t
    use check,only:check_condition,check_summary
    implicit none
    type(basis_p1_2d_t)::p1
    type(basis_p2_2d_t)::p2
    real(dp)::vertices(2,3),expected_jac(2,2),jac(2,2),det,expected_det,translation(2),scale,x,y,expected(2)
    integer::shift,scaling,orientation
    integer,parameter::shifts(4)=[0,20,30,40]
    do shift=1,4
        do scaling=0,2
            do orientation=1,2
                scale=2.0_dp**(-2*scaling)
                translation=[2.0_dp**shifts(shift),-2.0_dp**shifts(shift)]
                expected_jac=scale*reshape([1.0_dp,.5_dp,.25_dp,1.5_dp],[2,2])
                if(orientation==2)expected_jac(:,1)=-expected_jac(:,1)
                vertices(:,1)=translation
                vertices(:,2)=translation+expected_jac(:,1)
                vertices(:,3)=translation+expected_jac(:,2)
                expected_det=1.375_dp*scale*scale
                if(orientation==2)expected_det=-expected_det
                call p1%compute_jacobian(vertices,jac,det)
                call check_condition(maxval(abs(jac-expected_jac))<2e-13_dp*scale, &
                    'P1 Jacobian invariant under exact large translation and signed orientation')
                call check_condition(abs(det-expected_det)<2e-13_dp*scale*scale,'P1 signed determinant is exact affine area')
                call p2%compute_jacobian(vertices,jac,det)
                if(maxval(abs(jac-expected_jac))>=2e-13_dp*scale) &
                    print *,"P2 geometry failure: shift,scale,orientation,Jerr,deterr", &
                    shifts(shift),scale,orientation,maxval(abs(jac-expected_jac)),abs(det-expected_det)
                call check_condition(maxval(abs(jac-expected_jac))<2e-13_dp*scale, &
                    'P2 straight-sided Jacobian invariant under exact large translation')
                call check_condition(abs(det-expected_det)<2e-13_dp*scale*scale,'P2 signed determinant is exact affine area')
                expected=translation+matmul(expected_jac,[.25_dp,.125_dp])
                call p1%transform_to_physical(.25_dp,.125_dp,vertices,x,y)
                call check_condition(all([x,y]==expected),'P1 exactly representable forward affine map')
                call p2%transform_to_physical(.25_dp,.125_dp,vertices,x,y)
                call check_condition(all([x,y]==expected),'P2 exactly representable straight-sided map')
            end do
        end do
    end do
    call check_summary('Legacy scalar affine geometry invariant')
end program
