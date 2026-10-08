program test_triangle_vector_boundary_jets
 use check,only:check_condition,check_summary
 use fortfem_kinds,only:dp
 use fortfem_triangle_nedelec_arbitrary_order,only:triangle_nedelec_first_kind_t, &
 initialize_triangle_nedelec_first_kind,evaluate_triangle_nedelec_first_kind, &
 evaluate_triangle_nedelec_first_kind_jvp,evaluate_triangle_nedelec_first_kind_vjp,triangle_nedelec_dof_count
 use fortfem_triangle_nedelec_second_kind,only:triangle_nedelec_second_kind_t, &
 initialize_triangle_nedelec_second_kind,evaluate_triangle_nedelec_second_kind, &
 evaluate_triangle_nedelec_second_kind_jvp,evaluate_triangle_nedelec_second_kind_vjp, &
 triangle_nedelec_second_kind_dof_count
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 type(triangle_nedelec_first_kind_t)::first
 type(triangle_nedelec_second_kind_t)::second
 real(dp),allocatable::v(:,:),v1(:,:),v2(:,:),vd(:,:),vb(:,:),c(:),c1(:),c2(:),cd(:),cb(:)
 real(dp)::x(2),direction(2),bar(2),lhs,rhs,error,scale
 real(dp),parameter::h=1e-5_dp
 integer::family,degree,sample,n,status,plus_status,twice_status,i
 do family=1,2
 do degree=1,4
 if(family==1)then
 call initialize_triangle_nedelec_first_kind(degree,first,status)
 n=triangle_nedelec_dof_count(first)
 else
 call initialize_triangle_nedelec_second_kind(degree,second,status)
 n=triangle_nedelec_second_kind_dof_count(second)
 end if
 call check_condition(status==0,'Public vector basis initialization')
 if(status/=0)error stop 'Vector initialization failed'
 allocate(v(2,n),v1(2,n),v2(2,n),vd(2,n),vb(2,n),c(n),c1(n),c2(n),cd(n),cb(n))
 do i=1,n
 vb(:,i)=[real(i,dp)*.003_dp,-real(i,dp)*.005_dp];cb(i)=real(i,dp)*.007_dp
 end do
 do sample=1,6
 select case(sample)
 case(1)
 x=[0.0_dp,0.0_dp];direction=[.2_dp,.3_dp]
 case(2)
 x=[.5_dp,0.0_dp];direction=[-.2_dp,.3_dp]
 case(3)
 x=[0.0_dp,.5_dp];direction=[.2_dp,-.3_dp]
 case(4)
 x=[.2_dp,.3_dp];direction=[.17_dp,-.11_dp]
 case(5)
 x=[1.0_dp,0.0_dp];direction=[-.2_dp,.1_dp]
 case(6)
 x=[0.0_dp,1.0_dp];direction=[.1_dp,-.2_dp]
 end select
 if(family==1)then
 call evaluate_triangle_nedelec_first_kind(first,x(1),x(2),v,c,status)
 call evaluate_triangle_nedelec_first_kind(first,x(1)+h*direction(1),x(2)+h*direction(2),v1,c1,plus_status)
 call evaluate_triangle_nedelec_first_kind(first,x(1)+2*h*direction(1),x(2)+2*h*direction(2),v2,c2,twice_status)
 call evaluate_triangle_nedelec_first_kind_jvp(first,x(1),x(2),direction(1),direction(2),vd,cd,status)
 else
 call evaluate_triangle_nedelec_second_kind(second,x(1),x(2),v,c,status)
 call evaluate_triangle_nedelec_second_kind(second,x(1)+h*direction(1),x(2)+h*direction(2),v1,c1,plus_status)
 call evaluate_triangle_nedelec_second_kind(second,x(1)+2*h*direction(1),x(2)+2*h*direction(2),v2,c2,twice_status)
 call evaluate_triangle_nedelec_second_kind_jvp(second,x(1),x(2),direction(1),direction(2),vd,cd,status)
 end if
 call check_condition(status==0.and.plus_status==0.and.twice_status==0, &
 'Axis and vertex inward tangent accepted')
 call check_condition(all(ieee_is_finite(v)).and.all(ieee_is_finite(c)).and. &
 all(ieee_is_finite(vd)).and.all(ieee_is_finite(cd)), 'Public vector jets finite at axes and vertices')
 error=maxval(abs(vd-(-3*v+4*v1-v2)/(2*h)));scale=max(1.0_dp,maxval(abs(vd)))
 call check_condition(error/scale<2e-6_dp,'Vector tangent matches independent inward difference')
 error=maxval(abs(cd-(-3*c+4*c1-c2)/(2*h)));scale=max(1.0_dp,maxval(abs(cd)))
 call check_condition(error/scale<2e-6_dp,'Curl tangent matches independent inward difference')
 if(family==1)then
 call evaluate_triangle_nedelec_first_kind_vjp(first,x(1),x(2),vb,cb,bar(1),bar(2),status)
 else
 call evaluate_triangle_nedelec_second_kind_vjp(second,x(1),x(2),vb,cb,bar(1),bar(2),status)
 end if
 lhs=sum(vd*vb)+sum(cd*cb);rhs=dot_product(direction,bar)
 call check_condition(status==0.and.abs(lhs-rhs)<2e-11_dp*max(1.0_dp,abs(lhs)), &
 'Axis and vertex signed adjoint identity')
 end do
 deallocate(v,v1,v2,vd,vb,c,c1,c2,cd,cb)
 end do
 end do
 call check_summary('Triangle vector boundary polynomial jets')
end program
