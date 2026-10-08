program test_reference_basis_reproduction
    use fortfem_kinds, only: dp
    use basis_p1_2d_module, only: basis_p1_2d_t
    use basis_q1_quad_2d_module, only: q1_shape_functions, q1_shape_derivatives
    use fortfem_basis_1d, only: p1_basis, p1_basis_derivative
    use fortfem_basis_edge_2d, only: evaluate_edge_basis_2d, &
        evaluate_edge_basis_curl_2d, evaluate_edge_basis_div_2d
    use fortfem_basis_rt_2d, only: evaluate_rt_basis_2d, evaluate_rt_basis_div_2d
    use check, only: check_condition, check_summary
    implicit none
    type(basis_p1_2d_t) :: p1
    real(dp) :: q(4), qx(4), qy(4), node_x(4), node_y(4), g(2), value
    real(dp) :: vertices(2,3), v1(2,3), v2(2,3), tangent(2), normal(2)
    real(dp) :: divs(3), curls(3), exact, a, x,y,coefficients(4)
    integer :: i,j,k
    vertices(:,1)=[0.0_dp,0.0_dp]
    vertices(:,2)=[1.0_dp,0.0_dp]
    vertices(:,3)=[0.0_dp,1.0_dp]
    node_x=[-1.0_dp,1.0_dp,1.0_dp,-1.0_dp]
    node_y=[-1.0_dp,-1.0_dp,1.0_dp,1.0_dp]
    do k=0,2
        x=real(k,dp)/4; y=real(2-k,dp)/4
        value=0; g=0
        do i=1,3
            value=value+p1%eval(i,x,y)
            g=g+p1%grad(i,x,y)
        end do
        call check_condition(abs(value-1)<1e-14_dp,'P1 constant value')
        call check_condition(maxval(abs(g))<1e-14_dp,'P1 constant derivative')
        value=3*p1_basis(2,x)+2
        a=2*p1_basis_derivative(1,x)+5*p1_basis_derivative(2,x)
        call check_condition(abs(value-(3*x+2))<1e-14_dp,'1D affine value')
        call check_condition(abs(a-3)<1e-14_dp,'1D affine derivative')
        call q1_shape_functions(x,y,q)
        call q1_shape_derivatives(x,y,qx,qy)
        coefficients=node_x*node_y
        call check_condition(abs(sum(q)-1)<1e-14_dp,'Q1 constant value')
        call check_condition(abs(sum(q*coefficients)-x*y)<1e-14_dp, &
            'Q1 independent bilinear value')
        call check_condition(abs(sum(qx*coefficients)-y)<1e-14_dp, &
            'Q1 independent bilinear x derivative')
        call check_condition(abs(sum(qy*coefficients)-x)<1e-14_dp, &
            'Q1 independent bilinear y derivative')
    end do
    call evaluate_edge_basis_curl_2d(.2_dp,.3_dp,.5_dp,curls)
    call evaluate_edge_basis_div_2d(.2_dp,.3_dp,.5_dp,divs)
    call check_condition(maxval(abs(curls-2))<1e-14_dp,'Whitney Stokes curl')
    call check_condition(maxval(abs(divs))<1e-14_dp,'Whitney reference div')
    call evaluate_rt_basis_div_2d(.2_dp,.3_dp,.5_dp,divs)
    call check_condition(maxval(abs(divs-2))<1e-14_dp,'RT Gauss div')
    do i=1,3
        j=mod(i,3)+1
        tangent=vertices(:,j)-vertices(:,i)
        normal=[tangent(2),-tangent(1)]
        call evaluate_edge_basis_2d(vertices(1,i),vertices(2,i),.5_dp,v1)
        call evaluate_edge_basis_2d(vertices(1,j),vertices(2,j),.5_dp,v2)
        do k=1,3
            exact=0
            if(k==i) exact=1
            a=dot_product((v1(:,k)+v2(:,k))/2,tangent)
            call check_condition(abs(a-exact)<1e-14_dp,'Whitney edge circulation')
        end do
        call evaluate_rt_basis_2d(vertices(1,i),vertices(2,i),.5_dp,v1)
        call evaluate_rt_basis_2d(vertices(1,j),vertices(2,j),.5_dp,v2)
        do k=1,3
            exact=0
            if(k==i) exact=1
            a=dot_product((v1(:,k)+v2(:,k))/2,normal)
            call check_condition(abs(a-exact)<1e-14_dp,'RT normal flux')
        end do
    end do
    call check_summary('Generated reference element independent reproduction')
end program test_reference_basis_reproduction
