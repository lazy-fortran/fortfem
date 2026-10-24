program test_bem_shared_edge_charts
    use check, only: check_condition, check_summary
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortfem_generated_bem_edge_chart_1, only: bem_edge_chart_1
    use fortfem_generated_bem_edge_chart_2, only: bem_edge_chart_2
    use fortfem_generated_bem_edge_chart_3, only: bem_edge_chart_3
    use fortfem_generated_bem_edge_chart_4, only: bem_edge_chart_4
    use fortfem_generated_bem_edge_chart_5, only: bem_edge_chart_5
    use fortfem_generated_bem_edge_chart_6, only: bem_edge_chart_6
    use fortfem_generated_bem_edge_chart_7, only: bem_edge_chart_7
    implicit none
    real(dp) :: rn(5),rw(5),en(4),ew(4),zn(4),zw(4),qn(33),qw(33)
    real(dp) :: moments(10,10),volumes(7),out(5),rho,e,z,q,weight,expected
    real(dp) :: h,eta,x,y,aa,bb,angular,lo,hi,zz,ray,limit,delta,t,u,ss,vv
    real(dp) :: expected_map(5),maximum_map_error,maximum_moment_error
    real(dp) :: Xphys(3),Yphys(3),physical(3),scaled(3),domain_exact
    logical :: passed
    integer :: powers(2,10),sector,chart,a,b,c,d,i,j,k,n,degree,work
    passed=.true.
    h=40; n=0
    do degree=0,3
        do i=0,degree
            n=n+1; powers(:,n)=[i,degree-i]
        end do
    end do
    call gauss_legendre_ab(5,0.0_dp,1.0_dp,rn,rw)
    call gauss_legendre_ab(4,0.0_dp,1.0_dp,en,ew)
    call gauss_legendre_ab(4,0.0_dp,1.0_dp,zn,zw)
    call gauss_legendre_ab(33,0.0_dp,1.0_dp,qn,qw)
    moments=0; volumes=0; maximum_map_error=0; work=0
    do chart=1,7
        do a=1,5
            rho=rn(a)
            do b=1,4
                e=en(b)
                do c=1,4
                    z=zn(c)
                    do d=1,33
                        q=qn(d)
                        select case(chart)
                        case(1)
                            call bem_edge_chart_1(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        case(2)
                            call bem_edge_chart_2(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        case(3)
                            call bem_edge_chart_3(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        case(4)
                            call bem_edge_chart_4(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        case(5)
                            call bem_edge_chart_5(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        case(6)
                            call bem_edge_chart_6(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        case(7)
                            call bem_edge_chart_7(rho,e,z,q,h,out(1),out(2), &
                                out(3),out(4),out(5))
                        end select
                        if (.not. all(ieee_is_finite(out))) &
                            error stop 'nonfinite generated chart'
                        if (minval(out)<=0 .or. max(out(1)+out(2), &
                                out(3)+out(4))>=1) error stop 'chart domain'
                        ! Independent scalar chart equations, together with
                        ! exact simplex moments, check the emitted consumer.
                        select case(chart)
                        case(1)
                            sector=1; aa=q; bb=z; angular=1
                        case(2,3)
                            sector=3; eta=q*asinh(h); x=sinh(eta)/h
                            if (chart==3) x=-x
                            aa=max(0.0_dp,-x)+(1-abs(x))*z
                            bb=1-aa-x
                            angular=(1-abs(x))*cosh(eta)*asinh(h)/h
                        case(4:7)
                            sector=2
                            if (chart>=6) sector=4
                            lo=-asinh(1/h); hi=log(h)
                            if (mod(chart,2)==1) then
                                lo=log(h); hi=asinh(h)
                            end if
                            eta=lo+(hi-lo)*q
                            zz=(1+h*sinh(eta))/(h*h+1)
                            limit=1/max(zz,1-zz)
                            ray=1+(limit-1)*z
                            x=ray*zz; y=ray*(1-zz)
                            aa=ray-1
                            bb=(1-y)/(2-ray)
                            if (sector==4) bb=(1-x)/(2-ray)
                            angular=ray*(limit-1)*h*cosh(eta)*(hi-lo)/(h*h+1)
                        end select
                        delta=rho*(1-aa); u=rho*aa; vv=rho*bb
                        if (sector==2 .or. sector==4) then
                            delta=delta*bb; vv=rho
                        end if
                        t=(1-rho)*e+delta; ss=(1-rho)*e
                        if (sector>=3) then
                            ss=t; t=(1-rho)*e; x=u; u=vv; vv=x
                        end if
                        expected_map=[t,u,ss,vv,rho*rho*(1-rho)*angular]
                        maximum_map_error=max(maximum_map_error, &
                            maxval(abs(out-expected_map)))
                        ! Exact native vertex orders: target [C,A,B],
                        ! source [A,D,B]. Compare two physical constructions.
                        Xphys=[h*(1-out(1)-out(2)),out(1),0.0_dp]
                        Yphys=[h*(1-out(3)),out(3)+out(4),0.0_dp]
                        physical=Xphys-Yphys
                        scaled=[h*(ss-t-u),t-ss-vv,0.0_dp]
                        maximum_map_error=max(maximum_map_error, &
                            maxval(abs(physical-scaled)))
                        weight=out(5)*rw(a)*ew(b)*zw(c)*qw(d)
                        volumes(chart)=volumes(chart)+weight
                        do i=1,10
                            do j=1,10
                                moments(i,j)=moments(i,j)+weight* &
                                    out(1)**powers(1,i)*out(2)**powers(2,i)* &
                                    out(3)**powers(1,j)*out(4)**powers(2,j)
                                work=work+1
                            end do
                        end do
                    end do
                end do
            end do
        end do
        domain_exact=1.0_dp/48
        if (chart==1) domain_exact=1.0_dp/12
        if (chart==2 .or. chart==3) domain_exact=1.0_dp/24
        call record(abs(volumes(chart)-domain_exact)<=1.e-12_dp, &
            'transformed subdomain volume')
    end do
    maximum_moment_error=0
    do i=1,10
        do j=1,10
            expected=factorial(powers(1,i))*factorial(powers(2,i))/ &
                factorial(sum(powers(:,i))+2)* &
                factorial(powers(1,j))*factorial(powers(2,j))/ &
                factorial(sum(powers(:,j))+2)
            maximum_moment_error=max(maximum_moment_error, &
                abs(moments(i,j)-expected)/(1+abs(expected)))
        end do
    end do
    write(*,*) '100_NATIVE_MOMENTS_ERROR',maximum_moment_error
    write(*,*) 'MAP_PHYSICAL_ERROR',maximum_map_error
    write(*,*) 'SUBDOMAIN_VOLUMES',volumes
    write(*,*) 'CONTROL_WORK',work
    call record(maximum_moment_error<=1.e-12_dp,'100 exact simplex moments')
    call record(maximum_map_error<=1.e-12_dp,'independent affine map')
    call check_summary('shared-edge emitted charts')
    if (.not. passed) error stop 1
contains
    subroutine record(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed=passed .and. ok
        call check_condition(ok,label)
    end subroutine
    pure real(dp) function factorial(i) result(value)
        integer, intent(in) :: i
        integer :: j
        value=1
        do j=2,i
            value=value*j
        end do
    end function
end program
