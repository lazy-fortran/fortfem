program test_tetra_modal_vector_definitions
    use check,only:check_condition,check_summary
    use fortfem_kinds,only:dp
    use fortfem_generated_tetra_modal_vector_identities,only:evaluate_tetra_modal_vector_identities
    use fortfem_generated_tetra_modal_vector_identities_jvp,only:evaluate_tetra_modal_vector_identities_jvp
    use,intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    real(dp)::point(3),direction(3),phi,g(3),hessian(3,3),pd,gd(3)
    real(dp)::cc(3,3),values(3,3),curls(3,3),vd(3,3),cd(3,3),vp(3,3),cp(3,3),vm(3,3),cm(3,3)
    real(dp)::plus(3),minus(3),difference(3,3),fieldp(3),fieldm(3),curl_expected(3),err
    real(dp),parameter::step=2e-6_dp
    integer::powers(3),sample,p,q,r,family,axis,status
    direction=[.17_dp,-.23_dp,.11_dp]
    do sample=1,4
        select case(sample)
        case(1)
            point=0.0_dp
        case(2)
            point=[0.0_dp,.5_dp,-.25_dp]
        case(3)
            point=[.25_dp,.125_dp,.5_dp]
        case(4)
            point=[-.5_dp,0.0_dp,.25_dp]
        end select
        do p=0,3
            do q=0,3
                do r=0,3
                    powers=[p,q,r]
                    call scalar(point,powers,phi,g,hessian)
                    call evaluate_tetra_modal_vector_identities(point(1),point(2),point(3),phi,g(1),g(2),g(3),cc,values,curls)
                    call check_condition(all(ieee_is_finite(cc)).and.all(ieee_is_finite(values)).and. &
                        all(ieee_is_finite(curls)),'Modal fields and curls finite on zero axes')
                    do family=1,6
                        do axis=1,3
                            plus=point;minus=point;plus(axis)=plus(axis)+step;minus(axis)=minus(axis)-step
                            fieldp=field(plus,powers,family);fieldm=field(minus,powers,family)
                            difference(:,axis)=(fieldp-fieldm)/(2*step)
                        end do
                        curl_expected=[difference(3,2)-difference(2,3),difference(1,3)-difference(3,1), &
                            difference(2,1)-difference(1,2)]
                        if(family<=3)then
                            err=maxval(abs(cc(:,family)-curl_expected))
                        else
                            err=maxval(abs(curls(:,family-3)-curl_expected))
                            call check_condition(maxval(abs(values(:,family-3)-field(point,powers,family)))<1e-13_dp, &
                                'Generated cross field matches independent polynomial definition')
                        end if
                        call check_condition(err<2e-8_dp,'Generated curl matches independent coordinate differences')
                    end do
                    pd=dot_product(g,direction);gd=matmul(hessian,direction)
                    call evaluate_tetra_modal_vector_identities_jvp(point(1),point(2),point(3),phi,g(1),g(2),g(3), &
                        direction(1),direction(2),direction(3),pd,gd(1),gd(2),gd(3),vd,cd)
                    call mapped(point+step*direction,powers,vp,cp)
                    call mapped(point-step*direction,powers,vm,cm)
                    call check_condition(maxval(abs(vd-(vp-vm)/(2*step)))<2e-8_dp.and. &
                        maxval(abs(cd-(cp-cm)/(2*step)))<2e-8_dp,'Native directional products match independent total differences')
                end do
            end do
        end do
    end do
    call check_summary('Tetra modal fields native definition oracle')
contains
    function field(x,powers,family)result(v)
        real(dp),intent(in)::x(3)
        integer,intent(in)::powers(3),family
        real(dp)::v(3),phi
        phi=product(x**powers);v=0.0_dp
        if(family<=3)then
            v(family)=phi
        else
            select case(family)
            case(4)
                v=[-x(2),x(1),0.0_dp]*phi
            case(5)
                v=[-x(3),0.0_dp,x(1)]*phi
            case(6)
                v=[0.0_dp,-x(3),x(2)]*phi
            end select
        end if
    end function
    subroutine scalar(x,powers,phi,g,h)
        real(dp),intent(in)::x(3)
        integer,intent(in)::powers(3)
        real(dp),intent(out)::phi,g(3),h(3,3)
        integer::i,j,orders(3)
        orders=0;phi=partial(x,powers,orders)
        do i=1,3
            orders=0;orders(i)=1;g(i)=partial(x,powers,orders)
            do j=1,3
                orders=0;orders(i)=orders(i)+1;orders(j)=orders(j)+1;h(i,j)=partial(x,powers,orders)
            end do
        end do
    end subroutine
    function partial(x,powers,orders)result(v)
        real(dp),intent(in)::x(3)
        integer,intent(in)::powers(3),orders(3)
        real(dp)::v
        integer::i,k
        v=0.0_dp
        if(any(powers<orders))return
        v=1.0_dp
        do i=1,3
            do k=0,orders(i)-1
                v=v*real(powers(i)-k,dp)
            end do
            v=v*x(i)**(powers(i)-orders(i))
        end do
    end function
    subroutine mapped(x,powers,v,c)
        real(dp),intent(in)::x(3)
        integer,intent(in)::powers(3)
        real(dp),intent(out)::v(3,3),c(3,3)
        real(dp)::phi,g(3),h(3,3),cc(3,3)
        call scalar(x,powers,phi,g,h)
        call evaluate_tetra_modal_vector_identities(x(1),x(2),x(3),phi,g(1),g(2),g(3),cc,v,c)
    end subroutine
end program
