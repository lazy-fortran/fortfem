program test_polynomial_candidate_jets
    use check,only:check_condition,check_summary
    use fortfem_kinds,only:dp
    use fortfem_polynomial_candidate_jets,only:evaluate_monomial_jet2,evaluate_monomial_jet3
    use, intrinsic::ieee_arithmetic,only:ieee_is_finite
    implicit none
    real(dp)::point(3),value,g(3),h(3,3),expected,eg(3),eh(3,3)
    real(dp)::v2,g2(2),h2(2,2),vfirst,gfirst(3)
    integer::powers(3),orders(3),sample,p,q,r,i,j
    do sample=1,6
        select case(sample)
        case(1)
            point=0.0_dp
        case(2)
            point=[0.0_dp,0.5_dp,-0.25_dp]
        case(3)
            point=[1.0_dp,0.0_dp,0.5_dp]
        case(4)
            point=[-0.5_dp,0.25_dp,0.0_dp]
        case(5)
            point=[0.25_dp,0.5_dp,0.125_dp]
        case(6)
            point=[-1.0_dp,1.0_dp,-1.0_dp]
        end select
        do p=0,12
            do q=0,12
                powers=[p,q,0]
                call evaluate_monomial_jet2(point(1:2),powers(1:2),v2,g2,h2)
                call oracle(point,powers,expected,eg,eh)
                call check_condition(ieee_is_finite(v2).and.all(ieee_is_finite(g2)).and. &
                    all(ieee_is_finite(h2)).and.abs(v2-expected)<1e-12_dp.and. &
                    maxval(abs(g2-eg(1:2)))<1e-12_dp.and. &
                    maxval(abs(h2-eh(1:2,1:2)))<1e-12_dp,'2D analytic monomial jet, degree0-12 per axis')
                call evaluate_monomial_jet2(point(1:2),powers(1:2),vfirst,gfirst(1:2))
                call check_condition(abs(vfirst-v2)<1e-12_dp.and. &
                    maxval(abs(gfirst(1:2)-g2))<1e-12_dp,'2D first-order analytic monomial jet')
            end do
        end do
        do p=0,4
            do q=0,4
                do r=0,4
                    powers=[p,q,r]
                    call evaluate_monomial_jet3(point,powers,value,g,h)
                    call oracle(point,powers,expected,eg,eh)
                    call check_condition(ieee_is_finite(value).and.all(ieee_is_finite(g)).and. &
                        all(ieee_is_finite(h)).and.abs(value-expected)<1e-12_dp.and. &
                        maxval(abs(g-eg))<1e-12_dp.and.maxval(abs(h-eh))<1e-12_dp, &
                        '3D analytic monomial jet at interior, signed, zero-axis points')
                    call evaluate_monomial_jet3(point,powers,vfirst,gfirst)
                    call check_condition(abs(vfirst-value)<1e-12_dp.and. &
                        maxval(abs(gfirst-g))<1e-12_dp,'3D first-order analytic monomial jet')
                end do
            end do
        end do
    end do
    call check_summary('Independent analytic polynomial jets')
contains
    subroutine oracle(point,powers,value,gradient,hessian)
        real(dp),intent(in)::point(3)
        integer,intent(in)::powers(3)
        real(dp),intent(out)::value,gradient(3),hessian(3,3)
        integer::i,j,orders(3)
        orders=0;value=partial(point,powers,orders)
        do i=1,3
            orders=0;orders(i)=1;gradient(i)=partial(point,powers,orders)
            do j=1,3
                orders=0;orders(i)=orders(i)+1;orders(j)=orders(j)+1
                hessian(i,j)=partial(point,powers,orders)
            end do
        end do
    end subroutine oracle
    function partial(point,powers,orders)result(value)
        real(dp),intent(in)::point(3)
        integer,intent(in)::powers(3),orders(3)
        real(dp)::value
        integer::i,k
        value=0.0_dp
        if(any(powers<orders))return
        value=1.0_dp
        do i=1,3
            do k=0,orders(i)-1
                value=value*real(powers(i)-k,dp)
            end do
            value=value*point(i)**(powers(i)-orders(i))
        end do
    end function partial
end program test_polynomial_candidate_jets
