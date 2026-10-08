module fortfem_polynomial_candidate_jets
    use fortfem_kinds,only:dp
    implicit none
    private
    public::evaluate_monomial_jet2,evaluate_monomial_jet3,evaluate_integer_power
contains
    pure function evaluate_integer_power(coordinate, degree) result(power)
        real(dp), intent(in) :: coordinate
        integer, intent(in) :: degree
        real(dp) :: power, previous_power
        integer :: factor
        power=1.0_dp
        do factor=1,degree
            previous_power=power
            include '../generated/fortfem_power_multiply.inc'
        end do
    end function evaluate_integer_power
    pure subroutine power_jet(coordinate,degree,value,gradient,hessian)
        ! Zero histories keep degree-zero/one jets finite on coordinate axes.
        real(dp),intent(in)::coordinate
        integer,intent(in)::degree
        real(dp),intent(out)::value,gradient,hessian
        real(dp)::power,previous_power,second_previous_power
        integer::factor
        power=1.0_dp;previous_power=0.0_dp;second_previous_power=0.0_dp
        do factor=1,degree
            second_previous_power=previous_power
            previous_power=power
            include '../generated/fortfem_power_multiply.inc'
        end do
        include '../generated/fortfem_cached_power_jet.inc'
    end subroutine power_jet
    pure subroutine evaluate_monomial_jet2(point,powers,value,gradient,hessian)
        real(dp),intent(in)::point(2)
        integer,intent(in)::powers(2)
        real(dp),intent(out)::value,gradient(2)
        real(dp),intent(out),optional::hessian(2,2)
        real(dp)::factor_values(2),factor_gradients(2),factor_hessians(2)
        integer::coordinate
        do coordinate=1,2
            call power_jet(point(coordinate),powers(coordinate), &
                factor_values(coordinate),factor_gradients(coordinate),factor_hessians(coordinate))
        end do
        if(present(hessian))then
            include '../generated/fortfem_monomial_jet2_order2.inc'
        else
            include '../generated/fortfem_monomial_jet2_order1.inc'
        end if
    end subroutine evaluate_monomial_jet2
    pure subroutine evaluate_monomial_jet3(point,powers,value,gradient,hessian)
        real(dp),intent(in)::point(3)
        integer,intent(in)::powers(3)
        real(dp),intent(out)::value,gradient(3)
        real(dp),intent(out),optional::hessian(3,3)
        real(dp)::factor_values(3),factor_gradients(3),factor_hessians(3)
        integer::coordinate
        do coordinate=1,3
            call power_jet(point(coordinate),powers(coordinate), &
                factor_values(coordinate),factor_gradients(coordinate),factor_hessians(coordinate))
        end do
        if(present(hessian))then
            include '../generated/fortfem_monomial_jet3_order2.inc'
        else
            include '../generated/fortfem_monomial_jet3_order1.inc'
        end if
    end subroutine evaluate_monomial_jet3
end module fortfem_polynomial_candidate_jets
