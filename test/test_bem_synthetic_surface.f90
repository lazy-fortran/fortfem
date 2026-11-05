program test_bem_synthetic_surface
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_generated_bem_synthetic_surface, only: bem_synthetic_surface
    implicit none
    real(dp) :: theta,phi,x(3),dt(3),df(3),plus(3),minus(3),a(3),b(3)
    real(dp) :: rotated(3),angle,expected(3)
    real(dp), parameter :: step=2.e-6_dp
    integer :: i
    logical :: passed
    passed=.true.
    do i=1,4
        theta=.37_dp*i; phi=.23_dp*i
        call point(theta,phi,.1_dp,x,dt,df)
        call point(theta+step,phi,.1_dp,plus,a,b)
        call point(theta-step,phi,.1_dp,minus,a,b)
        call record(norm2(dt-(plus-minus)/(2*step))<1.e-8_dp, &
            'synthetic poloidal derivative matches independent differences')
        call point(theta,phi+step,.1_dp,plus,a,b)
        call point(theta,phi-step,.1_dp,minus,a,b)
        call record(norm2(df-(plus-minus)/(2*step))<1.e-8_dp, &
            'synthetic toroidal derivative matches independent differences')
    end do
    angle=2*acos(-1.0_dp)/5
    call point(theta,phi+angle,.1_dp,rotated,a,b)
    expected=[cos(angle)*x(1)-sin(angle)*x(2), &
        sin(angle)*x(1)+cos(angle)*x(2),x(3)]
    call record(norm2(rotated-expected)<1.e-12_dp, &
        'synthetic surface has exact five-sector rotational geometry')
    call point(theta,phi,0.0_dp,x,a,b)
    expected=[(5.5_dp+.55_dp*cos(theta))*cos(phi), &
        (5.5_dp+.55_dp*cos(theta))*sin(phi),.55_dp*sin(theta)]
    call record(norm2(x-expected)<1.e-12_dp,'zero ripple recovers the exact torus')
    call check_summary('Synthetic BEM scaling geometry')
    if (.not. passed) error stop 1
contains
    subroutine point(theta,phi,ripple,x,a,b)
        real(dp), intent(in) :: theta,phi,ripple
        real(dp), intent(out) :: x(3),a(3),b(3)
        call bem_synthetic_surface(theta,phi,5.5_dp,.55_dp,ripple,5.0_dp, &
            x(1),x(2),x(3),a(1),a(2),a(3),b(1),b(2),b(3))
    end subroutine
    subroutine record(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed=passed.and.ok
        call check_condition(ok,label)
    end subroutine
end program
