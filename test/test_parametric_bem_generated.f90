program test_parametric_bem_generated
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_bem_frame_jvp, only: bem_frame_jvp
    use fortfem_generated_bem_kernel, only: bem_kernel
    use fortfem_generated_bem_kernel_jvp, only: bem_kernel_jvp
    implicit none
    real(dp) :: x(6), dx(6), plus(4), minus(4), tangent(4), primal(4)
    real(dp), parameter :: h = 1.e-6_dp
    logical :: passed
    x = [1.0_dp,0.2_dp,0.1_dp,0.0_dp,1.1_dp,0.3_dp]
    dx = [0.1_dp,0.05_dp,-0.1_dp,0.2_dp,-0.1_dp,0.0_dp]
    call frame(x, primal)
    call bem_frame_jvp(x(1),x(2),x(3),x(4),x(5),x(6), &
        dx(1),dx(2),dx(3),dx(4),dx(5),dx(6), &
        tangent(1),tangent(2),tangent(3),tangent(4))
    call frame(x+h*dx, plus)
    call frame(x-h*dx, minus)
    passed = norm2(tangent-(plus-minus)/(2*h)) < 1.e-8_dp
    call check_condition(passed, 'generated frame JVP matches finite differences')
    call check_condition(abs(norm2(primal(2:))-1.0_dp) < 1.e-14_dp, &
        'generated normal is unit length')
    passed = passed .and. abs(norm2(primal(2:))-1.0_dp) < 1.e-14_dp
    call bem_kernel_jvp(x(1),x(2),x(3),x(4),x(5),x(6), &
        dx(1),dx(2),dx(3),dx(4),dx(5),dx(6),tangent(1),tangent(2))
    call kernel(x+h*dx, plus(:2))
    call kernel(x-h*dx, minus(:2))
    call check_condition(norm2(tangent(:2)-(plus(:2)-minus(:2))/(2*h)) < 1.e-8_dp, &
        'generated Green/double-layer JVP matches finite differences')
    passed = passed .and. &
        norm2(tangent(:2)-(plus(:2)-minus(:2))/(2*h)) < 1.e-8_dp
    call kernel(x, primal(:2))
    call check_condition(abs(primal(1)-1/norm2(x(:3))) < 1.e-14_dp, &
        'Green kernel matches independent inverse distance')
    passed = passed .and. abs(primal(1)-1/norm2(x(:3))) < 1.e-14_dp
    call check_summary('Generated parametric BEM leaves')
    if (.not. passed) error stop 1
contains
    subroutine frame(z, output)
        real(dp), intent(in) :: z(6)
        real(dp), intent(out) :: output(4)
        call bem_frame(z(1),z(2),z(3),z(4),z(5),z(6), &
            output(1),output(2),output(3),output(4))
    end subroutine
    subroutine kernel(z, output)
        real(dp), intent(in) :: z(6)
        real(dp), intent(out) :: output(2)
        call bem_kernel(z(1),z(2),z(3),z(4),z(5),z(6), output(1),output(2))
    end subroutine
end program
