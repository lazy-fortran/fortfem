program test_reference_p2_scalar_intervals
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval
    use fortfem_reference_scalar_intervals, only: generated_affine_p2_scalar_interval
    implicit none
    type(interval_t) :: value, gr, gz, det
    real(dp) :: vertices(2,3), nodes(6), reordered(6), shifted(2,3)
    real(dp), parameter :: gauge=1048576.0_dp
    real(dp) :: h, budget
    integer :: checks, j
    checks=0
    vertices=reshape([1.0_dp,0.0_dp,2.0_dp,0.0_dp,1.0_dp,1.0_dp],[2,3])
    ! Independent physical polynomial R^2+RZ+2Z^2+3R+5Z+7.
    nodes=[11.0_dp,17.0_dp,19.0_dp,13.75_dp,17.5_dp,14.5_dp]
    call evaluate(vertices,nodes,1.25_dp,0.25_dp)
    call contains(value,14.0_dp)
    call contains(gr,5.75_dp)
    call contains(gz,7.25_dp)
    call contains(det,1.0_dp)
    reordered=nodes([1,3,2,6,5,4])
    call evaluate(vertices(:,[1,3,2]),reordered,1.25_dp,0.25_dp)
    call contains(value,14.0_dp)
    call contains(gr,5.75_dp)
    call contains(gz,7.25_dp)
    call contains(det,-1.0_dp)
    shifted=vertices
    shifted(1,:)=shifted(1,:)+32
    shifted(2,:)=shifted(2,:)-16
    call evaluate(shifted,nodes,33.25_dp,-15.75_dp)
    call contains(value,14.0_dp)
    call contains(gr,5.75_dp)
    call contains(gz,7.25_dp)
    call evaluate(vertices,nodes+gauge,1.25_dp,0.25_dp)
    call contains(value,14.0_dp+gauge)
    call contains(gr,5.75_dp)
    call contains(gz,7.25_dp)
    ! Wrong edge ordering must not accidentally reproduce the physical oracle.
    reordered=nodes
    reordered([4,5])=nodes([5,4])
    call evaluate(vertices,reordered,1.25_dp,0.25_dp)
    if (value%lo <= 14 .and. value%hi >= 14) error stop 'wrong ordering unobserved'
    checks=checks+1
    ! Exact affine psi=R+2Z on shrinking dyadic cells. Polynomial grouping
    ! must avoid a fixed background-gradient interval floor at depth three.
    ! This machine-epsilon budget is for these inputs, not general conditioning.
    budget=128*epsilon(1.0_dp)*2
    do j=0,24,4
        h=2.0_dp**(-j)
        vertices=reshape([1.0_dp,0.0_dp,1+h,0.0_dp,1.0_dp,h],[2,3])
        nodes=[1.0_dp,1+h,1+2*h,1+h/2,1+1.5_dp*h,1+h]
        call evaluate_box(vertices,nodes,interval(1.0_dp,1+h/8),interval(0.0_dp,h/8))
        call contains(gr,1.0_dp)
        call contains(gz,2.0_dp)
        if (max(gr%hi-gr%lo,gz%hi-gz%lo)>budget) then
            print *, h,gr%hi-gr%lo,gz%hi-gz%lo,budget
            error stop 'affine cell-box enclosure retains a background-gradient floor'
        end if
        checks=checks+1
    end do
    print '(a,i0,a)', 'PASS: ',checks,' independent affine P2 interval controls'
contains
    subroutine evaluate(v,u,r,z)
        real(dp), intent(in) :: v(2,3), u(6), r,z
        call evaluate_box(v,u,interval(r),interval(z))
    end subroutine evaluate
    subroutine evaluate_box(v,u,r,z)
        real(dp), intent(in) :: v(2,3), u(6)
        type(interval_t), intent(in) :: r,z
        call generated_affine_p2_scalar_interval(r,z, &
            interval(v(1,1)),interval(v(2,1)),interval(v(1,2)),interval(v(2,2)), &
            interval(v(1,3)),interval(v(2,3)),interval(u(1)),interval(u(2)), &
            interval(u(3)),interval(u(4)),interval(u(5)),interval(u(6)),value,gr,gz,det)
    end subroutine evaluate_box
    subroutine contains(enclosure,expected)
        type(interval_t), intent(in) :: enclosure
        real(dp), intent(in) :: expected
        if (expected < enclosure%lo .or. expected > enclosure%hi) &
            error stop 'independent physical quadratic outside interval'
        checks=checks+1
    end subroutine contains
end program test_reference_p2_scalar_intervals
