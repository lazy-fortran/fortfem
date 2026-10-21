program test_parametric_bem_estimator
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t, &
        parametric_bem_refinement_estimate
    implicit none
    type(bem_operators_t) :: fine, empty
    real(dp) :: coarse_flux(2), fine_flux(4), indicator(2), energy
    integer :: parent(4), status, i
    logical :: passed

    passed = .true.
    allocate(fine%v(4, 4), fine%area(4))
    fine%v = 0.0_dp
    do i = 1, 4
        fine%v(i, i) = real(i, dp)
    end do
    fine%area = 0.5_dp
    parent = [1, 1, 2, 2]
    coarse_flux = [2.0_dp, -1.0_dp]
    fine_flux = [3.0_dp, 1.0_dp, 1.0_dp, -3.0_dp]
    call estimate()
    call record(status == 0, 'two-level estimator accepts nested P0 data')
    call record(abs(energy-sqrt(31.0_dp)) < 1.e-12_dp, &
        'energy uses the fine single-layer quadratic form')
    call record(maxval(abs(indicator-[1.0_dp, 2.0_dp])) < 1.e-12_dp, &
        'parent marking indicators include area and parent length')
    fine_flux = coarse_flux(parent)
    call estimate()
    call record(status == 0 .and. energy == 0.0_dp, &
        'identical functions on two meshes have zero energy difference')
    call record(all(indicator == 0.0_dp), 'identical functions need no marking')
    fine_flux = coarse_flux(parent) + 1.0_dp
    call estimate()
    call record(status == 0 .and. energy > 1.0_dp, &
        'refined flux change is detected independently of equation residual')
    parent(1) = 0
    call estimate()
    call record(status /= 0, 'zero parent is rejected before indexing')
    parent = 1
    call estimate()
    call record(status /= 0, 'each coarse panel needs children')
    parent = [1, 1, 2, 2]
    fine%area(1) = -1.0_dp
    call estimate()
    call record(status /= 0, 'negative area is rejected')
    fine%area(1) = 0.5_dp
    fine_flux(1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call estimate()
    call record(status /= 0, 'nonfinite flux is rejected')
    fine_flux = coarse_flux(parent) + 1.0_dp
    fine%v(1, 2) = 1.0_dp
    call estimate()
    call record(status /= 0, 'nonsymmetric energy operator is rejected')
    fine%v(1, 2) = 0.0_dp
    fine%v = -fine%v
    call estimate()
    call record(status /= 0, 'negative energy is rejected without absolute value')
    call parametric_bem_refinement_estimate(empty, coarse_flux, fine_flux, &
        parent, indicator, energy, status)
    call record(status /= 0, 'unallocated operator is rejected')
    call check_summary('Parametric BEM refined-mesh estimator')
    if (.not. passed) error stop 1
contains
    subroutine estimate()
        call parametric_bem_refinement_estimate(fine, coarse_flux, fine_flux, &
            parent, indicator, energy, status)
    end subroutine
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
