program test_spline_annulus_interface
    !! spline_interface_quadrature on the interface of layout_ogrid with an
    !! annulus (interface_fraction): the interface between the ruled blocks
    !! 2-5 and the periodic annulus (block 6) is found from the ruled side,
    !! once per Gauss point, with the annulus as the neighbour, and the
    !! neighbour parameters map to the same point (issue #278: no points were
    !! returned, because one side of the annulus borders four blocks).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL, LAYOUT_ANNULUS_INTERFACE
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build, SPLINE_C0
    use fortfem_spline_eval, only: spline_side_quadrature_t, spline_interface_quadrature, &
        spline_geometry_at
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: n = 3, np = 4
    type(block_layout_t) :: layout
    type(spline_space_options_t) :: options
    type(spline_space_t) :: space
    type(spline_side_quadrature_t) :: side
    real(dp) :: wall(2, 96), iface(2, 96), angle, jet(2, 6), gap
    integer :: i, status, failures

    failures = 0
    do i = 1, 96
        angle = 2*pi*(i - 1)/96
        wall(:, i) = [3.0_dp + 1.2_dp*cos(angle), 1.5_dp*sin(angle)]
        iface(:, i) = [3.1_dp + 0.6_dp*cos(angle + 0.2_dp*sin(angle)), 0.8_dp*sin(angle)]
    end do
    call layout_ogrid([3.0_dp, 0.0_dp], 1.2_dp, 0.3_dp, layout, status, interface_fraction=0.6_dp)
    call require(status == 0, 'layout')
    call layout_set_curve(layout, LAYOUT_WALL, wall, status)
    call require(status == 0, 'wall')
    call layout_set_curve(layout, LAYOUT_ANNULUS_INTERFACE, iface, status)
    call require(status == 0, 'interface')
    options%degree = 2
    options%regularity = SPLINE_C0
    allocate (options%knots(2, layout%block_count))
    options%knots(:, :)%element_count = n
    options%knots(1, 6)%element_count = 4*n
    call spline_space_build(layout, options, space, status)
    call require(status == 0, 'space')
    call spline_interface_quadrature(space, LAYOUT_ANNULUS_INTERFACE, np, side, status)
    call require(status == 0, 'interface quadrature status')
    print '(a,i5)', 'interface points', side%point_count
    call require(side%point_count == 4*n*np, 'one point per Gauss point of the interface')
    call require(all(side%block >= 2 .and. side%block <= 5), 'from the ruled blocks')
    call require(all(side%neighbour_block == 6), 'the annulus is the neighbour')
    gap = 0
    do i = 1, side%point_count
        call spline_geometry_at(space, side%neighbour_block(i), side%neighbour_xi(:, i), jet, &
                                status)
        gap = max(gap, norm2(jet(:, 1) - side%s(:, i)))
    end do
    print '(a,es10.2)', 'largest distance of the neighbour point', gap
    call require(gap < 1e-10_dp, 'the neighbour parameters give the same point')
    if (failures > 0) error stop 'annulus interface quadrature'
    print '(a)', 'PASS annulus interface quadrature'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) return
        print '(2a)', 'FAIL ', what
        failures = failures + 1
    end subroutine require

end program test_spline_annulus_interface
