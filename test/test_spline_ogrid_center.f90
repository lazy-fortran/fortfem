program test_spline_ogrid_center
    !! The core of an O-grid stays centred on its center when the outer curve
    !! is set: the mean of the four inner corners (the image of the centre of
    !! the bilinear core block) is the center, also for a wall that is not
    !! centred on it (a magnetic axis shifted outward from the boundary
    !! centroid). The core keeps the shape of the outer corners, scaled by
    !! core_scale about their mean.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, layout_set_curve, &
        LAYOUT_WALL
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), tol = 1e-13_dp
    real(dp), parameter :: center(2) = [1.2_dp, 0.05_dp]
    type(block_layout_t) :: layout
    real(dp) :: samples(2, 128), mean_outer(2), mean_inner(2), shape_error, angle
    integer :: status, i, k, failures

    failures = 0

    ! a wall centred on the center leaves the corners where they were
    call layout_ogrid(center, 0.4_dp, 0.5_dp, layout, status)
    call require(status == 0, 'centred O-grid')
    do i = 1, size(samples, 2)
        angle = 2*pi*(i - 1)/size(samples, 2)
        samples(:, i) = center + 0.4_dp*[cos(angle), sin(angle)]
    end do
    call layout_set_curve(layout, LAYOUT_WALL, samples, status)
    call require(status == 0, 'centred wall')
    mean_inner = sum(layout%vertices(:, 1:4), 2)/4
    call require(norm2(mean_inner - center) < tol, 'core centred for a centred wall')

    ! an elongated, triangular wall whose centroid lies 0.08 inboard of the
    ! center and 0.03 below it
    call layout_ogrid(center, 0.4_dp, 0.5_dp, layout, status)
    do i = 1, size(samples, 2)
        angle = 2*pi*(i - 1)/size(samples, 2)
        samples(:, i) = center + [-0.08_dp, -0.03_dp] + &
                        [0.4_dp*cos(angle + 0.3_dp*sin(angle)), 0.6_dp*sin(angle)]
    end do
    call layout_set_curve(layout, LAYOUT_WALL, samples, status)
    call require(status == 0, 'shifted wall')
    mean_inner = sum(layout%vertices(:, 1:4), 2)/4
    mean_outer = sum(layout%vertices(:, 5:8), 2)/4
    print '(a,2es11.3)', 'core centre - center', mean_inner - center
    call require(norm2(mean_inner - center) < tol, 'core centred for a shifted wall')
    shape_error = 0
    do k = 1, 4
        shape_error = max(shape_error, norm2(layout%vertices(:, k) - center - &
            layout%core_scale*(layout%vertices(:, k + 4) - mean_outer)))
    end do
    call require(shape_error < tol, 'core similar to the outer corners')

    if (failures > 0) error stop 'O-grid core centre'
    print '(a)', 'PASS O-grid core stays on its center'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) return
        print '(2a)', 'FAIL ', what
        failures = failures + 1
    end subroutine require

end program test_spline_ogrid_center
