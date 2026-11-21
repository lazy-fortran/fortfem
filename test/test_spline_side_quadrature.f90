program test_spline_side_quadrature
    !! Boundary line integrals and point location on the O-grid disc: the
    !! counter-clockwise boundary integral of (s1 ds2 - s2 ds1)/2 equals the
    !! area from the cell quadrature (Green), the wall length converges to
    !! 2 pi, and spline_locate recovers block and parameters of quadrature
    !! points.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_spline_layout, only: block_layout_t, layout_ogrid, LAYOUT_WALL
    use fortfem_spline_space, only: spline_space_t, spline_space_options_t, &
        spline_space_build
    use fortfem_spline_eval, only: spline_quadrature_t, spline_quadrature_build, &
        spline_side_quadrature_t, spline_side_quadrature, spline_locate
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    type(block_layout_t) :: layout
    type(spline_space_t) :: space
    type(spline_quadrature_t) :: quad
    type(spline_side_quadrature_t) :: side
    real(dp) :: green, length(2), xi(2), worst
    integer :: status, q, block, level, misses

    call layout_ogrid([0.3_dp, -0.2_dp], 1.0_dp, 0.5_dp, layout, status)
    do level = 1, 2
        call build(4*level, space, status)
        call spline_side_quadrature(space, LAYOUT_WALL, 4, side, status)
        length(level) = sum(side%weight)
    end do
    call spline_quadrature_build(space, 4, quad, status)
    green = 0.0_dp
    do q = 1, side%point_count
        green = green + 0.5_dp*side%weight(q)*((side%s(1, q) - 0.3_dp)* &
            side%tangent(2, q) - (side%s(2, q) + 0.2_dp)*side%tangent(1, q))
    end do
    call check_condition(status == 0 .and. side%point_count == 4*4*8 .and. &
        abs(green - sum(quad%weight)) < 1.0e-13_dp .and. &
        maxval(abs(norm2(side%tangent, dim=1) - 1.0_dp)) < 1.0e-14_dp, &
        "counter-clockwise wall integral reproduces the area (Green)")
    call check_condition(abs(length(2) - 2*pi) < abs(length(1) - 2*pi)/12.0_dp, &
        "wall length converges to 2 pi at order p + 1")

    misses = 0
    worst = 0.0_dp
    do q = 1, quad%point_count, 7
        call spline_locate(space, 1, quad%s(:, q), block, xi, status)
        if (status /= 0 .or. block /= quad%block(q)) then
            misses = misses + 1
        else
            worst = max(worst, maxval(abs(xi - quad%xi(:, q))))
        end if
    end do
    call check_condition(misses == 0 .and. worst < 1.0e-12_dp, &
        "spline_locate recovers block and parameters of quadrature points")
    call spline_locate(space, 1, [5.0_dp, 5.0_dp], block, xi, status)
    call check_condition(status /= 0, "points outside the layout are rejected")
    call check_summary("spline_side_quadrature")

contains

    subroutine build(cells, space, status)
        integer, intent(in) :: cells
        type(spline_space_t), intent(out) :: space
        integer, intent(out) :: status

        type(spline_space_options_t) :: options

        options%degree = 3
        allocate(options%knots(2, 5))
        options%knots(:, :)%element_count = cells
        call spline_space_build(layout, options, space, status)
    end subroutine build

end program test_spline_side_quadrature
