program test_tetra_whitney_reproduction
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_tetra_nedelec_first_order, only: evaluate_tetra_nedelec_first_order
    implicit none
    real(dp), parameter :: vertices(3, 4) = reshape( &
        [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
    integer, parameter :: edges(2, 6) = reshape( &
        [1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])
    integer, parameter :: faces(3, 4) = reshape( &
        [1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4], [3, 4])
    real(dp), parameter :: a(3) = [2.0_dp, -3.0_dp, 5.0_dp]
    real(dp), parameter :: b(3) = [-7.0_dp, 11.0_dp, 13.0_dp]
    real(dp) :: values(3, 6), curls(3, 6), point(3), tangent(3)
    real(dp) :: moments(6), exact(3), circulation(6), area_vector(3)
    real(dp) :: delta1(3), delta2(3), weights(4), reconstructed(3)
    integer :: edge, face, side, first, second, sample, status

    ! Exact oriented line moments of a + b cross r. Midpoints integrate
    ! every affine vector field exactly; no basis formula is used here.
    do edge = 1, 6
        first = edges(1, edge)
        second = edges(2, edge)
        point = (vertices(:, first) + vertices(:, second))/2.0_dp
        tangent = vertices(:, second) - vertices(:, first)
        exact = a + cross(b, point)
        moments(edge) = dot_product(exact, tangent)
    end do
    do sample = 1, 9
        weights = real([sample, 10 - sample, sample + 1, 11 - sample], dp)/22.0_dp
        point = matmul(vertices, weights)
        call evaluate_tetra_nedelec_first_order(point, values, curls, status)
        call check_condition(status == 0, "Whitney field accepts interior point")
        reconstructed = matmul(values, moments)
        exact = a + cross(b, point)
        call check_condition(maxval(abs(reconstructed - exact)) < 2.0e-13_dp, &
            "Oriented moments reproduce constant plus rigid rotation")
        call check_condition(maxval(abs(matmul(curls, moments) - 2.0_dp*b)) &
            < 2.0e-13_dp, "Reproduced rigid rotation has exact curl")
    end do

    ! Stokes theorem on all four oriented faces, separately for each basis.
    do face = 1, 4
        circulation = 0.0_dp
        do side = 1, 3
            first = faces(side, face)
            second = faces(mod(side, 3) + 1, face)
            point = (vertices(:, first) + vertices(:, second))/2.0_dp
            tangent = vertices(:, second) - vertices(:, first)
            call evaluate_tetra_nedelec_first_order(point, values, curls, status)
            circulation = circulation + matmul(transpose(values), tangent)
        end do
        delta1 = vertices(:, faces(2, face)) - vertices(:, faces(1, face))
        delta2 = vertices(:, faces(3, face)) - vertices(:, faces(1, face))
        area_vector = cross(delta1, delta2)/2.0_dp
        point = sum(vertices(:, faces(:, face)), dim=2)/3.0_dp
        call evaluate_tetra_nedelec_first_order(point, values, curls, status)
        do edge = 1, 6
            call check_condition(abs(circulation(edge) - &
                dot_product(curls(:, edge), area_vector)) < 2.0e-14_dp, &
                "Face circulation equals signed curl flux")
        end do
    end do
    call evaluate_tetra_nedelec_first_order([-0.1_dp, 0.2_dp, 0.3_dp], &
        values, curls, status)
    call check_condition(status /= 0, "Whitney wrapper rejects exterior point")
    call check_condition(maxval(abs(values)) == 0.0_dp .and. &
        maxval(abs(curls)) == 0.0_dp, "Exterior point returns zero arrays")
    call check_summary("Tetrahedral Whitney independent reproduction")
contains
    pure function cross(x, y) result(z)
        real(dp), intent(in) :: x(3), y(3)
        real(dp) :: z(3)
        z = [x(2)*y(3) - x(3)*y(2), x(3)*y(1) - x(1)*y(3), &
            x(1)*y(2) - x(2)*y(1)]
    end function cross
end program test_tetra_whitney_reproduction
