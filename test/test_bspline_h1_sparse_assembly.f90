program test_bspline_h1_sparse_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_bspline_h1_operator_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_x(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_y(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    type(csc_t) :: mass, stiffness
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: coefficients(:), control_points(:, :, :)
    real(dp), allocatable :: product(:), weights(:, :), x_points(:), y_points(:)
    real(dp) :: energy, integral
    integer :: ix, iy

    x_points = greville_abscissae(knots_x, 2)
    y_points = greville_abscissae(knots_y, 2)
    allocate(control_points(2, size(x_points), size(y_points)))
    allocate(weights(size(x_points), size(y_points)))
    weights = 1.0_dp
    do iy = 1, size(y_points)
        do ix = 1, size(x_points)
            control_points(:, ix, iy) = [x_points(ix), y_points(iy)]
        end do
    end do

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 4, mass, &
        sparse_status, stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    allocate(coefficients(size(x_points)*size(y_points)))
    coefficients = 1.0_dp
    product = csc_matvec(mass, coefficients)
    integral = dot_product(coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(integral - 1.0_dp) < 3.0e-13_dp, &
        "Sparse spline mass integrates a constant exactly")

    call assemble_bspline_h1_operator_csc( &
        knots_x, knots_y, 2, 2, control_points, weights, 4, stiffness, &
        sparse_status, stiffness_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    do iy = 1, size(y_points)
        coefficients(1 + (iy - 1)*size(x_points): &
            iy*size(x_points)) = x_points
    end do
    product = csc_matvec(stiffness, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 1.0_dp) < 4.0e-12_dp, &
        "Sparse isogeometric stiffness gives exact affine-field energy")
    call check_condition(abs(sum(product)) < 3.0e-13_dp, &
        "Spline stiffness annihilates constants in the weak form")

    call check_summary("Sparse isogeometric H1 assembly")

contains

    pure function greville_abscissae(knots, degree) result(points)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable :: points(:)
        integer :: basis

        allocate(points(size(knots) - degree - 1))
        do basis = 1, size(points)
            points(basis) = sum(knots(basis + 1:basis + degree))/real(degree, dp)
        end do
    end function greville_abscissae

end program test_bspline_h1_sparse_assembly
