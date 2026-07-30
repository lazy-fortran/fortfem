program test_bspline_3d_sparse_assembly
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_bspline_h1_operator_3d_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_x(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.35_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_y(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    type(csc_t) :: mass, stiffness
    type(fortsparse_status_t) :: sparse_status
    real(dp), allocatable :: coefficients(:), controls(:, :, :, :)
    real(dp), allocatable :: product(:), weights(:, :, :)
    real(dp), allocatable :: x_points(:), y_points(:), z_points(:)
    real(dp) :: energy
    integer :: ix, iy, iz

    x_points = greville_abscissae(knots_x, 2)
    y_points = greville_abscissae(knots_y, 2)
    z_points = greville_abscissae(knots_y, 2)
    allocate( &
        controls(3, size(x_points), size(y_points), size(z_points)), &
        weights(size(x_points), size(y_points), size(z_points)))
    weights = 1.0_dp
    do iz = 1, size(z_points)
        do iy = 1, size(y_points)
            do ix = 1, size(x_points)
                controls(:, ix, iy, iz) = [ &
                    1.0_dp + 2.0_dp*x_points(ix), &
                    -1.0_dp + 3.0_dp*y_points(iy), &
                    0.5_dp + 4.0_dp*z_points(iz)]
            end do
        end do
    end do

    call assemble_bspline_h1_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, mass, &
        sparse_status, stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    allocate(coefficients(size(x_points)*size(y_points)*size(z_points)))
    coefficients = 1.0_dp
    product = csc_matvec(mass, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(sparse_status%code == 0 .and. &
        abs(energy - 24.0_dp) < 2.0e-11_dp, &
        "Sparse 3D spline mass integrates affine-box volume exactly")

    call assemble_bspline_h1_operator_3d_csc( &
        knots_x, knots_y, knots_y, 2, 2, 2, controls, weights, 4, stiffness, &
        sparse_status, stiffness_coefficient=1.0_dp, mass_coefficient=0.0_dp)
    do iz = 1, size(z_points)
        do iy = 1, size(y_points)
            do ix = 1, size(x_points)
                coefficients(index_3d( &
                    ix, iy, iz, size(x_points), size(y_points))) = &
                    1.0_dp + 2.0_dp*x_points(ix)
            end do
        end do
    end do
    product = csc_matvec(stiffness, coefficients)
    energy = dot_product(coefficients, product)
    call check_condition(abs(energy - 24.0_dp) < 3.0e-10_dp, &
        "Sparse 3D spline diffusion gives exact coordinate-field energy")

    call check_summary("Sparse 3D isogeometric H1 assembly")

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

    pure integer function index_3d( &
            ix, iy, iz, count_x, count_y) result(index)
        integer, intent(in) :: ix, iy, iz, count_x, count_y

        index = ix + (iy - 1)*count_x + (iz - 1)*count_x*count_y
    end function index_3d

end program test_bspline_3d_sparse_assembly
