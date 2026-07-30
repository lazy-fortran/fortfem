program test_bspline_jorek_harmonics
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_bspline_toroidal_poloidal_bracket, &
        apply_toroidal_fourier_derivative, &
        assemble_bspline_h1_weighted_mass_csc, &
        assemble_bspline_poloidal_bracket_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_r(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_z(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    integer, parameter :: modes(4) = [0, 1, 2, 3]
    complex(dp), allocatable :: advecting(:, :), derivative(:, :)
    complex(dp), allocatable :: residual(:, :), transported(:, :)
    real(dp), allocatable :: control_points(:, :, :), reference(:)
    real(dp), allocatable :: r_points(:), weights(:, :), z_points(:)
    type(csc_t) :: bracket, weighted_mass
    type(fortsparse_status_t) :: sparse_status
    integer :: dof, ix, iy, status

    r_points = greville_abscissae(knots_r, 2)
    z_points = greville_abscissae(knots_z, 2)
    allocate( &
        control_points(2, size(r_points), size(z_points)), &
        weights(size(r_points), size(z_points)))
    weights = 1.0_dp
    do iy = 1, size(z_points)
        do ix = 1, size(r_points)
            control_points(:, ix, iy) = [ &
                1.0_dp + r_points(ix), z_points(iy)]
        end do
    end do
    dof = size(r_points)*size(z_points)
    allocate(advecting(dof, 4), transported(dof, 4))
    advecting = cmplx(0.0_dp, 0.0_dp, dp)
    transported = cmplx(0.0_dp, 0.0_dp, dp)
    do iy = 1, size(z_points)
        do ix = 1, size(r_points)
            dof = ix + (iy - 1)*size(r_points)
            advecting(dof, 2) = cmplx(1.0_dp + r_points(ix), 0.0_dp, dp)
            transported(dof, 3) = cmplx(0.0_dp, z_points(iy), dp)
        end do
    end do

    call apply_bspline_toroidal_poloidal_bracket( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, advecting, &
        transported, 5, residual, sparse_status)
    call assemble_bspline_poloidal_bracket_csc( &
        knots_r, knots_z, 2, 2, control_points, weights, &
        real(advecting(:, 2), dp), 5, bracket, sparse_status)
    reference = csc_matvec(bracket, aimag(transported(:, 3)))
    call check_condition(sparse_status%code == 0 .and. maxval(abs( &
        residual(:, 4) - cmplx(0.0_dp, 1.0_dp, dp)*reference)) < 2.0e-14_dp, &
        "JOREK bracket convolution sends the (1,2) triad exactly to mode 3")
    call check_condition(maxval(abs(residual(:, :3))) < 2.0e-14_dp, &
        "JOREK bracket convolution leaves nonresonant modes empty")

    call apply_toroidal_fourier_derivative( &
        modes, transported, derivative, status)
    call check_condition(status == 0 .and. maxval(abs( &
        derivative(:, 3) + 2.0_dp*aimag(transported(:, 3)))) < 2.0e-14_dp, &
        "Toroidal derivative applies the exact i*n Fourier multiplier")

    call assemble_bspline_h1_weighted_mass_csc( &
        knots_r, knots_z, 2, 2, control_points, weights, &
        real(advecting(:, 2), dp), 5, weighted_mass, sparse_status)
    reference = csc_matvec(weighted_mass, [(1.0_dp, ix = 1, size(reference))])
    call check_condition(sparse_status%code == 0 .and. &
        abs(sum(reference) - 1.5_dp) < 3.0e-12_dp, &
        "Coefficient-weighted spline mass integrates the affine density")

    call check_summary("JOREK isogeometric toroidal harmonics")

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

end program test_bspline_jorek_harmonics
