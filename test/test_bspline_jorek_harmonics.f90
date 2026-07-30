program test_bspline_jorek_harmonics
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_bspline_jorek_flux_rhs, &
        apply_bspline_jorek_flux_jvp, &
        advance_bspline_jorek_poloidal_flux_midpoint_steps, &
        apply_bspline_toroidal_poloidal_bracket, &
        apply_toroidal_fourier_derivative, &
        assemble_bspline_h1_weighted_mass_csc, &
        assemble_bspline_h1_operator_csc, &
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
    complex(dp), allocatable :: flux_rhs(:, :)
    complex(dp), allocatable :: flux_jvp(:, :), flux_minus(:, :), flux_plus(:, :)
    complex(dp), allocatable :: current_direction(:, :), flux_direction(:, :)
    complex(dp), allocatable :: potential_direction(:, :)
    real(dp), allocatable :: control_points(:, :, :), reference(:)
    real(dp), allocatable :: flux_history(:, :)
    real(dp), allocatable :: flux_initial(:), flux_state(:), potential(:)
    real(dp), allocatable :: r_points(:), weights(:, :), z_points(:)
    real(dp), parameter :: finite_difference_step = 1.0e-6_dp
    type(csc_t) :: bracket, mass, weighted_mass
    type(fortsparse_status_t) :: sparse_status
    integer :: dof, ix, iy, status
    real(dp) :: initial_norm

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

    advecting = cmplx(0.0_dp, 0.0_dp, dp)
    transported = cmplx(0.0_dp, 0.0_dp, dp)
    residual = cmplx(0.0_dp, 0.0_dp, dp)
    transported(:, 2) = cmplx(1.0_dp, 0.0_dp, dp)
    advecting(:, 3) = cmplx(1.0_dp, 0.0_dp, dp)
    residual(:, 1) = cmplx(1.0_dp, 0.0_dp, dp)
    call apply_bspline_jorek_flux_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, advecting, &
        transported, residual, 3.0_dp, 2.0_dp, 5, flux_rhs, sparse_status)
    call check_condition(abs(sum(flux_rhs(:, 1)) - 3.0_dp) < 3.0e-12_dp, &
        "JOREK flux residual applies the resistive current source")
    call check_condition(abs(sum(flux_rhs(:, 2)) + &
        cmplx(0.0_dp, 2.0_dp, dp)) < 3.0e-12_dp, &
        "JOREK flux residual applies the exact -F0*dphi(u) coupling")
    call check_condition(abs(sum(flux_rhs(:, 3)) + 6.0_dp) < 2.0e-10_dp, &
        "JOREK flux residual applies eta*R^-2*dphi_phi(psi)")

    advecting = cmplx(0.0_dp, 0.0_dp, dp)
    transported = cmplx(0.0_dp, 0.0_dp, dp)
    residual = cmplx(0.0_dp, 0.0_dp, dp)
    do iy = 1, size(z_points)
        do ix = 1, size(r_points)
            dof = ix + (iy - 1)*size(r_points)
            advecting(dof, 1) = cmplx(1.0_dp + r_points(ix), 0.0_dp, dp)
            transported(dof, 1) = cmplx(z_points(iy), 0.0_dp, dp)
        end do
    end do
    call apply_bspline_jorek_flux_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, advecting, &
        transported, residual, 0.0_dp, 0.0_dp, 5, flux_rhs, sparse_status)
    call check_condition(abs(sum(flux_rhs(:, 1)) - 0.75_dp) < 3.0e-12_dp, &
        "JOREK flux residual reproduces the cylindrical R*[R,Z] weak action")

    allocate( &
        flux_direction(size(advecting, 1), size(modes)), &
        potential_direction(size(advecting, 1), size(modes)), &
        current_direction(size(advecting, 1), size(modes)))
    do iy = 1, size(modes)
        do ix = 1, size(advecting, 1)
            advecting(ix, iy) = cmplx( &
                sin(real(3*ix + iy, dp)), cos(real(ix + 2*iy, dp)), dp)
            transported(ix, iy) = cmplx( &
                cos(real(2*ix - iy, dp)), sin(real(ix + iy, dp)), dp)
            residual(ix, iy) = cmplx( &
                sin(real(ix - 3*iy, dp)), cos(real(2*ix + iy, dp)), dp)
            flux_direction(ix, iy) = cmplx( &
                cos(real(ix + 4*iy, dp)), sin(real(2*ix - iy, dp)), dp)
            potential_direction(ix, iy) = cmplx( &
                sin(real(2*ix + 3*iy, dp)), cos(real(ix - iy, dp)), dp)
            current_direction(ix, iy) = cmplx( &
                cos(real(3*ix - 2*iy, dp)), sin(real(ix + 2*iy, dp)), dp)
        end do
    end do
    call apply_bspline_jorek_flux_jvp( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, advecting, &
        transported, flux_direction, potential_direction, current_direction, &
        0.07_dp, 1.3_dp, 5, flux_jvp, sparse_status)
    call apply_bspline_jorek_flux_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, &
        advecting + finite_difference_step*flux_direction, &
        transported + finite_difference_step*potential_direction, &
        residual + finite_difference_step*current_direction, 0.07_dp, &
        1.3_dp, 5, flux_plus, sparse_status)
    call apply_bspline_jorek_flux_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, &
        advecting - finite_difference_step*flux_direction, &
        transported - finite_difference_step*potential_direction, &
        residual - finite_difference_step*current_direction, 0.07_dp, &
        1.3_dp, 5, flux_minus, sparse_status)
    call check_condition(maxval(abs( &
        (flux_plus - flux_minus)/(2.0_dp*finite_difference_step) - &
        flux_jvp)) < 2.0e-8_dp, &
        "Analytical JOREK flux JVP matches an independent central difference")

    allocate(flux_initial(size(reference)), flux_state(size(reference)))
    allocate(potential(size(reference)))
    do iy = 1, size(z_points)
        do ix = 1, size(r_points)
            dof = ix + (iy - 1)*size(r_points)
            flux_initial(dof) = sin(real(2*ix + 3*iy, dp))
            potential(dof) = &
                (1.0_dp + r_points(ix))*z_points(iy)
        end do
    end do
    call assemble_bspline_h1_operator_csc( &
        knots_r, knots_z, 2, 2, control_points, weights, 5, mass, &
        sparse_status, stiffness_coefficient=0.0_dp, mass_coefficient=1.0_dp)
    initial_norm = dot_product(flux_initial, csc_matvec(mass, flux_initial))
    flux_state = flux_initial
    call advance_bspline_jorek_poloidal_flux_midpoint_steps( &
        knots_r, knots_z, 2, 2, control_points, weights, potential, 5, &
        0.02_dp, 100, flux_state, sparse_status, flux_history)
    call check_condition( &
        size(flux_history, 2) == 101 .and. &
        maxval(abs(flux_history(:, 1) - flux_initial)) < 2.0e-14_dp .and. &
        maxval(abs(flux_history(:, 101) - flux_state)) < 2.0e-14_dp, &
        "JOREK midpoint trajectory contains both exact endpoints")
    call check_condition(abs(dot_product( &
        flux_state, csc_matvec(mass, flux_state)) - initial_norm) < 2.0e-11_dp, &
        "JOREK midpoint bracket propagator preserves the spatial mass norm")
    call advance_bspline_jorek_poloidal_flux_midpoint_steps( &
        knots_r, knots_z, 2, 2, control_points, weights, potential, 5, &
        -0.02_dp, 100, flux_state, sparse_status)
    call check_condition(maxval(abs(flux_state - flux_initial)) < 2.0e-11_dp, &
        "JOREK midpoint bracket propagator is exactly time reversible")

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
