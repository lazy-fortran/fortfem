program test_bspline_jorek_thermodynamics
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_bspline_jorek_thermodynamic_jvp, &
        apply_bspline_jorek_thermodynamic_rhs
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: knots_r(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_z(7) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 1.0_dp, 1.0_dp, 1.0_dp]
    integer, parameter :: modes(3) = [0, 1, 2]
    real(dp), parameter :: gamma = 5.0_dp/3.0_dp
    real(dp), parameter :: finite_difference_step = 1.0e-6_dp
    complex(dp), allocatable :: field(:, :), field_direction(:, :)
    complex(dp), allocatable :: jvp(:, :), minus(:, :), plus(:, :)
    complex(dp), allocatable :: potential(:, :), potential_direction(:, :)
    complex(dp), allocatable :: rhs(:, :)
    real(dp), allocatable :: cylindrical_test(:)
    real(dp), allocatable :: control_points(:, :, :), r(:), weights(:, :), z(:)
    type(fortsparse_status_t) :: status
    integer :: dof, ix, iz

    r = greville_abscissae(knots_r, 2)
    z = greville_abscissae(knots_z, 2)
    allocate(control_points(2, size(r), size(z)), weights(size(r), size(z)))
    weights = 1.0_dp
    do iz = 1, size(z)
        do ix = 1, size(r)
            control_points(:, ix, iz) = [1.0_dp + r(ix), z(iz)]
        end do
    end do
    dof = size(r)*size(z)
    allocate( &
        field(dof, size(modes)), potential(dof, size(modes)), &
        field_direction(dof, size(modes)), &
        potential_direction(dof, size(modes)))
    field = cmplx(0.0_dp, 0.0_dp, dp)
    potential = cmplx(0.0_dp, 0.0_dp, dp)
    field(:, 1) = cmplx(1.0_dp, 0.0_dp, dp)
    do iz = 1, size(z)
        do ix = 1, size(r)
            dof = ix + (iz - 1)*size(r)
            potential(dof, 1) = cmplx(z(iz), 0.0_dp, dp)
        end do
    end do
    call apply_bspline_jorek_thermodynamic_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, field, &
        potential, gamma, 5, rhs, status)
    call check_condition(status%code == 0 .and. &
        abs(sum(rhs(:, 1)) - 2.0_dp*gamma) < 3.0e-12_dp, &
        "JOREK pressure residual exactly reproduces 2*gamma*p*dZ(u)")
    call check_condition(maxval(abs(rhs(:, 2:))) < 2.0e-14_dp, &
        "JOREK pressure residual does not populate nonresonant modes")

    allocate(cylindrical_test(size(field, 1)))
    field = cmplx(0.0_dp, 0.0_dp, dp)
    potential = cmplx(0.0_dp, 0.0_dp, dp)
    field(:, 1) = cmplx(1.0_dp, 0.0_dp, dp)
    do iz = 1, size(z)
        do ix = 1, size(r)
            dof = ix + (iz - 1)*size(r)
            cylindrical_test(dof) = 1.0_dp + r(ix)
            potential(dof, 1) = cmplx( &
                r(ix)*(1.0_dp - r(ix))*z(iz)*(1.0_dp - z(iz)), &
                0.0_dp, dp)
        end do
    end do
    call apply_bspline_jorek_thermodynamic_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, field, &
        potential, 1.0_dp, 5, rhs, status)
    call check_condition(abs(dot_product( &
        cylindrical_test, real(rhs(:, 1), dp))) < 3.0e-12_dp, &
        "JOREK density transport conserves cylindrical mass")

    call seed_random_numbers()
    call random_complex_field(field)
    call random_complex_field(potential)
    call random_complex_field(field_direction)
    call random_complex_field(potential_direction)
    call apply_bspline_jorek_thermodynamic_jvp( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, field, &
        potential, field_direction, potential_direction, gamma, 5, jvp, &
        status)
    call apply_bspline_jorek_thermodynamic_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, &
        field + finite_difference_step*field_direction, &
        potential + finite_difference_step*potential_direction, gamma, 5, &
        plus, status)
    call apply_bspline_jorek_thermodynamic_rhs( &
        knots_r, knots_z, 2, 2, control_points, weights, modes, &
        field - finite_difference_step*field_direction, &
        potential - finite_difference_step*potential_direction, gamma, 5, &
        minus, status)
    call check_condition(maxval(abs( &
        (plus - minus)/(2.0_dp*finite_difference_step) - jvp)) < 3.0e-8_dp, &
        "Analytical JOREK thermodynamic JVP matches a central difference")

    call check_summary("JOREK isogeometric thermodynamic transport")

contains

    subroutine random_complex_field(values)
        complex(dp), intent(out) :: values(:, :)
        real(dp) :: imaginary_part(size(values, 1), size(values, 2))
        real(dp) :: real_part(size(values, 1), size(values, 2))

        call random_number(real_part)
        call random_number(imaginary_part)
        values = cmplx(real_part, imaginary_part, dp)
    end subroutine random_complex_field

    subroutine seed_random_numbers()
        integer, allocatable :: seed(:)
        integer :: entry

        call random_seed(size=entry)
        allocate(seed(entry))
        seed = [(86028121 + 104395301*entry, entry = 1, size(seed))]
        call random_seed(put=seed)
    end subroutine seed_random_numbers

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

end program test_bspline_jorek_thermodynamics
