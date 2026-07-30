program test_maxwell_pec_sphere_mie
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_maxwell_efie_far_field_rwg_3d, &
        generate_sphere_surface_mesh, solve_maxwell_pec_efie_rwg_3d, &
        solve_maxwell_pec_regularized_cfie_rwg_3d
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    complex(dp), allocatable :: density(:)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: cfie_cross_section, cfie_error
    real(dp) :: cross_sections(0:2), errors(0:2), exact_cross_section
    integer :: level, status
    logical :: all_passed

    all_passed = .true.
    exact_cross_section = mie_pec_cross_section(0.8_dp, 1.0_dp)
    do level = 0, 2
        call generate_sphere_surface_mesh(1.0_dp, level, vertices, triangles)
        call solve_maxwell_pec_efie_rwg_3d( &
            vertices, triangles, [0.0_dp, 0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
            cmplx(0.0_dp, 0.0_dp, dp)], 0.8_dp, 1.0_dp, 6, 1.0e-3_dp, 1, &
            density, status)
        cross_sections(level) = numerical_cross_section( &
            vertices, triangles, density, 0.8_dp, status)
        errors(level) = abs(cross_sections(level) - exact_cross_section)
        call record_condition(status == 0 .and. cross_sections(level) > 0.0_dp, &
            "PEC sphere EFIE produces positive radiated power")
    end do
    if (errors(1) >= errors(0) .or. errors(2) >= errors(1) .or. &
        errors(2)/exact_cross_section >= 0.55_dp) then
        write (*, '(A,ES14.5,A,2ES14.5,A,2ES14.5)') &
            "Mie cross section ", exact_cross_section, " numerical ", &
            cross_sections(1:2), " errors ", errors(1:2)
    end if
    call record_condition(errors(1) < errors(0) .and. errors(2) < errors(1), &
        "PEC sphere EFIE cross section converges toward the Mie series")
    call record_condition(errors(2)/exact_cross_section < 0.55_dp, &
        "Refined PEC sphere reaches the coarse-mesh Mie accuracy target")
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call solve_maxwell_pec_regularized_cfie_rwg_3d( &
        vertices, triangles, [0.0_dp, 0.0_dp, 1.0_dp], &
        [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)], 0.8_dp, 1.0_dp, 3, 1.0e-4_dp, 1, &
        density, status)
    cfie_cross_section = numerical_cross_section( &
        vertices, triangles, density, 0.8_dp, status)
    cfie_error = abs(cfie_cross_section - exact_cross_section)
    if (status /= 0 .or. cfie_cross_section <= 0.0_dp .or. &
        cfie_error >= exact_cross_section) then
        write (*, *) "CFIE/Mie cross sections", &
            cfie_cross_section, exact_cross_section
    end if
    call record_condition(status == 0 .and. cfie_cross_section > 0.0_dp .and. &
        cfie_error < exact_cross_section, &
        "regularized CFIE sphere solve agrees with Mie scattering on coarse mesh")
    call check_summary("PEC sphere Maxwell EFIE versus Mie series")
    if (.not. all_passed) error stop 1

contains

    function numerical_cross_section( &
            vertices, triangles, density, wave_number, status) result(value)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :)
        complex(dp), intent(in) :: density(:)
        integer, intent(out) :: status
        real(dp) :: value

        real(dp) :: azimuth, direction(3), radial, mu_nodes(8), mu_weights(8)
        complex(dp) :: far_field(3)
        integer :: azimuth_index, mu_index

        call gauss_legendre_ab( &
            8, -1.0_dp, 1.0_dp, mu_nodes, mu_weights)
        value = 0.0_dp
        do mu_index = 1, 8
            radial = sqrt(max(0.0_dp, 1.0_dp - mu_nodes(mu_index)**2))
            do azimuth_index = 0, 15
                azimuth = 2.0_dp*acos(-1.0_dp)* &
                    real(azimuth_index, dp)/16.0_dp
                direction = [ &
                    radial*cos(azimuth), radial*sin(azimuth), &
                    mu_nodes(mu_index)]
                call evaluate_maxwell_efie_far_field_rwg_3d( &
                    vertices, triangles, density, direction, wave_number, &
                    1.0_dp, 8, far_field, status)
                if (status /= 0) return
                value = value + mu_weights(mu_index)* &
                    2.0_dp*acos(-1.0_dp)/16.0_dp*sum(abs(far_field)**2)
            end do
        end do
    end function numerical_cross_section

    pure function mie_pec_cross_section(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        real(dp) :: value

        complex(dp) :: electric, hankel, hankel_minus
        complex(dp) :: magnetic, riccati_derivative
        real(dp) :: argument, bessel, bessel_minus, derivative
        integer :: degree, maximum_degree

        argument = wave_number*radius
        maximum_degree = max(12, int(argument + 4.0_dp*argument**(1.0_dp/3.0_dp) &
            + 2.0_dp))
        value = 0.0_dp
        do degree = 1, maximum_degree
            bessel = bessel_jn(degree, argument)
            bessel_minus = bessel_jn(degree - 1, argument)
            hankel = cmplx(bessel, bessel_yn(degree, argument), dp)
            hankel_minus = cmplx( &
                bessel_minus, bessel_yn(degree - 1, argument), dp)
            derivative = argument*bessel_minus - real(degree, dp)*bessel
            riccati_derivative = argument*hankel_minus - &
                real(degree, dp)*hankel
            electric = -derivative/riccati_derivative
            magnetic = -bessel/hankel
            value = value + real(2*degree + 1, dp)*( &
                abs(electric)**2 + abs(magnetic)**2)
        end do
        value = 2.0_dp*acos(-1.0_dp)*value/wave_number**2
    end function mie_pec_cross_section

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_pec_sphere_mie
