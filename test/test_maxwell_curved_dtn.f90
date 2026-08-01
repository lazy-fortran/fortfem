program test_maxwell_curved_dtn
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_maxwell_trace_to_flux, apply_maxwell_trace_to_flux_jvp, &
        apply_maxwell_trace_to_flux_map, apply_maxwell_trace_to_flux_vjp, &
        assemble_maxwell_trace_to_flux_map, &
        assemble_maxwell_torus_curved_dtn_rwg_3d, generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp) :: electric_form(2, 2), flux_form(2, 2), trace_mass(2, 2)
    complex(dp) :: electric_form_dot(2, 2), flux_form_dot(2, 2)
    complex(dp) :: trace_mass_dot(2, 2), trace(2), trace_dot(2)
    complex(dp) :: flux(2), flux_dot(2), flux_plus(2), flux_minus(2)
    complex(dp) :: trace_bar(2), electric_form_bar(2, 2)
    complex(dp) :: flux_form_bar(2, 2), trace_mass_bar(2, 2)
    complex(dp) :: flux_bar(2), trace_seed(2)
    complex(dp), allocatable :: generic_map(:, :)
    complex(dp), allocatable :: torus_map(:, :), torus_trace(:), torus_flux(:)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), vertices(:, :)
    real(dp) :: lhs, rhs, epsilon, finite_difference_error
    integer :: status, torus_size
    logical :: all_passed

    all_passed = .true.
    electric_form = reshape([ &
        cmplx(2.0_dp, 0.2_dp, dp), cmplx(0.25_dp, -0.1_dp, dp), &
        cmplx(0.5_dp, -0.3_dp, dp), cmplx(1.5_dp, 0.15_dp, dp)], [2, 2])
    flux_form = reshape([ &
        cmplx(1.0_dp, -0.2_dp, dp), cmplx(0.2_dp, 0.1_dp, dp), &
        cmplx(-0.3_dp, 0.05_dp, dp), cmplx(0.7_dp, 0.25_dp, dp)], [2, 2])
    trace_mass = reshape([ &
        cmplx(1.4_dp, 0.0_dp, dp), cmplx(0.1_dp, -0.03_dp, dp), &
        cmplx(0.1_dp, 0.03_dp, dp), cmplx(1.1_dp, 0.0_dp, dp)], [2, 2])
    trace = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.35_dp, 0.4_dp, dp)]
    electric_form_dot = reshape([ &
        cmplx(0.03_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
        cmplx(0.02_dp, 0.01_dp, dp), cmplx(-0.02_dp, 0.03_dp, dp)], [2, 2])
    flux_form_dot = reshape([ &
        cmplx(-0.02_dp, 0.01_dp, dp), cmplx(0.03_dp, -0.04_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(0.04_dp, 0.01_dp, dp)], [2, 2])
    trace_mass_dot = reshape([ &
        cmplx(0.02_dp, 0.0_dp, dp), cmplx(-0.01_dp, 0.01_dp, dp), &
        cmplx(-0.01_dp, -0.01_dp, dp), cmplx(0.03_dp, 0.0_dp, dp)], [2, 2])
    trace_dot = [cmplx(-0.2_dp, 0.1_dp, dp), cmplx(0.3_dp, -0.15_dp, dp)]
    call apply_maxwell_trace_to_flux( &
        electric_form, flux_form, trace_mass, trace, flux, status)
    call record_condition(status == 0, &
        "curved Maxwell trace-to-flux map accepts a nonsingular form")
    call assemble_maxwell_trace_to_flux_map( &
        electric_form, flux_form, trace_mass, generic_map, status)
    call record_condition(status == 0, &
        "curved Maxwell trace-to-flux map assembles its Galerkin matrix")
    if (status == 0) then
        call apply_maxwell_trace_to_flux_map( &
            generic_map, trace, flux_plus, status)
        call record_condition(status == 0 .and. &
            maxval(abs(flux_plus - flux)) < 2.0e-12_dp, &
            "assembled and matrix-free trace maps agree")
    end if
    call apply_maxwell_trace_to_flux_jvp( &
        electric_form, flux_form, trace_mass, trace, electric_form_dot, &
        flux_form_dot, trace_mass_dot, trace_dot, flux_dot, status)
    call record_condition(status == 0, &
        "curved Maxwell trace-to-flux map JVP succeeds")

    epsilon = 2.0e-6_dp
    call apply_maxwell_trace_to_flux( &
        electric_form + epsilon*electric_form_dot, &
        flux_form + epsilon*flux_form_dot, &
        trace_mass + epsilon*trace_mass_dot, trace + epsilon*trace_dot, &
        flux_plus, status)
    call apply_maxwell_trace_to_flux( &
        electric_form - epsilon*electric_form_dot, &
        flux_form - epsilon*flux_form_dot, &
        trace_mass - epsilon*trace_mass_dot, trace - epsilon*trace_dot, &
        flux_minus, status)
    finite_difference_error = maxval(abs( &
        flux_dot - (flux_plus - flux_minus)/(2.0_dp*epsilon)))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "curved Maxwell trace-to-flux JVP matches a complete difference")

    flux_bar = [cmplx(0.3_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp)]
    call apply_maxwell_trace_to_flux_vjp( &
        electric_form, flux_form, trace_mass, trace, flux_bar, trace_bar, &
        electric_form_bar, flux_form_bar, trace_mass_bar, status)
    call record_condition(status == 0, &
        "curved Maxwell trace-to-flux map VJP succeeds")
    trace_seed = [cmplx(0.11_dp, -0.07_dp, dp), cmplx(-0.09_dp, 0.13_dp, dp)]
    call apply_maxwell_trace_to_flux_jvp( &
        electric_form, flux_form, trace_mass, trace, 0.7_dp*electric_form_dot, &
        0.7_dp*flux_form_dot, 0.7_dp*trace_mass_dot, trace_seed, flux_dot, &
        status)
    lhs = real(sum(conjg(flux_bar)*flux_dot), dp)
    rhs = real(sum(conjg(electric_form_bar)* &
        (0.7_dp*electric_form_dot)) + &
        sum(conjg(flux_form_bar)*(0.7_dp*flux_form_dot)) + &
        sum(conjg(trace_mass_bar)*(0.7_dp*trace_mass_dot)) + &
        sum(conjg(trace_bar)*trace_seed), dp)
    call record_condition(abs(lhs - rhs) < 2.0e-10_dp, &
        "curved Maxwell trace-to-flux products satisfy the complex adjoint identity")

    call generate_torus_surface_mesh(2.0_dp, 0.6_dp, 3, 3, vertices, &
        triangles, parameters)
    call assemble_maxwell_torus_curved_dtn_rwg_3d( &
        vertices, triangles, parameters, 2.0_dp, 0.6_dp, 0.4_dp, 1.7_dp, &
        3, 3.0e-4_dp, 1, 0.12_dp, torus_map, status)
    torus_size = 0
    if (status == 0) torus_size = size(torus_map, 1)
    call record_condition(status == 0, &
        "exact-curved torus Maxwell DtN assembles")
    if (status == 0) then
        call record_condition(torus_size > 0 .and. &
            torus_size == size(torus_map, 2), &
            "exact-curved torus Maxwell DtN is a square RWG/BC map")
        allocate(torus_trace(torus_size), torus_flux(torus_size))
        torus_trace = cmplx(0.0_dp, 0.0_dp, dp)
        torus_trace(1) = cmplx(1.0_dp, -0.2_dp, dp)
        call apply_maxwell_trace_to_flux_map( &
            torus_map, torus_trace, torus_flux, status)
        call record_condition(status == 0 .and. &
            maxval(abs(torus_flux)) > tiny(1.0_dp), &
            "curved Maxwell DtN map action is nontrivial")
    end if
    call check_summary("Curved Maxwell trace-to-flux DtN")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_curved_dtn
