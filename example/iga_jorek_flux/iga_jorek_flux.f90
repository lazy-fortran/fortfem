program iga_jorek_flux
    use fortfem_api, only: &
        advance_bspline_jorek_poloidal_flux_midpoint_steps, &
        assemble_bspline_h1_operator_csc, evaluate_bspline_basis
    use fortfem_kinds, only: dp
    use fortplot, only: &
        colorbar, figure, plot, pcolormesh, savefig, set_yscale, title, &
        xlabel, ylabel
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 2
    integer, parameter :: quadrature_order = 5
    integer, parameter :: step_count = 200
    real(dp), parameter :: time_step = 0.01_dp
    real(dp), parameter :: knots_r(9) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.2_dp, 0.5_dp, 0.8_dp, &
        1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: knots_z(8) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, 0.7_dp, &
        1.0_dp, 1.0_dp, 1.0_dp]
    character(*), parameter :: output_directory = &
        "output/example/iga_jorek_flux"
    real(dp), allocatable :: control_points(:, :, :), flux(:)
    real(dp), allocatable :: flux_history(:, :), flux_initial(:)
    real(dp), allocatable :: norm_drift(:), potential(:), r(:), z(:)
    real(dp), allocatable :: r_edges(:), z_edges(:), weights(:, :)
    real(dp), allocatable :: initial_map(:, :), final_map(:, :), time(:)
    type(csc_t) :: mass
    type(fortsparse_status_t) :: status
    real(dp) :: elapsed, initial_norm, start_time, reversibility_error
    integer :: ix, iz, sample, unit

    call execute_command_line("mkdir -p "//output_directory)
    r = greville_abscissae(knots_r, degree)
    z = greville_abscissae(knots_z, degree)
    allocate(r_edges(41), z_edges(41))
    do ix = 1, size(r_edges)
        r_edges(ix) = real(ix - 1, dp)/real(size(r_edges) - 1, dp)
        z_edges(ix) = r_edges(ix)
    end do
    allocate( &
        control_points(2, size(r), size(z)), &
        weights(size(r), size(z)), &
        flux(size(r)*size(z)), potential(size(r)*size(z)))
    weights = 1.0_dp
    do iz = 1, size(z)
        do ix = 1, size(r)
            control_points(:, ix, iz) = [1.0_dp + r(ix), z(iz)]
            sample = ix + (iz - 1)*size(r)
            flux(sample) = sin(acos(-1.0_dp)*r(ix))* &
                sin(2.0_dp*acos(-1.0_dp)*z(iz))
            potential(sample) = (1.0_dp + r(ix))*z(iz)
        end do
    end do
    flux_initial = flux

    call assemble_bspline_h1_operator_csc( &
        knots_r, knots_z, degree, degree, control_points, weights, &
        quadrature_order, mass, status, stiffness_coefficient=0.0_dp, &
        mass_coefficient=1.0_dp)
    if (status%code /= 0) error stop "JOREK gallery mass assembly failed"
    initial_norm = dot_product(flux, csc_matvec(mass, flux))

    call cpu_time(start_time)
    call advance_bspline_jorek_poloidal_flux_midpoint_steps( &
        knots_r, knots_z, degree, degree, control_points, weights, &
        potential, quadrature_order, time_step, step_count, flux, status, &
        flux_history)
    call cpu_time(elapsed)
    elapsed = elapsed - start_time
    if (status%code /= 0) error stop "JOREK gallery propagation failed"
    allocate(norm_drift(step_count + 1), time(step_count + 1))
    do sample = 1, step_count + 1
        time(sample) = real(sample - 1, dp)*time_step
        norm_drift(sample) = abs( &
            dot_product(flux_history(:, sample), &
            csc_matvec(mass, flux_history(:, sample)))/initial_norm - 1.0_dp)
    end do

    call advance_bspline_jorek_poloidal_flux_midpoint_steps( &
        knots_r, knots_z, degree, degree, control_points, weights, &
        potential, quadrature_order, -time_step, step_count, flux, status)
    if (status%code /= 0) error stop "JOREK gallery reverse solve failed"
    reversibility_error = maxval(abs(flux - flux_initial))
    if (maxval(norm_drift) > 5.0e-11_dp) &
        error stop "JOREK gallery mass invariant regression"
    if (reversibility_error > 5.0e-11_dp) &
        error stop "JOREK gallery reversibility regression"

    allocate(initial_map(40, 40), final_map(40, 40))
    call sample_flux_map(flux_history(:, 1), initial_map)
    call sample_flux_map(flux_history(:, step_count + 1), final_map)
    call render_plots()
    call write_benchmark()

contains

    subroutine render_plots()
        call figure(figsize=[8.0_dp, 5.5_dp])
        call pcolormesh(r_edges, z_edges, initial_map, cmap="coolwarm")
        call colorbar(label="poloidal magnetic flux psi")
        call xlabel("minor radial coordinate")
        call ylabel("vertical coordinate")
        call title("JOREK ideal flux subflow: initial state")
        call savefig(output_directory//"/jorek_flux_initial_2d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, max(norm_drift, epsilon(1.0_dp)))
        call set_yscale("log")
        call xlabel("time")
        call ylabel("relative spline mass-norm drift")
        call title("JOREK ideal flux subflow: geometric invariant")
        call savefig(output_directory//"/jorek_flux_invariant_1d.png")

        call figure(figsize=[8.0_dp, 5.5_dp])
        call pcolormesh(r_edges, z_edges, final_map, cmap="coolwarm")
        call colorbar(label="poloidal magnetic flux psi")
        call xlabel("minor radial coordinate")
        call ylabel("vertical coordinate")
        call title("JOREK ideal flux subflow: final state")
        call savefig(output_directory//"/jorek_flux_final_2d.png")
    end subroutine render_plots

    subroutine sample_flux_map(coefficients, map)
        real(dp), intent(in) :: coefficients(:)
        real(dp), intent(out) :: map(:, :)

        real(dp), allocatable :: basis_r(:), basis_z(:)
        real(dp), allocatable :: derivative_r(:), derivative_z(:)
        real(dp) :: coordinate_r, coordinate_z
        integer :: basis_r_count, basis_z_count, ir, ix, iz, j, status

        basis_r_count = size(r)
        basis_z_count = size(z)
        if (size(coefficients) /= basis_r_count*basis_z_count) then
            error stop "JOREK plot coefficient shape mismatch"
        end if
        do ix = 1, size(map, 2)
            coordinate_r = (real(ix, dp) - 0.5_dp)/real(size(map, 2), dp)
            call evaluate_bspline_basis( &
                knots_r, degree, coordinate_r, basis_r, derivative_r, status)
            if (status /= 0) error stop "JOREK radial plot basis failed"
            do j = 1, size(map, 1)
                coordinate_z = (real(j, dp) - 0.5_dp)/real(size(map, 1), dp)
                call evaluate_bspline_basis( &
                    knots_z, degree, coordinate_z, basis_z, derivative_z, &
                    status)
                if (status /= 0) error stop "JOREK vertical plot basis failed"
                map(j, ix) = 0.0_dp
                do iz = 1, basis_z_count
                    do ir = 1, size(r)
                        map(j, ix) = map(j, ix) + basis_r(ir) * &
                            basis_z(iz)*coefficients( &
                            ir + (iz - 1)*size(r))
                    end do
                end do
            end do
        end do
    end subroutine sample_flux_map

    subroutine write_benchmark()
        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a,i0)") "spline degrees of freedom: ", size(flux)
        write (unit, "(a,i0)") "implicit midpoint steps: ", step_count
        write (unit, "(a,es14.6)") "forward propagation seconds: ", elapsed
        write (unit, "(a,es14.6)") "steps per second: ", &
            real(step_count, dp)/max(elapsed, tiny(1.0_dp))
        write (unit, "(a,es14.6)") "maximum relative mass-norm drift: ", &
            maxval(norm_drift)
        write (unit, "(a,es14.6)") "forward-backward coefficient error: ", &
            reversibility_error
        close (unit)
    end subroutine write_benchmark

    pure function greville_abscissae(knots, polynomial_degree) result(points)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: polynomial_degree
        real(dp), allocatable :: points(:)
        integer :: basis

        allocate(points(size(knots) - polynomial_degree - 1))
        do basis = 1, size(points)
            points(basis) = sum( &
                knots(basis + 1:basis + polynomial_degree))/ &
                real(polynomial_degree, dp)
        end do
    end function greville_abscissae

end program iga_jorek_flux
