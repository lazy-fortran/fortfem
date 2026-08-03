program iga_jorek_flux
    use fortfem_feec, only: &
        assemble_bspline_h1_operator_csc, evaluate_bspline_basis
    use fortfem_time, only: advance_bspline_jorek_poloidal_flux_midpoint_steps
    use fortfem_kinds, only: dp
    use fortplot, only: &
        colorbar, figure, plot, pcolormesh, quiver, savefig, set_yscale, &
        title, xlabel, ylabel
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
    real(dp), allocatable :: initial_br(:, :), initial_bz(:, :)
    type(csc_t) :: mass
    type(fortsparse_status_t) :: status
    real(dp) :: elapsed, initial_norm, start_time, reversibility_error
    integer :: ix, iz, sample, unit

    call execute_command_line("mkdir -p "//output_directory)
    r = greville_abscissae(knots_r, degree)
    z = greville_abscissae(knots_z, degree)
    ! Keep the solution view substantially finer than the spline control
    ! mesh.  The physical field is reconstructed from the B-spline basis,
    ! so a 40-by-40 raster made the first gallery image look like a coarse
    ! source-cell diagnostic rather than a resolved solution.
    allocate(r_edges(81), z_edges(81))
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

    allocate( &
        initial_map(80, 80), final_map(80, 80), &
        initial_br(80, 80), initial_bz(80, 80))
    call sample_flux_field( &
        flux_history(:, 1), initial_map, initial_br, initial_bz)
    call sample_flux_map(flux_history(:, step_count + 1), final_map)
    call write_initial_field_csv()
    call render_plots()
    call write_benchmark()

contains

    subroutine render_plots()
        integer, parameter :: vector_stride = 8
        real(dp) :: arrow_r(10*10), arrow_z(10*10)
        real(dp) :: arrow_br(10*10), arrow_bz(10*10)
        integer :: arrow_count, ix, iz

        arrow_count = 0
        do ix = 1, size(initial_map, 2), vector_stride
            do iz = 1, size(initial_map, 1), vector_stride
                arrow_count = arrow_count + 1
                arrow_r(arrow_count) = 1.0_dp + &
                    (real(ix, dp) - 0.5_dp)/real(size(initial_map, 2), dp)
                arrow_z(arrow_count) = &
                    (real(iz, dp) - 0.5_dp)/real(size(initial_map, 1), dp)
                arrow_br(arrow_count) = initial_br(iz, ix)
                arrow_bz(arrow_count) = initial_bz(iz, ix)
            end do
        end do
        call figure(figsize=[8.0_dp, 5.5_dp])
        call pcolormesh(1.0_dp + r_edges, z_edges, initial_map, cmap="coolwarm")
        call quiver( &
            arrow_r(:arrow_count), arrow_z(:arrow_count), &
            arrow_br(:arrow_count), arrow_bz(:arrow_count), &
            scale=3.5_dp, scale_units="xy", angles="xy", color="black", &
            width=0.0025_dp)
        call colorbar(label="poloidal magnetic flux psi")
        call xlabel("major-radius coordinate R")
        call ylabel("vertical coordinate")
        call title("JOREK ideal flux subflow: flux solution and poloidal field")
        call savefig(output_directory//"/jorek_flux_initial_2d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(time, max(norm_drift, epsilon(1.0_dp)))
        call set_yscale("log")
        call xlabel("time")
        call ylabel("relative spline mass-norm drift")
        call title("JOREK ideal flux subflow: geometric invariant")
        call savefig(output_directory//"/jorek_flux_invariant_1d.png")

        call figure(figsize=[8.0_dp, 5.5_dp])
        call pcolormesh(1.0_dp + r_edges, z_edges, final_map, cmap="coolwarm")
        call colorbar(label="poloidal magnetic flux psi")
        call xlabel("major-radius coordinate R")
        call ylabel("vertical coordinate")
        call title("JOREK ideal flux subflow: final state")
        call savefig(output_directory//"/jorek_flux_final_2d.png")
    end subroutine render_plots

    subroutine sample_flux_field(coefficients, map, br, bz)
        real(dp), intent(in) :: coefficients(:)
        real(dp), intent(out) :: map(:, :), br(:, :), bz(:, :)

        real(dp), allocatable :: basis_r(:), basis_z(:)
        real(dp), allocatable :: derivative_r(:), derivative_z(:)
        real(dp) :: coordinate_r, coordinate_z, flux_value
        real(dp) :: flux_r, flux_z, major_radius
        integer :: basis_r_count, basis_z_count, ir, ix, iz, j, status

        basis_r_count = size(r)
        basis_z_count = size(z)
        if (size(coefficients) /= basis_r_count*basis_z_count) then
            error stop "JOREK field coefficient shape mismatch"
        end if
        do ix = 1, size(map, 2)
            coordinate_r = (real(ix, dp) - 0.5_dp)/real(size(map, 2), dp)
            call evaluate_bspline_basis( &
                knots_r, degree, coordinate_r, basis_r, derivative_r, status)
            if (status /= 0) error stop "JOREK radial field basis failed"
            do j = 1, size(map, 1)
                coordinate_z = (real(j, dp) - 0.5_dp)/real(size(map, 1), dp)
                call evaluate_bspline_basis( &
                    knots_z, degree, coordinate_z, basis_z, derivative_z, &
                    status)
                if (status /= 0) error stop "JOREK vertical field basis failed"
                flux_value = 0.0_dp
                flux_r = 0.0_dp
                flux_z = 0.0_dp
                do iz = 1, basis_z_count
                    do ir = 1, basis_r_count
                        flux_value = flux_value + basis_r(ir)*basis_z(iz)* &
                            coefficients(ir + (iz - 1)*basis_r_count)
                        flux_r = flux_r + derivative_r(ir)*basis_z(iz)* &
                            coefficients(ir + (iz - 1)*basis_r_count)
                        flux_z = flux_z + basis_r(ir)*derivative_z(iz)* &
                            coefficients(ir + (iz - 1)*basis_r_count)
                    end do
                end do
                major_radius = 1.0_dp + coordinate_r
                map(j, ix) = flux_value
                br(j, ix) = flux_z/major_radius
                bz(j, ix) = -flux_r/major_radius
            end do
        end do
    end subroutine sample_flux_field

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

    subroutine write_initial_field_csv()
        integer :: local_unit, local_status, ix, iz
        real(dp) :: magnitude

        open (newunit=local_unit, &
            file=output_directory//"/jorek_flux_initial_field.csv", &
            status="replace", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot open JOREK field CSV"
        write (local_unit, "(a)", iostat=local_status) &
            "R,Z,psi,Br,Bz,B_magnitude"
        if (local_status /= 0) error stop "cannot write JOREK field CSV header"
        do ix = 1, size(initial_map, 2)
            do iz = 1, size(initial_map, 1)
                magnitude = sqrt(initial_br(iz, ix)**2 + &
                    initial_bz(iz, ix)**2)
                write (local_unit, "(6(es24.16,','))", iostat=local_status) &
                    1.0_dp + (real(ix, dp) - 0.5_dp)/ &
                    real(size(initial_map, 2), dp), &
                    (real(iz, dp) - 0.5_dp)/real(size(initial_map, 1), dp), &
                    initial_map(iz, ix), initial_br(iz, ix), &
                    initial_bz(iz, ix), magnitude
                if (local_status /= 0) error stop "cannot write JOREK field CSV"
            end do
        end do
        close (local_unit, iostat=local_status)
        if (local_status /= 0) error stop "cannot close JOREK field CSV"
    end subroutine write_initial_field_csv

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
