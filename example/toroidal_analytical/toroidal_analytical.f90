program toroidal_analytical
    use fortfem_api, only: &
        cartesian_to_toroidal, &
        evaluate_helmholtz_representation_torus_curved_3d, &
        evaluate_laplace_representation_torus_curved_3d, &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p, &
        generate_torus_surface_mesh, &
        solve_helmholtz_bem_dtn_torus_curved_3d, &
        solve_laplace_bem_dtn_torus_curved_3d, toroidal_point_to_cartesian
    use fortfem_kinds, only: dp
    use fortplot, only: figure, plot, pcolormesh, add_scatter, &
        xlabel, ylabel, title, legend, savefig
    implicit none

    integer, parameter :: trace_points = 61
    integer, parameter :: theta_cells = 24
    integer, parameter :: phi_cells = 16
    integer, parameter :: surface_points = theta_cells*phi_cells
    integer, parameter :: degree_index = 2
    integer, parameter :: order = 1
    real(dp), parameter :: scale = 1.0_dp
    real(dp), parameter :: eta = &
        1.3169578969248167086250463473079684_dp
    real(dp), parameter :: target_eta = 0.35_dp
    real(dp), parameter :: wave_number = 0.2_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/toroidal_analytical"
    real(dp) :: theta_trace(trace_points), potential(trace_points)
    real(dp) :: axial_z(trace_points)
    real(dp) :: bem_potential(trace_points)
    real(dp) :: field_norm(trace_points), bem_field_norm(trace_points)
    real(dp) :: theta_edges(theta_cells + 1), phi_edges(phi_cells + 1)
    real(dp) :: surface_value(phi_cells, theta_cells)
    real(dp) :: bem_error(phi_cells, theta_cells)
    real(dp) :: x(surface_points), y(surface_points), z(surface_points)
    real(dp), allocatable :: boundary_trace(:), boundary_flux(:)
    real(dp), allocatable :: parameters(:, :)
    complex(dp), allocatable :: helmholtz_trace(:), helmholtz_flux(:)
    complex(dp) :: helmholtz_exact(trace_points)
    complex(dp) :: helmholtz_bem(trace_points)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: field(3), theta_value, phi_value, denominator, target(3)
    real(dp) :: bem_gradient(3), start_time, end_time
    real(dp) :: bem_seconds, exact_seconds, helmholtz_seconds
    real(dp) :: poisson_solve_seconds, helmholtz_solve_seconds
    integer :: i, j, point, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call generate_torus_surface_mesh( &
        2.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), 13, 15, &
        vertices, triangles, parameters)
    call build_boundary_data()
    call cpu_time(start_time)
    call solve_laplace_bem_dtn_torus_curved_3d( &
        parameters, triangles, 2.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), &
        boundary_trace, 7, boundary_flux, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "Toroidal Laplace DtN solve failed"
    poisson_solve_seconds = end_time - start_time
    call generate_trace()
    call generate_helmholtz_trace()
    call generate_surface_map()
    call generate_geometry()
    call write_trace_data()

contains

    subroutine generate_trace()
        call cpu_time(start_time)
        do i = 1, trace_points
            theta_trace(i) = 2.0_dp*pi*real(i - 1, dp)/ &
                real(trace_points - 1, dp)
            call evaluate_toroidal_harmonic_p( &
                degree_index, order, target_eta, theta_trace(i), 0.4_dp, &
                potential(i), status)
            if (status /= 0) error stop "Toroidal potential evaluation failed"
            call evaluate_toroidal_ampere_field_p( &
                degree_index, order, scale, target_eta, theta_trace(i), 0.4_dp, &
                field, status)
            if (status /= 0) error stop "Toroidal Ampere evaluation failed"
            field_norm(i) = norm2(field)
        end do
        call cpu_time(end_time)
        exact_seconds = end_time - start_time

        call cpu_time(start_time)
        do i = 1, trace_points
            call toroidal_point_to_cartesian( &
                scale, target_eta, theta_trace(i), 0.4_dp, target)
            call evaluate_laplace_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, target, &
                7, bem_potential(i), status)
            if (status /= 0) error stop "Toroidal BEM evaluation failed"
            call numerical_laplace_gradient(target, bem_gradient)
            bem_field_norm(i) = norm2(bem_gradient)
        end do
        call cpu_time(end_time)
        bem_seconds = end_time - start_time

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(theta_trace, potential, label="Poisson analytical", &
            linestyle="-")
        call plot(theta_trace, bem_potential, label="Poisson BEM", &
            linestyle="--")
        call plot(theta_trace, field_norm, label="Ampere analytical |H|", &
            linestyle="-.")
        call plot(theta_trace, bem_field_norm, label="Ampere BEM |H|", &
            linestyle=":")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("normalized mode value")
        call title("Toroidal Poisson and Ampere: analytical vs BEM")
        call legend()
        call savefig(output_directory//"/toroidal_trace_1d.png")
    end subroutine generate_trace

    subroutine generate_surface_map()
        do i = 1, theta_cells + 1
            theta_edges(i) = 2.0_dp*pi*real(i - 1, dp)/ &
                real(theta_cells, dp)
        end do
        do j = 1, phi_cells + 1
            phi_edges(j) = 2.0_dp*pi*real(j - 1, dp)/ &
                real(phi_cells, dp)
        end do
        do i = 1, theta_cells
            theta_value = 0.5_dp*(theta_edges(i) + theta_edges(i + 1))
            do j = 1, phi_cells
                phi_value = 0.5_dp*(phi_edges(j) + phi_edges(j + 1))
                call evaluate_toroidal_harmonic_p( &
                    degree_index, order, target_eta, theta_value, phi_value, &
                    surface_value(j, i), status)
                if (status /= 0) error stop "Toroidal map evaluation failed"
                call toroidal_point_to_cartesian( &
                    scale, target_eta, theta_value, phi_value, target)
                call evaluate_laplace_representation_torus_curved_3d( &
                    parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                    1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, &
                    target, 7, bem_potential(1), status)
                if (status /= 0) error stop "Toroidal BEM map failed"
                bem_error(j, i) = &
                    abs(bem_potential(1) - surface_value(j, i))
            end do
        end do
        bem_error = bem_error/max(maxval(abs(surface_value)), tiny(1.0_dp))

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(theta_edges, phi_edges, surface_value, cmap="coolwarm")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle phi [rad]")
        call title("Poisson analytical toroidal surface mode")
        call savefig(output_directory//"/toroidal_surface_2d.png")

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(theta_edges, phi_edges, bem_error, cmap="viridis")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle phi [rad]")
        call title("Poisson BEM normalized absolute error")
        call savefig(output_directory//"/toroidal_bem_error_2d.png")
    end subroutine generate_surface_map

    subroutine generate_geometry()
        point = 0
        do j = 1, phi_cells
            phi_value = 2.0_dp*pi*real(j - 1, dp)/real(phi_cells, dp)
            do i = 1, theta_cells
                theta_value = 2.0_dp*pi*real(i - 1, dp)/ &
                    real(theta_cells, dp)
                denominator = cosh(eta) - cos(theta_value)
                point = point + 1
                x(point) = scale*sinh(eta)*cos(phi_value)/denominator
                y(point) = scale*sinh(eta)*sin(phi_value)/denominator
                z(point) = scale*sin(theta_value)/denominator
            end do
        end do

        call figure(figsize=[7.5_dp, 6.5_dp])
        call add_scatter(x, y, z, label="eta=acosh(2) torus", marker=".")
        call title("Toroidal exterior interface")
        call savefig(output_directory//"/toroidal_geometry_3d.png")
    end subroutine generate_geometry

    subroutine generate_helmholtz_trace()
        real(dp) :: radius, source(3)

        source = [2.0_dp/sqrt(3.0_dp), 0.0_dp, 0.0_dp]
        allocate(helmholtz_trace(size(vertices, 2)))
        call build_helmholtz_boundary_data(source)
        call cpu_time(start_time)
        call solve_helmholtz_bem_dtn_torus_curved_3d( &
            parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
            1.0_dp/sqrt(3.0_dp), wave_number, helmholtz_trace, 7, &
            helmholtz_flux, status)
        call cpu_time(end_time)
        if (status /= 0) error stop "Toroidal Helmholtz DtN solve failed"
        helmholtz_solve_seconds = end_time - start_time
        call cpu_time(start_time)
        do i = 1, trace_points
            axial_z(i) = &
                -1.5_dp + 3.0_dp*real(i - 1, dp)/real(trace_points - 1, dp)
            target = [0.0_dp, 0.0_dp, axial_z(i)]
            radius = norm2(target - source)
            helmholtz_exact(i) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
            call evaluate_helmholtz_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), helmholtz_trace, helmholtz_flux, &
                target, wave_number, 7, helmholtz_bem(i), status)
            if (status /= 0) error stop "Toroidal Helmholtz BEM failed"
        end do
        call cpu_time(end_time)
        helmholtz_seconds = end_time - start_time

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(axial_z, abs(helmholtz_exact), &
            label="Helmholtz analytical |u|", linestyle="-")
        call plot(axial_z, abs(helmholtz_bem), &
            label="Helmholtz BEM |u|", linestyle="--")
        call xlabel("axial coordinate z")
        call ylabel("field magnitude")
        call title("Outgoing Helmholtz point source through toroidal BEM")
        call legend()
        call savefig(output_directory//"/toroidal_helmholtz_1d.png")
    end subroutine generate_helmholtz_trace

    subroutine write_trace_data()
        open (newunit=unit, &
            file=output_directory//"/toroidal_trace.csv", &
            status="replace", action="write")
        write (unit, "(a)") &
            "theta,potential_exact,potential_bem,H_norm_exact,H_norm_bem"
        do i = 1, trace_points
            write (unit, "(*(es24.16,:,','))") &
                theta_trace(i), potential(i), bem_potential(i), &
                field_norm(i), bem_field_norm(i)
        end do
        close (unit)
        open (newunit=unit, file=output_directory//"/benchmark.txt", &
            status="replace", action="write")
        write (unit, "(a,i0)") "boundary triangles: ", size(triangles, 2)
        write (unit, "(a,i0)") "trace targets: ", trace_points
        write (unit, "(a,es14.6)") "analytical seconds: ", exact_seconds
        write (unit, "(a,es14.6)") &
            "Poisson BEM DtN solve seconds: ", poisson_solve_seconds
        write (unit, "(a,es14.6)") "BEM seconds: ", bem_seconds
        write (unit, "(a,es14.6)") &
            "Helmholtz BEM DtN solve seconds: ", helmholtz_solve_seconds
        write (unit, "(a,es14.6)") &
            "Helmholtz BEM seconds: ", helmholtz_seconds
        write (unit, "(a,es14.6)") "Poisson max relative error: ", maxval( &
            abs(bem_potential - potential))/maxval(abs(potential))
        write (unit, "(a,es14.6)") "Ampere max relative error: ", maxval( &
            abs(bem_field_norm - field_norm))/maxval(field_norm)
        write (unit, "(a,es14.6)") "Helmholtz max relative error: ", maxval( &
            abs(helmholtz_bem - helmholtz_exact))/ &
            maxval(abs(helmholtz_exact))
        close (unit)
        if (maxval(abs(bem_potential - potential))/ &
            maxval(abs(potential)) >= 3.0e-1_dp) &
            error stop "Poisson BEM gallery accuracy regression"
        if (maxval(abs(bem_field_norm - field_norm))/ &
            maxval(field_norm) >= 3.0e-1_dp) &
            error stop "Ampere BEM gallery accuracy regression"
        if (maxval(abs(helmholtz_bem - helmholtz_exact))/ &
            maxval(abs(helmholtz_exact)) >= 3.0e-2_dp) &
            error stop "Helmholtz BEM gallery accuracy regression"
    end subroutine write_trace_data

    subroutine build_boundary_data()
        real(dp) :: eta_at_point

        allocate(boundary_trace(size(vertices, 2)))
        do point = 1, size(vertices, 2)
            call cartesian_to_toroidal( &
                vertices(:, point), scale, eta_at_point, theta_value, phi_value)
            call evaluate_toroidal_harmonic_p( &
                degree_index, order, eta, theta_value, phi_value, &
                boundary_trace(point), status)
            if (status /= 0) error stop "Toroidal boundary trace failed"
        end do
    end subroutine build_boundary_data

    subroutine build_helmholtz_boundary_data(source)
        real(dp), intent(in) :: source(3)
        real(dp) :: displacement(3), radius

        do point = 1, size(vertices, 2)
            displacement = vertices(:, point) - source
            radius = norm2(displacement)
            helmholtz_trace(point) = &
                exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
        end do
    end subroutine build_helmholtz_boundary_data

    subroutine numerical_laplace_gradient(evaluation_point, gradient)
        real(dp), intent(in) :: evaluation_point(3)
        real(dp), intent(out) :: gradient(3)

        real(dp), parameter :: step = 1.0e-4_dp
        real(dp) :: backward_point(3), backward_value
        real(dp) :: forward_point(3), forward_value
        integer :: coordinate

        do coordinate = 1, 3
            forward_point = evaluation_point
            backward_point = evaluation_point
            forward_point(coordinate) = forward_point(coordinate) + step
            backward_point(coordinate) = backward_point(coordinate) - step
            call evaluate_laplace_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, &
                forward_point, 7, forward_value, status)
            if (status /= 0) error stop "Forward BEM gradient failed"
            call evaluate_laplace_representation_torus_curved_3d( &
                parameters, triangles, 2.0_dp/sqrt(3.0_dp), &
                1.0_dp/sqrt(3.0_dp), boundary_trace, boundary_flux, &
                backward_point, 7, backward_value, status)
            if (status /= 0) error stop "Backward BEM gradient failed"
            gradient(coordinate) = -(forward_value - backward_value)/(2*step)
        end do
    end subroutine numerical_laplace_gradient

end program toroidal_analytical
