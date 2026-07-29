program toroidal_analytical
    use fortfem_api, only: &
        cartesian_to_toroidal, evaluate_laplace_representation_triangles_3d, &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p, &
        generate_torus_surface_mesh, toroidal_point_to_cartesian, &
        toroidal_poisson_exterior_dtn_p, toroidal_vector_to_cartesian
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
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/toroidal_analytical"
    real(dp) :: theta_trace(trace_points), potential(trace_points)
    real(dp) :: bem_potential(trace_points)
    real(dp) :: field_eta(trace_points), field_theta(trace_points)
    real(dp) :: field_phi(trace_points)
    real(dp) :: field_norm(trace_points), bem_field_norm(trace_points)
    real(dp) :: theta_edges(theta_cells + 1), phi_edges(phi_cells + 1)
    real(dp) :: surface_value(phi_cells, theta_cells)
    real(dp) :: bem_error(phi_cells, theta_cells)
    real(dp) :: x(surface_points), y(surface_points), z(surface_points)
    real(dp), allocatable :: boundary_trace(:), boundary_flux(:)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: field(3), theta_value, phi_value, denominator, target(3)
    real(dp) :: bem_gradient(3), cartesian_field(3), start_time, end_time
    real(dp) :: bem_seconds, exact_seconds
    integer :: i, j, point, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call generate_torus_surface_mesh( &
        2.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp), 72, 56, &
        vertices, triangles)
    call build_boundary_data()
    call generate_trace()
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
            field_eta(i) = field(1)
            field_theta(i) = field(2)
            field_phi(i) = field(3)
            field_norm(i) = norm2(field)
        end do
        call cpu_time(end_time)
        exact_seconds = end_time - start_time

        call cpu_time(start_time)
        do i = 1, trace_points
            call toroidal_point_to_cartesian( &
                scale, target_eta, theta_trace(i), 0.4_dp, target)
            call evaluate_laplace_representation_triangles_3d( &
                vertices, triangles, boundary_trace, boundary_flux, target, &
                6, bem_potential(i), status, bem_gradient)
            if (status /= 0) error stop "Toroidal BEM evaluation failed"
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
                call evaluate_laplace_representation_triangles_3d( &
                    vertices, triangles, boundary_trace, boundary_flux, &
                    target, 6, bem_potential(1), status)
                if (status /= 0) error stop "Toroidal BEM map failed"
                bem_error(j, i) = abs(bem_potential(1) - surface_value(j, i))/ &
                    max(abs(surface_value(j, i)), 1.0e-8_dp)
            end do
        end do

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
        call title("Poisson BEM relative error")
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
        write (unit, "(a,es14.6)") "BEM seconds: ", bem_seconds
        write (unit, "(a,es14.6)") "Poisson max relative error: ", maxval( &
            abs(bem_potential - potential))/maxval(abs(potential))
        write (unit, "(a,es14.6)") "Ampere max relative error: ", maxval( &
            abs(bem_field_norm - field_norm))/maxval(field_norm)
        close (unit)
    end subroutine write_trace_data

    subroutine build_boundary_data()
        real(dp) :: centroid(3), eta_at_point, normal_derivative
        real(dp) :: dtn_value, value
        integer :: element

        allocate(boundary_trace(size(vertices, 2)))
        allocate(boundary_flux(size(triangles, 2)))
        do point = 1, size(vertices, 2)
            call cartesian_to_toroidal( &
                vertices(:, point), scale, eta_at_point, theta_value, phi_value)
            call evaluate_toroidal_harmonic_p( &
                degree_index, order, eta, theta_value, phi_value, &
                boundary_trace(point), status)
            if (status /= 0) error stop "Toroidal boundary trace failed"
        end do
        do element = 1, size(triangles, 2)
            centroid = sum(vertices(:, triangles(:, element)), dim=2)/3.0_dp
            call cartesian_to_toroidal( &
                centroid, scale, eta_at_point, theta_value, phi_value)
            call toroidal_poisson_exterior_dtn_p( &
                degree_index, order, scale, eta, theta_value, phi_value, &
                value, normal_derivative, dtn_value, status)
            if (status /= 0) error stop "Toroidal boundary flux failed"
            boundary_flux(element) = normal_derivative
        end do
    end subroutine build_boundary_data

end program toroidal_analytical
