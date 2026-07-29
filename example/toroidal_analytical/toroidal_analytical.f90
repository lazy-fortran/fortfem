program toroidal_analytical
    use fortfem_api, only: &
        evaluate_toroidal_harmonic_p, evaluate_toroidal_ampere_field_p
    use fortfem_kinds, only: dp
    use fortplot, only: figure, plot, pcolormesh, add_scatter, &
        xlabel, ylabel, title, legend, savefig
    implicit none

    integer, parameter :: trace_points = 181
    integer, parameter :: theta_cells = 72
    integer, parameter :: phi_cells = 48
    integer, parameter :: surface_points = theta_cells*phi_cells
    integer, parameter :: degree_index = 2
    integer, parameter :: order = 1
    real(dp), parameter :: scale = 1.0_dp
    real(dp), parameter :: eta = &
        1.3169578969248167086250463473079684_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/toroidal_analytical"
    real(dp) :: theta_trace(trace_points), potential(trace_points)
    real(dp) :: field_eta(trace_points), field_theta(trace_points)
    real(dp) :: field_phi(trace_points)
    real(dp) :: theta_edges(theta_cells + 1), phi_edges(phi_cells + 1)
    real(dp) :: surface_value(phi_cells, theta_cells)
    real(dp) :: x(surface_points), y(surface_points), z(surface_points)
    real(dp) :: field(3), theta_value, phi_value, denominator
    integer :: i, j, point, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    call generate_trace()
    call generate_surface_map()
    call generate_geometry()
    call write_trace_data()

contains

    subroutine generate_trace()
        do i = 1, trace_points
            theta_trace(i) = 2.0_dp*pi*real(i - 1, dp)/ &
                real(trace_points - 1, dp)
            call evaluate_toroidal_harmonic_p( &
                degree_index, order, eta, theta_trace(i), 0.4_dp, &
                potential(i), status)
            if (status /= 0) error stop "Toroidal potential evaluation failed"
            call evaluate_toroidal_ampere_field_p( &
                degree_index, order, scale, eta, theta_trace(i), 0.4_dp, &
                field, status)
            if (status /= 0) error stop "Toroidal Ampere evaluation failed"
            field_eta(i) = field(1)
            field_theta(i) = field(2)
            field_phi(i) = field(3)
        end do

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(theta_trace, potential, label="Poisson potential", &
            linestyle="-")
        call plot(theta_trace, field_eta, label="Ampere H_eta", &
            linestyle="--")
        call plot(theta_trace, field_theta, label="Ampere H_theta", &
            linestyle="-.")
        call plot(theta_trace, field_phi, label="Ampere H_phi", &
            linestyle=":")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("normalized mode value")
        call title("Toroidal n=2, m=1 analytical reference")
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
                    degree_index, order, eta, theta_value, phi_value, &
                    surface_value(j, i), status)
                if (status /= 0) error stop "Toroidal map evaluation failed"
            end do
        end do

        call figure(figsize=[8.5_dp, 5.5_dp])
        call pcolormesh(theta_edges, phi_edges, surface_value, cmap="coolwarm")
        call xlabel("poloidal angle theta [rad]")
        call ylabel("toroidal angle phi [rad]")
        call title("Poisson toroidal surface mode")
        call savefig(output_directory//"/toroidal_surface_2d.png")
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
        write (unit, "(a)") "theta,potential,H_eta,H_theta,H_phi"
        do i = 1, trace_points
            write (unit, "(*(es24.16,:,','))") &
                theta_trace(i), potential(i), field_eta(i), &
                field_theta(i), field_phi(i)
        end do
        close (unit)
    end subroutine write_trace_data

end program toroidal_analytical
