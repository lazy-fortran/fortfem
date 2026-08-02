program direct_force_campaign_gallery
    !! Physical-first neutral direct-force objective campaign on a torus.
    use fortfem_feec, only: &
        evaluate_force_balance_objective, &
        evaluate_force_balance_objective_jvp
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_parametric_surface, add_scatter, &
        colorbar, figure, legend, pcolormesh, plot, savefig, title, view_init, &
        xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: theta_count = 25, phi_count = 33
    integer, parameter :: sample_count = theta_count*phi_count
    integer, parameter :: arrow_theta_stride = 4, arrow_phi_stride = 4
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: major_radius = 2.4_dp, minor_radius = 0.7_dp
    real(dp), parameter :: fd_step = 2.0e-7_dp
    character(*), parameter :: output_directory = &
        "output/example/direct_force_campaign_gallery"

    real(dp) :: theta(theta_count), phi(phi_count)
    real(dp) :: surface_x(theta_count, phi_count)
    real(dp) :: surface_y(theta_count, phi_count)
    real(dp) :: surface_z(theta_count, phi_count)
    real(dp) :: state(theta_count, phi_count)
    real(dp) :: residual(sample_count), residual_dot(sample_count)
    real(dp) :: weights(sample_count), weights_dot(sample_count)
    real(dp) :: force_x(theta_count, phi_count), force_y(theta_count, phi_count)
    real(dp) :: force_z(theta_count, phi_count)
    real(dp) :: objective, objective_dot, objective_plus, objective_minus
    real(dp) :: fd_error, radius_error, start_time, elapsed_seconds
    real(dp) :: line_x(2), line_y(2), line_z(2)
    real(dp) :: force_scale
    integer :: i, j, point, unit, command_status
    type(fortsparse_status_t) :: status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, &
        exitstat=command_status)
    if (command_status /= 0) error stop "cannot create direct-force output directory"
    call initialize_gallery_sequence()

    do j = 1, phi_count
        phi(j) = 2.0_dp*pi*real(j - 1, dp)/real(phi_count - 1, dp)
    end do
    do i = 1, theta_count
        theta(i) = 2.0_dp*pi*real(i - 1, dp)/real(theta_count - 1, dp)
    end do

    point = 0
    radius_error = 0.0_dp
    do j = 1, phi_count
        do i = 1, theta_count
            point = point + 1
            surface_x(i, j) = (major_radius + minor_radius*cos(theta(i)))*cos(phi(j))
            surface_y(i, j) = (major_radius + minor_radius*cos(theta(i)))*sin(phi(j))
            surface_z(i, j) = minor_radius*sin(theta(i))
            state(i, j) = 0.5_dp + 0.20_dp*cos(theta(i) - 2.0_dp*phi(j)) + &
                0.10_dp*sin(2.0_dp*theta(i) + phi(j))
            force_x(i, j) = 0.08_dp*sin(theta(i))*cos(phi(j))
            force_y(i, j) = 0.06_dp*cos(theta(i))*sin(phi(j))
            force_z(i, j) = 0.05_dp*sin(theta(i) + phi(j))
            residual(point) = force_x(i, j) + 0.5_dp*force_y(i, j) - &
                0.25_dp*force_z(i, j)
            residual_dot(point) = 0.015_dp*cos(theta(i) - phi(j))
            weights(point) = 1.0_dp + 0.15_dp*cos(theta(i))**2
            weights_dot(point) = 0.02_dp*sin(phi(j))
            radius_error = max(radius_error, abs(sqrt(surface_x(i, j)**2 + &
                surface_y(i, j)**2) - major_radius - minor_radius*cos(theta(i))))
        end do
    end do

    call cpu_time(start_time)
    call evaluate_force_balance_objective( &
        residual, weights, objective, status)
    if (status%code /= 0) error stop "direct-force objective failed"
    call evaluate_force_balance_objective_jvp( &
        residual, weights, residual_dot, weights_dot, objective_dot, status)
    if (status%code /= 0) error stop "direct-force objective JVP failed"
    call evaluate_force_balance_objective( &
        residual + fd_step*residual_dot, weights + fd_step*weights_dot, &
        objective_plus, status)
    call evaluate_force_balance_objective( &
        residual - fd_step*residual_dot, weights - fd_step*weights_dot, &
        objective_minus, status)
    fd_error = abs(objective_dot - &
        (objective_plus - objective_minus)/(2.0_dp*fd_step))
    call cpu_time(elapsed_seconds)
    elapsed_seconds = elapsed_seconds - start_time
    if (fd_error > 2.0e-9_dp .or. radius_error > 2.0e-12_dp) &
        error stop "direct-force analytical oracle failed"

    ! Physical solution first: the torus is the state carrier and short force
    ! vectors make the supplied residual visible without implying a closure.
    call figure(figsize=[9.0_dp, 7.0_dp])
    call add_parametric_surface(surface_x, surface_y, surface_z, &
        color="lightsteelblue", alpha=0.70_dp, linewidth=0.0_dp, &
        filled=.true., label="manufactured toroidal surface")
    call add_scatter(reshape(surface_x, [sample_count]), &
        reshape(surface_y, [sample_count]), reshape(surface_z, [sample_count]), &
        c=reshape(state, [sample_count]), cmap="viridis", marker=".", &
        markersize=3.0_dp, label="manufactured state")
    force_scale = 0.38_dp/max(1.0e-12_dp, maxval(sqrt(force_x**2 + &
        force_y**2 + force_z**2)))
    do j = 1, phi_count - 1, arrow_phi_stride
        do i = 1, theta_count - 1, arrow_theta_stride
            line_x = [surface_x(i, j), surface_x(i, j) + force_scale*force_x(i, j)]
            line_y = [surface_y(i, j), surface_y(i, j) + force_scale*force_y(i, j)]
            line_z = [surface_z(i, j), surface_z(i, j) + force_scale*force_z(i, j)]
            call add_3d_plot(line_x, line_y, line_z, color="black", linewidth=1.2_dp)
        end do
    end do
    call colorbar(label="manufactured state u")
    call xlabel("x [m]")
    call ylabel("y [m]")
    call title("Direct-force foundation: toroidal state and supplied force vectors")
    call view_init(elev=27.0_dp, azim=-52.0_dp)
    call savefig(output_directory//"/direct_force_torus_solution_3d.png")

    call figure(figsize=[8.5_dp, 5.8_dp])
    call pcolormesh(theta, phi, transpose(state), cmap="viridis")
    call colorbar(label="manufactured state u")
    call xlabel("poloidal angle theta [rad]")
    call ylabel("toroidal angle phi [rad]")
    call title("Same toroidal state in parameter coordinates")
    call savefig(output_directory//"/direct_force_state_2d.png")

    call figure(figsize=[8.5_dp, 5.8_dp])
    call plot([0.0_dp, 1.0_dp], [objective, objective], &
        label="objective", linewidth=2.2_dp)
    call plot([0.0_dp, 1.0_dp], [objective_dot, objective_dot], &
        label="JVP directional derivative", linestyle="--", linewidth=2.0_dp)
    call plot([0.0_dp, 1.0_dp], [fd_error, fd_error], &
        label="centered-difference error", linestyle=":", linewidth=2.0_dp)
    call xlabel("diagnostic coordinate")
    call ylabel("value")
    call title("Direct-force objective and derivative diagnostics")
    call legend()
    call savefig(output_directory//"/direct_force_objective_diagnostics_1d.png")

    call write_csv()
    open (newunit=unit, file=output_directory//"/benchmark.json", &
        status="replace", action="write")
    write (unit, "(a)") "{" // &
        '"samples":'//trim(adjustl(itoa(sample_count)))//","// &
        '"objective":'//trim(adjustl(rtoa(objective)))//","// &
        '"jvp":'//trim(adjustl(rtoa(objective_dot)))//","// &
        '"fd_error":'//trim(adjustl(rtoa(fd_error)))//","// &
        '"geometry_error":'//trim(adjustl(rtoa(radius_error)))//","// &
        '"elapsed_seconds":'//trim(adjustl(rtoa(elapsed_seconds)))//","// &
        '"provenance":"analytic-torus-direct-force-v1",'// &
        '"primary_plot":"direct_force_torus_solution_3d.png",'// &
        '"closure":"neutral-caller-owned-residual"}'
    close (unit)

contains

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit
        open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        write (sequence_unit, "(a)") "physical_solution"
        write (sequence_unit, "(a)") "parameter_solution"
        write (sequence_unit, "(a)") "diagnostics"
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine write_csv()
        integer :: csv_unit, local_i, local_j, local_point
        open (newunit=csv_unit, file=output_directory//"/direct_force_surface.csv", &
            status="replace", action="write")
        write (csv_unit, "(a)") "theta,phi,x,y,z,state,force_x,force_y,force_z,residual,weight"
        local_point = 0
        do local_j = 1, phi_count
            do local_i = 1, theta_count
                local_point = local_point + 1
                write (csv_unit, "(11(es24.16e3,:,','))") theta(local_i), phi(local_j), &
                    surface_x(local_i, local_j), surface_y(local_i, local_j), &
                    surface_z(local_i, local_j), state(local_i, local_j), &
                    force_x(local_i, local_j), force_y(local_i, local_j), &
                    force_z(local_i, local_j), residual(local_point), weights(local_point)
            end do
        end do
        close (csv_unit)
        open (newunit=csv_unit, file=output_directory//"/direct_force_diagnostics.csv", &
            status="replace", action="write")
        write (csv_unit, "(a)") "objective,jvp,fd_error,geometry_error,elapsed_seconds"
        write (csv_unit, "(5(es24.16e3,:,','))") objective, objective_dot, fd_error, &
            radius_error, elapsed_seconds
        close (csv_unit)
    end subroutine write_csv

    function itoa(value) result(text)
        integer, intent(in) :: value
        character(32) :: text
        write (text, "(i0)") value
    end function itoa

    function rtoa(value) result(text)
        real(dp), intent(in) :: value
        character(64) :: text
        write (text, "(es24.16e3)") value
    end function rtoa

end program direct_force_campaign_gallery
