program eulerian_island_gallery
    !! Closure-neutral Eulerian island/separatrix foundation gallery.
    !!
    !! A fixed Cartesian sample carries an analytic slab-island flux.  The
    !! caller owns the force and divergence residual samples passed to the
    !! Eulerian non-nested composition; no plasma closure is selected here.
    use fortfem_feec, only: &
        assemble_eulerian_nonnested_residual, &
        assemble_eulerian_nonnested_residual_jvp
    use fortfem_time, only: CONTINUATION_EVENT_SIGN_CROSSING
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, contourf, figure, legend, plot, quiver, &
        savefig, title, xlabel, ylabel
    use fortsparse, only: FORTSPARSE_OK, fortsparse_status_t
    implicit none

    integer, parameter :: nx = 61, ny = 81
    integer, parameter :: point_count = nx*ny
    integer, parameter :: force_count = 2*point_count
    integer, parameter :: total_count = force_count + point_count
    integer, parameter :: quiver_nx = 11, quiver_ny = 13
    integer, parameter :: quiver_count = quiver_nx*quiver_ny
    real(dp), parameter :: white(3) = [1.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: cyan(3) = [0.0_dp, 1.0_dp, 1.0_dp]
    real(dp), parameter :: island_parameter = 0.40_dp
    real(dp), parameter :: x_min = -1.2_dp, x_max = 1.2_dp
    real(dp), parameter :: y_min = -0.8_dp, y_max = 0.8_dp
    character(*), parameter :: output_directory = &
        "output/example/eulerian_island_gallery"

    real(dp) :: x_grid(nx), y_grid(ny)
    real(dp) :: psi_map(nx, ny), bx_map(nx, ny), by_map(nx, ny)
    real(dp) :: force_norm_map(nx, ny), divergence_map(nx, ny)
    real(dp) :: force_residual(force_count), divergence_residual(point_count)
    real(dp) :: force_dot(force_count), divergence_dot(point_count)
    real(dp) :: residual(total_count), residual_dot(total_count)
    real(dp) :: previous_margin(2), current_margin(2), minimum_margin
    integer, parameter :: separatrix_count = 101
    real(dp) :: separatrix_x(separatrix_count), separatrix_y(separatrix_count)
    real(dp) :: separatrix_x_lower(separatrix_count)
    real(dp) :: separatrix_y_lower(separatrix_count)
    real(dp) :: quiver_x(quiver_count), quiver_y(quiver_count)
    real(dp) :: quiver_u(quiver_count), quiver_v(quiver_count)
    real(dp) :: maximum_divergence, residual_l2, residual_dot_l2
    real(dp) :: force_x, force_y, current_density, x, y, psi_value
    real(dp) :: field_norm, event_tolerance
    integer :: ix, iy, point, arrow, unit, event_code, event_index
    integer :: command_status
    type(fortsparse_status_t) :: status

    command_status = 0
    call execute_command_line("mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create island gallery output directory"
    call initialize_gallery_sequence()

    do iy = 1, ny
        y = y_min + (y_max - y_min)*real(iy - 1, dp)/real(ny - 1, dp)
        y_grid(iy) = y
        do ix = 1, nx
            x = x_min + (x_max - x_min)*real(ix - 1, dp)/real(nx - 1, dp)
            x_grid(ix) = x
            point = point_index(ix, iy)
            psi_value = island_flux(x, y)
            psi_map(ix, iy) = psi_value
            bx_map(ix, iy) = 2.0_dp*y
            by_map(ix, iy) = -4.0_dp*x*(x**2 - island_parameter)

            ! B is a rotated flux gradient, so its divergence is analytic zero.
            divergence_map(ix, iy) = 0.0_dp
            divergence_residual(point) = divergence_map(ix, iy)

            ! The caller-owned force sample is a manufactured J x B vector.
            ! It is diagnostic data only; this gallery does not choose a closure.
            current_density = -12.0_dp*x**2 + 4.0_dp*island_parameter - 2.0_dp
            force_x = current_density*by_map(ix, iy)
            force_y = -current_density*bx_map(ix, iy)
            force_residual(point) = force_x
            force_residual(point_count + point) = force_y
            force_norm_map(ix, iy) = sqrt(force_x**2 + force_y**2)
        end do
    end do

    force_dot = 0.025_dp*force_residual
    divergence_dot = 0.0_dp
    previous_margin = [-0.08_dp, 0.12_dp]
    current_margin = [0.04_dp, 0.12_dp]
    event_tolerance = 0.01_dp

    call assemble_eulerian_nonnested_residual( &
        force_residual, divergence_residual, residual, status, &
        previous_margin=previous_margin, current_margin=current_margin, &
        event_tolerance=event_tolerance, event_code=event_code, &
        event_index=event_index, minimum_margin=minimum_margin)
    if (status%code /= FORTSPARSE_OK) error stop "island residual composition failed"
    call assemble_eulerian_nonnested_residual_jvp( &
        force_residual, divergence_residual, force_dot, divergence_dot, &
        residual_dot, status, previous_margin=previous_margin, &
        current_margin=current_margin, event_tolerance=event_tolerance)
    if (status%code /= FORTSPARSE_OK) error stop "island residual JVP failed"

    maximum_divergence = maxval(abs(divergence_map))
    residual_l2 = sqrt(dot_product(residual, residual))
    residual_dot_l2 = sqrt(dot_product(residual_dot, residual_dot))
    if (event_code /= CONTINUATION_EVENT_SIGN_CROSSING .or. event_index /= 1) &
        error stop "island separatrix event oracle failed"

    do ix = 1, separatrix_count
        x = -sqrt(2.0_dp*island_parameter) + &
            2.0_dp*sqrt(2.0_dp*island_parameter)* &
            real(ix - 1, dp)/real(separatrix_count - 1, dp)
        separatrix_x(ix) = x
        separatrix_y(ix) = sqrt(max(0.0_dp, island_parameter**2 - &
            (x**2 - island_parameter)**2))
        separatrix_x_lower(ix) = x
        separatrix_y_lower(ix) = -separatrix_y(ix)
    end do
    arrow = 0
    do iy = 1, quiver_ny
        do ix = 1, quiver_nx
            arrow = arrow + 1
            quiver_x(arrow) = x_min + (x_max - x_min)*real(ix - 1, dp)/ &
                real(quiver_nx - 1, dp)
            quiver_y(arrow) = y_min + (y_max - y_min)*real(iy - 1, dp)/ &
                real(quiver_ny - 1, dp)
            quiver_u(arrow) = 2.0_dp*quiver_y(arrow)
            quiver_v(arrow) = -4.0_dp*quiver_x(arrow)* &
                (quiver_x(arrow)**2 - island_parameter)
            field_norm = sqrt(quiver_u(arrow)**2 + quiver_v(arrow)**2)
            if (field_norm > epsilon(1.0_dp)) then
                quiver_u(arrow) = 0.12_dp*quiver_u(arrow)/field_norm
                quiver_v(arrow) = 0.12_dp*quiver_v(arrow)/field_norm
            end if
        end do
    end do

    ! Physical state first: flux contours, separatrix, and tangent field.
    call figure(figsize=[8.0_dp, 6.5_dp])
    call contourf(x_grid, y_grid, psi_map, cmap="viridis", show_colorbar=.true.)
    call colorbar(label="island flux psi")
    call plot(separatrix_x, separatrix_y, color=white, linewidth=2.4_dp, &
        label="separatrix psi = psi_X")
    call plot(separatrix_x_lower, separatrix_y_lower, color=white, linewidth=2.4_dp)
    call quiver(quiver_x, quiver_y, quiver_u, quiver_v, scale=1.0_dp, &
        scale_units="xy", angles="xy", color="black", width=0.0025_dp, &
        headwidth=3.0_dp)
    call xlabel("slab x")
    call ylabel("slab y")
    call title("Eulerian island flux and separatrix field")
    call legend()
    call savefig(output_directory//"/island_flux_solution_2d.png")
    call record_gallery_stage("physical_solution")

    ! Diagnostics second: caller-owned force residual and event margin.
    call figure(figsize=[8.0_dp, 6.5_dp])
    call contourf(x_grid, y_grid, force_norm_map, cmap="magma", &
        show_colorbar=.true.)
    call colorbar(label="|caller force residual|")
    call plot(separatrix_x, separatrix_y, color=cyan, linewidth=2.0_dp)
    call plot(separatrix_x_lower, separatrix_y_lower, color=cyan, linewidth=2.0_dp)
    call xlabel("slab x")
    call ylabel("slab y")
    call title("Eulerian island residual diagnostic")
    call savefig(output_directory//"/island_flux_diagnostics_2d.png")
    call record_gallery_stage("diagnostics")

    call write_field_csv()
    call write_provenance_json(event_code, event_index, minimum_margin, &
        maximum_divergence, residual_l2, residual_dot_l2)

contains

    integer function point_index(i, j)
        integer, intent(in) :: i, j

        point_index = i + (j - 1)*nx
    end function point_index

    pure real(dp) function island_flux(x_value, y_value)
        real(dp), intent(in) :: x_value, y_value

        island_flux = (x_value**2 - island_parameter)**2 + y_value**2
    end function island_flux

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit, local_status

        open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot initialize gallery sequence"
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: sequence_unit, local_status

        open (newunit=sequence_unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot open gallery sequence"
        write (sequence_unit, "(a)", iostat=local_status) stage
        close (sequence_unit)
        if (local_status /= 0) error stop "cannot record gallery sequence"
    end subroutine record_gallery_stage

    subroutine write_field_csv()
        integer :: csv_unit, i, j, local_point

        open (newunit=csv_unit, file=output_directory//"/island_flux.csv", &
            status="replace", action="write")
        write (csv_unit, "(a)") "x,y,psi,bx,by,divergence,force_norm"
        do j = 1, ny
            do i = 1, nx
                local_point = point_index(i, j)
                write (csv_unit, "(*(es24.16,:,','))") x_grid(i), y_grid(j), &
                    psi_map(i, j), bx_map(i, j), by_map(i, j), &
                    divergence_residual(local_point), force_norm_map(i, j)
            end do
        end do
        close (csv_unit)
    end subroutine write_field_csv

    subroutine write_provenance_json(local_event_code, local_event_index, &
            local_minimum_margin, local_maximum_divergence, local_residual_l2, &
            local_residual_dot_l2)
        integer, intent(in) :: local_event_code, local_event_index
        real(dp), intent(in) :: local_minimum_margin, local_maximum_divergence
        real(dp), intent(in) :: local_residual_l2, local_residual_dot_l2
        integer :: json_unit

        open (newunit=json_unit, file=output_directory//"/provenance.json", &
            status="replace", action="write")
        write (json_unit, "(a)") "{"
        write (json_unit, "(a)") &
            '  "schema": "fortfem-eulerian-island-gallery-v1",'
        write (json_unit, "(a)") &
            '  "closure": "neutral-caller-owned-force-divergence",'
        write (json_unit, "(a)") &
            '  "primary_plot": "island_flux_solution_2d.png",'
        write (json_unit, "(a,i0,a)") '  "grid_x": ', nx, ","
        write (json_unit, "(a,i0,a)") '  "grid_y": ', ny, ","
        write (json_unit, "(a,es24.16,a)") '  "island_amplitude": ', &
            island_parameter, ","
        write (json_unit, "(a,es24.16,a)") '  "separatrix_level": ', &
            island_parameter**2, ","
        write (json_unit, "(a,i0,a)") '  "event_code": ', local_event_code, ","
        write (json_unit, "(a,i0,a)") '  "event_index": ', local_event_index, ","
        write (json_unit, "(a,es24.16,a)") '  "minimum_margin": ', &
            local_minimum_margin, ","
        write (json_unit, "(a,es24.16,a)") '  "maximum_divergence": ', &
            local_maximum_divergence, ","
        write (json_unit, "(a,es24.16,a)") '  "residual_l2": ', &
            local_residual_l2, ","
        write (json_unit, "(a,es24.16)") '  "residual_jvp_l2": ', &
            local_residual_dot_l2
        write (json_unit, "(a)") "}"
        close (json_unit)
    end subroutine write_provenance_json

end program eulerian_island_gallery
