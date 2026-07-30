program curved_acoustic_ntd
    use fortfem_api, only: &
        apply_curved_acoustic_displacement_ntd_2d, &
        circular_helmholtz_dtn_eigenvalue
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, set_yscale, title, &
        xlabel, ylabel
    implicit none

    integer, parameter :: mesh_sizes(4) = [12, 24, 48, 96]
    character(*), parameter :: output_directory = &
        "output/example/curved_acoustic_ntd"
    real(dp) :: errors(size(mesh_sizes)), mesh_widths(size(mesh_sizes))
    real(dp), allocatable :: angle(:)
    complex(dp), allocatable :: exact(:), pressure(:)
    integer :: mesh

    call execute_command_line("mkdir -p "//output_directory)
    do mesh = 1, size(mesh_sizes)
        call solve_circle_mode( &
            mesh_sizes(mesh), angle, pressure, exact, errors(mesh))
        mesh_widths(mesh) = 2.0_dp*acos(-1.0_dp)/real(mesh_sizes(mesh), dp)
    end do

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(angle, real(exact, dp), label="analytical Hankel mode", &
        linestyle="-")
    call plot(angle, real(pressure, dp), label="polygonal BEM NtD", &
        linestyle="--")
    call xlabel("boundary angle")
    call ylabel("Re(p)")
    call title("Curved acoustic displacement-to-pressure map")
    call legend()
    call savefig(output_directory//"/curved_acoustic_pressure.png")

    call figure(figsize=[8.0_dp, 5.0_dp])
    call plot(mesh_widths, errors, label="relative boundary L2 error")
    call set_yscale("log")
    call xlabel("angular panel width")
    call ylabel("relative error")
    call title("Curved acoustic NtD convergence")
    call legend()
    call savefig(output_directory//"/curved_acoustic_convergence.png")

contains

    subroutine solve_circle_mode( &
            panel_count, angles, numerical, analytical, relative_error)
        integer, intent(in) :: panel_count
        real(dp), allocatable, intent(out) :: angles(:)
        complex(dp), allocatable, intent(out) :: numerical(:), analytical(:)
        real(dp), intent(out) :: relative_error

        real(dp), parameter :: angular_frequency = 2.3_dp
        real(dp), parameter :: density = 1.7_dp
        real(dp), parameter :: radius = 0.8_dp
        real(dp), parameter :: wavenumber = 1.4_dp
        complex(dp), allocatable :: displacement(:)
        integer, allocatable :: panel_nodes(:, :)
        real(dp), allocatable :: panel_end(:, :), panel_start(:, :)
        complex(dp) :: dtn_eigenvalue, pressure_amplitude
        real(dp) :: midpoint_angle, pi
        integer :: panel, status

        allocate(angles(panel_count), numerical(panel_count))
        allocate(analytical(panel_count), displacement(panel_count))
        allocate(panel_nodes(2, panel_count))
        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        pi = acos(-1.0_dp)
        do panel = 1, panel_count
            angles(panel) = 2.0_dp*pi*real(panel - 1, dp)/ &
                real(panel_count, dp)
            panel_start(1, panel) = radius*cos(angles(panel))
            panel_start(2, panel) = radius*sin(angles(panel))
            panel_end(1, panel) = radius*cos( &
                2.0_dp*pi*real(panel, dp)/real(panel_count, dp))
            panel_end(2, panel) = radius*sin( &
                2.0_dp*pi*real(panel, dp)/real(panel_count, dp))
            panel_nodes(1, panel) = panel
            panel_nodes(2, panel) = modulo(panel, panel_count) + 1
            midpoint_angle = angles(panel) + pi/real(panel_count, dp)
            displacement(panel) = exp(cmplx(0.0_dp, midpoint_angle, dp))
        end do

        call circular_helmholtz_dtn_eigenvalue( &
            1, wavenumber, radius, dtn_eigenvalue, status)
        if (status /= 0) error stop "analytical circular DtN failed"
        pressure_amplitude = density*angular_frequency**2/dtn_eigenvalue
        analytical = pressure_amplitude*exp(cmplx(0.0_dp, angles, dp))
        call apply_curved_acoustic_displacement_ntd_2d( &
            panel_start, panel_end, panel_nodes, wavenumber, density, &
            angular_frequency, 16, displacement, numerical, status)
        if (status /= 0) error stop "curved acoustic NtD solve failed"
        relative_error = sqrt( &
            sum(abs(numerical - analytical)**2)/sum(abs(analytical)**2))
    end subroutine solve_circle_mode

end program curved_acoustic_ntd
