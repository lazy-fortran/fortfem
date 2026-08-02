program circular_dtn_modes
    use fortfem_boundary, only: apply_circular_helmholtz_dtn
    use fortfem_kinds, only: dp
    use fortplot, only: colorbar, figure, legend, pcolormesh, plot, savefig, &
        title, xlabel, ylabel
    implicit none

    integer, parameter :: point_count = 32
    character(*), parameter :: output_directory = &
        "output/example/circular_dtn_modes"
    complex(dp) :: normal_derivative(point_count), trace(point_count)
    real(dp) :: angle, discarded_relative_norm, theta(point_count)
    integer :: command_status, point, status

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"

    ! Boundary data on the unit circle:
    ! u(theta) = exp(2 i theta) + 0.1 exp(9 i theta).
    do point = 1, point_count
        angle = 2.0_dp * acos(-1.0_dp) * &
            real(point - 1, dp) / real(point_count, dp)
        theta(point) = angle
        trace(point) = exp(cmplx(0.0_dp, 2.0_dp * angle, dp)) + &
            0.1_dp * exp(cmplx(0.0_dp, 9.0_dp * angle, dp))
    end do

    ! Retaining |mode| <= 4 removes the small ninth harmonic. The diagnostic
    ! is the relative discrete L2 norm of the discarded boundary trace.
    call apply_circular_helmholtz_dtn( &
        trace, 3.0_dp, 1.0_dp, 4, normal_derivative, &
        discarded_relative_norm, status)
    if (status /= 0) error stop "Circular DtN application failed"

    write(*, '(a, es12.4)') &
        "discarded relative trace norm: ", discarded_relative_norm
    write(*, '(a, 2es14.5)') &
        "normal derivative at theta=0: ", normal_derivative(1)

    ! Render the physical field before one-dimensional diagnostics.  This is
    ! also the gallery's primary view of the example.
    call render_field()

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(theta, real(trace), label="Re(trace)", linestyle="-")
    call plot(theta, aimag(trace), label="Im(trace)", linestyle="--")
    call xlabel("boundary angle theta [rad]")
    call ylabel("boundary trace")
    call title("Circular Helmholtz DtN input trace")
    call legend()
    call savefig(output_directory//"/circular_dtn_trace_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(theta, real(normal_derivative), label="Re(DtN trace)", &
        marker="o")
    call plot(theta, aimag(normal_derivative), label="Im(DtN trace)", &
        marker="s")
    call xlabel("boundary angle theta [rad]")
    call ylabel("normal derivative")
    call title("Circular Helmholtz DtN response")
    call savefig(output_directory//"/circular_dtn_response_1d.png")

contains

    subroutine render_field()
        integer, parameter :: field_nx = 64, field_ny = 64
        real(dp) :: x_edges(field_nx + 1), y_edges(field_ny + 1)
        real(dp) :: field(field_ny, field_nx)
        real(dp) :: x_value, y_value, radius, angle
        real(dp) :: circle_x(field_nx + 1), circle_y(field_nx + 1)
        integer :: i, j

        do i = 1, field_nx + 1
            x_edges(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/field_nx
            circle_x(i) = cos(2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/field_nx)
            circle_y(i) = sin(2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/field_nx)
        end do
        y_edges = x_edges
        do j = 1, field_ny
            y_value = 0.5_dp*(y_edges(j) + y_edges(j + 1))
            do i = 1, field_nx
                x_value = 0.5_dp*(x_edges(i) + x_edges(i + 1))
                radius = sqrt(x_value*x_value + y_value*y_value)
                if (radius <= 1.0_dp) then
                    angle = atan2(y_value, x_value)
                    field(j, i) = real(exp(cmplx(0.0_dp, 2.0_dp*angle, dp)) &
                        *radius**2 + 0.1_dp*exp( &
                        cmplx(0.0_dp, 9.0_dp*angle, dp))*radius**9, dp)
                else
                    field(j, i) = 0.0_dp
                end if
            end do
        end do

        call figure(figsize=[8.0_dp, 6.5_dp])
        call pcolormesh(x_edges, y_edges, field, cmap="coolwarm")
        call colorbar(label="Re(u) in unit disk")
        call plot(circle_x, circle_y, color=[0.0_dp, 0.0_dp, 0.0_dp], &
            linewidth=1.5_dp)
        call xlabel("x")
        call ylabel("y")
        call title("Circular Helmholtz DtN input field")
        call savefig(output_directory//"/circular_dtn_field_2d.png")
    end subroutine render_field
end program circular_dtn_modes
