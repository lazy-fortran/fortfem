program maxwell_open_boundary_comparison
    use fortfem_api, only: apply_planar_maxwell_dtn
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    implicit none

    integer, parameter :: nx = 24, ny = 16, mode_count = 7
    real(dp), parameter :: length_x = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: length_y = 4.0_dp
    real(dp), parameter :: wave_number = 7.5_dp
    character(*), parameter :: output_directory = &
        "output/example/maxwell_open_boundary_comparison"

    complex(dp) :: derivative(2, nx, ny), trace(2, nx, ny)
    real(dp) :: exact_te(mode_count), exact_tm(mode_count)
    real(dp) :: numerical_te(mode_count), numerical_tm(mode_count)
    real(dp) :: modes(mode_count), reflection(3), reflection_log(3)
    real(dp) :: method_index(3)
    real(dp) :: beta, end_time, phase_angle, seconds, start_time
    real(dp) :: max_modal_error
    integer :: i, j, mode, status, unit

    call execute_command_line("mkdir -p "//output_directory)
    do mode = 0, mode_count - 1
        beta = sqrt(wave_number**2 - real(mode, dp)**2)
        modes(mode + 1) = real(mode, dp)
        exact_te(mode + 1) = beta
        exact_tm(mode + 1) = wave_number**2/beta

        trace = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, ny
            do i = 1, nx
                phase_angle = 2.0_dp*acos(-1.0_dp)*real(mode*(i - 1), dp)/ &
                    real(nx, dp)
                trace(2, i, j) = exp(cmplx(0.0_dp, phase_angle, dp))
            end do
        end do
        call apply_planar_maxwell_dtn( &
            trace, wave_number, length_x, length_y, derivative, status)
        if (status /= 0) error stop "Maxwell TE DtN gallery failed"
        numerical_te(mode + 1) = abs(derivative(2, 1, 1))

        trace = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, ny
            do i = 1, nx
                phase_angle = 2.0_dp*acos(-1.0_dp)*real(mode*(i - 1), dp)/ &
                    real(nx, dp)
                trace(1, i, j) = exp(cmplx(0.0_dp, phase_angle, dp))
            end do
        end do
        call apply_planar_maxwell_dtn( &
            trace, wave_number, length_x, length_y, derivative, status)
        if (status /= 0) error stop "Maxwell TM DtN gallery failed"
        numerical_tm(mode + 1) = abs(derivative(1, 1, 1))
    end do
    max_modal_error = max( &
        maxval(abs(numerical_te - exact_te)), &
        maxval(abs(numerical_tm - exact_tm)))
    if (max_modal_error >= 1.0e-11_dp) &
        error stop "Maxwell DtN modal accuracy regression"

    trace = cmplx(0.0_dp, 0.0_dp, dp)
    do j = 1, ny
        do i = 1, nx
            trace(1, i, j) = exp(cmplx(0.0_dp, &
                2.0_dp*acos(-1.0_dp)*real(2*(i - 1), dp)/real(nx, dp), dp))
            trace(2, i, j) = exp(cmplx(0.0_dp, &
                2.0_dp*acos(-1.0_dp)*real(3*(i - 1), dp)/real(nx, dp), dp))
        end do
    end do
    call cpu_time(start_time)
    call apply_planar_maxwell_dtn( &
        trace, wave_number, length_x, length_y, derivative, status)
    call cpu_time(end_time)
    if (status /= 0) error stop "Maxwell DtN benchmark failed"
    seconds = end_time - start_time

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(modes, exact_te, label="TE analytical", linestyle="-")
    call plot(modes, numerical_te, label="TE DtN", linestyle="--")
    call plot(modes, exact_tm, label="TM analytical", linestyle="-.")
    call plot(modes, numerical_tm, label="TM DtN", linestyle=":")
    call xlabel("tangential Fourier mode")
    call ylabel("capacity eigenvalue magnitude")
    call title("Biperiodic Maxwell DtN: analytical TE/TM spectrum")
    call legend()
    call savefig(output_directory//"/maxwell_dtn_modes_1d.png")

    method_index = [1.0_dp, 2.0_dp, 3.0_dp]
    reflection = [ &
        max_modal_error, exp(-2.0_dp*12.0_dp/3.0_dp), 1.0_dp]
    reflection_log = log10(max(reflection, tiny(1.0_dp)))
    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(method_index, reflection_log, label="reflection magnitude", &
        marker="o")
    call xlabel("method: 1=DtN, 2=PML, 3=far wall")
    call ylabel("log10 predicted/measured reflection")
    call title("Maxwell open-boundary reflection scale")
    call legend()
    call savefig(output_directory//"/maxwell_reflection_1d.png")

    open (newunit=unit, file=output_directory//"/benchmark.txt", &
        status="replace", action="write")
    write (unit, "(a,i0,a,i0)") "boundary samples: ", nx, " x ", ny
    write (unit, "(a,i0)") "tested TE/TM modes: ", mode_count
    write (unit, "(a,es14.6)") "DtN application seconds: ", seconds
    write (unit, "(a,es14.6)") &
        "maximum TE/TM eigenvalue error: ", max_modal_error
    write (unit, "(a,es14.6)") &
        "quadratic PML predicted reflection: ", reflection(2)
    write (unit, "(a,es14.6)") &
        "ordinary far-wall reflection: ", reflection(3)
    close (unit)
end program maxwell_open_boundary_comparison
