program maxwell_open_boundary_comparison
    use fortfem_api, only: &
        apply_planar_maxwell_dtn, &
        assemble_planar_nedelec_maxwell_dtn_form, &
        build_planar_nedelec_trace_sampling, &
        generate_structured_tetra_box_mesh, solve_tetra_nedelec_curl_mass
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: nx = 24, ny = 16, mode_count = 7
    integer, parameter :: trace_nx = 5, trace_ny = 5
    real(dp), parameter :: length_x = 2.0_dp*acos(-1.0_dp)
    real(dp), parameter :: length_y = 4.0_dp
    real(dp), parameter :: wave_number = 7.5_dp
    character(*), parameter :: output_directory = &
        "output/example/maxwell_open_boundary_comparison"

    complex(dp) :: derivative(2, nx, ny), trace(2, nx, ny)
    complex(dp), allocatable :: boundary_form(:, :), sampling(:, :)
    complex(dp), allocatable :: sampled_trace(:), weak_exact(:), weak_result(:)
    integer, allocatable :: boundary_dofs(:), tetrahedra(:, :)
    real(dp), allocatable :: coefficients(:), vertices(:, :)
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: exact_te(mode_count), exact_tm(mode_count)
    real(dp) :: numerical_te(mode_count), numerical_tm(mode_count)
    real(dp) :: coupling_error(4), orders(4), trace_error(4)
    real(dp) :: modes(mode_count), reflection(3), reflection_log(3)
    real(dp) :: method_index(3)
    real(dp) :: beta, end_time, phase_angle, seconds, start_time
    real(dp) :: max_modal_error
    real(dp) :: bounds(3, 2), origin(3), periods(3, 2), trace_weight
    integer :: boundary_count(4), i, j, mode, order, status, unit

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

    bounds(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    bounds(:, 2) = [length_x, length_y, 1.0_dp]
    origin = bounds(:, 1)
    periods(:, 1) = [length_x, 0.0_dp, 0.0_dp]
    periods(:, 2) = [0.0_dp, length_y, 0.0_dp]
    trace_weight = length_x*length_y/real(trace_nx*trace_ny, dp)
    call generate_structured_tetra_box_mesh( &
        bounds, [1, 1, 1], vertices, tetrahedra, status)
    if (status /= 0) error stop "Maxwell trace mesh failed"
    do order = 1, 4
        orders(order) = real(order, dp)
        call solve_tetra_nedelec_curl_mass( &
            vertices, tetrahedra, order, constant_source, 1.0_dp, 1.0_dp, &
            coefficients, sparse_status)
        if (sparse_status%code /= 0) error stop "Maxwell trace solve failed"
        call build_planar_nedelec_trace_sampling( &
            vertices, tetrahedra, order, origin, periods, trace_nx, trace_ny, &
            boundary_dofs, sampling, status)
        if (status /= 0) error stop "Maxwell trace sampling failed"
        sampled_trace = matmul( &
            sampling, cmplx(coefficients(boundary_dofs), 0.0_dp, dp))
        trace_error(order) = max( &
            maxval(abs(sampled_trace(1::2) - 1.25_dp)), &
            maxval(abs(sampled_trace(2::2) + 0.75_dp)))
        call assemble_planar_nedelec_maxwell_dtn_form( &
            vertices, tetrahedra, order, origin, periods, trace_nx, trace_ny, &
            wave_number, boundary_dofs, boundary_form, status)
        if (status /= 0) error stop "Maxwell trace DtN form failed"
        allocate(weak_exact(size(boundary_dofs)))
        weak_exact = matmul( &
            transpose(sampling), &
            -cmplx(0.0_dp, wave_number*trace_weight, dp)*sampled_trace)
        weak_result = matmul( &
            boundary_form, &
            cmplx(coefficients(boundary_dofs), 0.0_dp, dp))
        coupling_error(order) = maxval(abs(weak_result - weak_exact))
        boundary_count(order) = size(boundary_dofs)
        if (trace_error(order) >= 3.0e-10_dp .or. &
            coupling_error(order) >= 3.0e-8_dp) &
            error stop "Maxwell Nedelec-DtN coupling regression"
        deallocate( &
            coefficients, boundary_dofs, sampling, sampled_trace, &
            boundary_form, weak_exact, weak_result)
    end do

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

    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(orders, log10(max(trace_error, tiny(1.0_dp))), &
        label="constant trace error", marker="o")
    call plot(orders, log10(max(coupling_error, tiny(1.0_dp))), &
        label="pulled-back DtN error", marker="s")
    call xlabel("tetrahedral Nedelec order")
    call ylabel("log10 absolute error")
    call title("Automatic H(curl)-to-DtN coupling")
    call legend()
    call savefig(output_directory//"/maxwell_nedelec_dtn_1d.png")

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
    do order = 1, 4
        write (unit, "(a,i0,a,i0,a,es14.6,a,es14.6)") &
            "Nedelec order ", order, " boundary dofs: ", &
            boundary_count(order), " trace error: ", trace_error(order), &
            " pulled-back DtN error: ", coupling_error(order)
    end do
    close (unit)

contains

    pure subroutine constant_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)

        associate(unused => x + y + z)
        end associate
        value = [1.25_dp, -0.75_dp, 0.5_dp]
    end subroutine constant_source

end program maxwell_open_boundary_comparison
