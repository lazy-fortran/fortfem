program regularized_surface_current_layer
    !! Resolved slab view of a tangential sheet-current regularization.
    use fortfem_api, only: evaluate_regularized_surface_current_integral, &
        evaluate_regularized_surface_current_layer
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, set_xscale, &
        set_yscale, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 3
    integer, parameter :: profile_point_count = 241
    integer, parameter :: convergence_count = 5
    integer, parameter :: maximum_point_count = 129
    integer, parameter :: benchmark_repetitions = 2000
    integer, parameter :: point_counts(convergence_count) = &
        [9, 17, 33, 65, 129]
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: thickness = 0.015_dp
    real(dp), parameter :: half_width = 5.0_dp*thickness
    real(dp), parameter :: sheet_current_vector(component_count) = &
        [0.0_dp, 1200.0_dp, -450.0_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/regularized_surface_current_layer"

    real(dp) :: signed_distance(maximum_point_count)
    real(dp) :: normal_weights(maximum_point_count)
    real(dp) :: sheet_current(maximum_point_count, component_count)
    real(dp) :: volume_current(maximum_point_count, component_count)
    real(dp) :: profile_distance(profile_point_count)
    real(dp) :: profile_sheet_current(profile_point_count, component_count)
    real(dp) :: profile_volume_current(profile_point_count, component_count)
    real(dp) :: analytical_current(profile_point_count, component_count)
    real(dp) :: point_count_axis(convergence_count)
    real(dp) :: integral_error(convergence_count)
    real(dp) :: integrated_current(component_count), normalization
    real(dp) :: spacing, gaussian, current_norm, maximum_profile_error
    real(dp) :: start_time, end_time, action_seconds
    integer :: level, point, unit, repetition, point_count
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    current_norm = sqrt(sum(sheet_current_vector**2))

    do point = 1, profile_point_count
        profile_distance(point) = -half_width + 2.0_dp*half_width* &
            real(point - 1, dp)/real(profile_point_count - 1, dp)
        profile_sheet_current(point, :) = sheet_current_vector
        gaussian = exp(-(profile_distance(point)/thickness)**2)/ &
            (sqrt(pi)*thickness)
        analytical_current(point, :) = gaussian*sheet_current_vector
    end do
    call evaluate_regularized_surface_current_layer( &
        profile_distance, profile_sheet_current, thickness, &
        profile_volume_current, status)
    if (status%code /= 0) error stop "regularized sheet profile failed"
    maximum_profile_error = maxval(abs( &
        profile_volume_current - analytical_current))/ &
        maxval(abs(analytical_current))
    if (maximum_profile_error > 2.0e-14_dp) &
        error stop "regularized sheet profile disagrees with Gaussian oracle"

    do level = 1, convergence_count
        point_count = point_counts(level)
        point_count_axis(level) = real(point_count, dp)
        spacing = 2.0_dp*half_width/real(point_count - 1, dp)
        do point = 1, point_count
            signed_distance(point) = -half_width + real(point - 1, dp)*spacing
            normal_weights(point) = spacing
            sheet_current(point, :) = sheet_current_vector
        end do
        normal_weights(1) = 0.5_dp*spacing
        normal_weights(point_count) = 0.5_dp*spacing
        call evaluate_regularized_surface_current_integral( &
            signed_distance(:point_count), normal_weights(:point_count), &
            sheet_current(:point_count, :), thickness, normalization, &
            integrated_current, status)
        if (status%code /= 0) error stop "regularized sheet integral failed"
        integral_error(level) = sqrt(sum( &
            (integrated_current - sheet_current_vector)**2))/current_norm
    end do
    do level = 2, convergence_count
        if (integral_error(level) >= integral_error(level - 1)) &
            error stop "regularized sheet integral did not converge"
    end do
    if (integral_error(convergence_count) > 5.0e-12_dp) &
        error stop "regularized sheet integral accuracy regression"

    call cpu_time(start_time)
    do repetition = 1, benchmark_repetitions
        call evaluate_regularized_surface_current_layer( &
            profile_distance, profile_sheet_current, thickness, &
            profile_volume_current, status)
        if (status%code /= 0) error stop "regularized sheet benchmark failed"
    end do
    call cpu_time(end_time)
    action_seconds = (end_time - start_time)/real(benchmark_repetitions, dp)

    call render_plots()
    call write_profile_data()
    call write_convergence_data()
    call write_benchmark_metadata()

contains

    subroutine render_plots()
        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(profile_distance, profile_volume_current(:, 2), &
            label="FortFEM J_y", color=blue, linestyle="-", linewidth=2.0_dp)
        call plot(profile_distance(1:profile_point_count:12), &
            analytical_current(1:profile_point_count:12, 2), &
            label="Gaussian oracle J_y", color=blue, linestyle="None", &
            marker="o")
        call plot(profile_distance, profile_volume_current(:, 3), &
            label="FortFEM J_z", color=orange, linestyle="-", linewidth=2.0_dp)
        call plot(profile_distance(1:profile_point_count:12), &
            analytical_current(1:profile_point_count:12, 3), &
            label="Gaussian oracle J_z", color=orange, linestyle="None", &
            marker="s")
        call xlabel("signed normal distance d [m]")
        call ylabel("volume current density K delta_epsilon [A/m^2]")
        call title("Resolved tangential surface-current layer")
        call legend()
        call savefig(output_directory//"/regularized_sheet_profile_1d.png")

        call figure(figsize=[9.0_dp, 5.5_dp])
        call plot(point_count_axis, integral_error, &
            label="relative integrated-current error", color=blue, marker="o")
        call set_xscale("log")
        call set_yscale("log")
        call xlabel("normal-grid points")
        call ylabel("relative error in integral of J_epsilon")
        call title("Regularized sheet recovers the prescribed surface current")
        call legend()
        call savefig( &
            output_directory//"/regularized_sheet_integral_convergence.png")
    end subroutine render_plots

    subroutine write_profile_data()
        open (newunit=unit, file=output_directory// &
            "/regularized_sheet_profile.csv", status="replace", &
            action="write")
        write (unit, "(a)") &
            "distance_m,Jx_A_per_m2,Jy_A_per_m2,Jz_A_per_m2,"// &
            "oracle_Jx_A_per_m2,oracle_Jy_A_per_m2,oracle_Jz_A_per_m2"
        do point = 1, profile_point_count
            write (unit, "(*(es24.16,:,','))") profile_distance(point), &
                profile_volume_current(point, :), analytical_current(point, :)
        end do
        close (unit)
    end subroutine write_profile_data

    subroutine write_convergence_data()
        open (newunit=unit, file=output_directory// &
            "/regularized_sheet_integral_convergence.csv", status="replace", &
            action="write")
        write (unit, "(a)") "normal_grid_points,relative_integrated_current_error"
        do level = 1, convergence_count
            write (unit, "(i0,',',es24.16)") &
                point_counts(level), integral_error(level)
        end do
        close (unit)
    end subroutine write_convergence_data

    subroutine write_benchmark_metadata()
        open (newunit=unit, file=output_directory//"/benchmark.json", &
            status="replace", action="write")
        write (unit, "(a)") "{"
        write (unit, "(a)") &
            '  "schema": "fortfem-example-benchmark-v1",'
        write (unit, "(a)") &
            '  "kernel": "evaluate_regularized_surface_current_layer",'
        write (unit, "(a)") &
            '  "oracle": "exp(-(d/epsilon)^2)/(sqrt(pi)*epsilon)",'
        write (unit, "(a,es24.16,a)") &
            '  "thickness_m": ', thickness, ','
        write (unit, "(a,i0,a)") &
            '  "profile_points": ', profile_point_count, ','
        write (unit, "(a,i0,a)") &
            '  "repetitions": ', benchmark_repetitions, ','
        write (unit, "(a,es24.16,a)") &
            '  "seconds_per_profile": ', action_seconds, ','
        write (unit, "(a,es24.16,a)") &
            '  "maximum_relative_profile_error": ', &
            maximum_profile_error, ','
        write (unit, "(a,es24.16)") &
            '  "finest_relative_integrated_current_error": ', &
            integral_error(convergence_count)
        write (unit, "(a)") "}"
        close (unit)
    end subroutine write_benchmark_metadata

end program regularized_surface_current_layer
