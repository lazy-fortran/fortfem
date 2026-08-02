program sheet_current_parity
    !! Physical-first slab fixture for explicit versus resolved sheet ledgers.
    use fortfem_api, only: evaluate_regularized_surface_current_layer, &
        evaluate_sheet_current_parity
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: component_count = 3, profile_count = 241
    integer, parameter :: parity_count = 401
    real(dp), parameter :: pi = acos(-1.0_dp), thickness = 0.02_dp
    real(dp), parameter :: half_width = 5.0_dp*thickness
    real(dp), parameter :: surface_measure = 2.5_dp
    real(dp), parameter :: sheet(component_count) = [0.0_dp, 1.25_dp, -0.5_dp]
    real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
    real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
    character(*), parameter :: output_directory = &
        "output/example/sheet_current_parity"

    real(dp) :: profile_distance(profile_count), profile_current(profile_count, 3)
    real(dp) :: profile_volume(profile_count, 3), profile_oracle(profile_count, 3)
    real(dp) :: distance(parity_count), weights(parity_count), current(parity_count, 3)
    real(dp) :: regularized(3), explicit(3), relative_error, spacing, gaussian
    real(dp) :: start_time, end_time, seconds_per_profile
    integer :: point, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    do point = 1, profile_count
        profile_distance(point) = -half_width + 2.0_dp*half_width* &
            real(point - 1, dp)/real(profile_count - 1, dp)
        profile_current(point, :) = sheet
        gaussian = exp(-(profile_distance(point)/thickness)**2)/ &
            (sqrt(pi)*thickness)
        profile_oracle(point, :) = gaussian*sheet
    end do
    call cpu_time(start_time)
    call evaluate_regularized_surface_current_layer( &
        profile_distance, profile_current, thickness, profile_volume, status)
    call cpu_time(end_time)
    if (status%code /= 0) error stop "sheet profile evaluation failed"
    seconds_per_profile = end_time - start_time

    spacing = 2.0_dp*half_width/real(parity_count - 1, dp)
    do point = 1, parity_count
        distance(point) = -half_width + spacing*real(point - 1, dp)
        weights(point) = spacing
        if (point == 1 .or. point == parity_count) weights(point) = 0.5_dp*spacing
        current(point, :) = sheet
    end do
    call evaluate_sheet_current_parity( &
        distance, weights, current, thickness, surface_measure, sheet, &
        regularized, explicit, relative_error, status)
    if (status%code /= 0 .or. relative_error > 2.0e-12_dp) &
        error stop "sheet representation parity failed"

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(profile_distance, profile_volume(:, 2), label="resolved J_y", &
        color=blue, linewidth=2.0_dp)
    call plot(profile_distance(1:profile_count:12), profile_oracle(1:profile_count:12, 2), &
        label="Gaussian oracle J_y", color=blue, linestyle="None", marker="o")
    call plot(profile_distance, profile_volume(:, 3), label="resolved J_z", &
        color=orange, linewidth=2.0_dp)
    call plot(profile_distance(1:profile_count:12), profile_oracle(1:profile_count:12, 3), &
        label="Gaussian oracle J_z", color=orange, linestyle="None", marker="s")
    call xlabel("signed normal distance d [m]")
    call ylabel("volume current density [A/m^2]")
    call title("Explicit sheet and resolved layer: physical slab profile")
    call legend()
    call savefig(output_directory//"/sheet_current_parity_profile_1d.png")

    open (newunit=unit, file=output_directory//"/sheet_current_parity_profile.csv", &
        status="replace", action="write")
    write (unit, "(a)") "distance_m,Jy_resolved,Jz_resolved,Jy_oracle,Jz_oracle"
    do point = 1, profile_count
        write (unit, "(*(es24.16,:,','))") profile_distance(point), &
            profile_volume(point, 2), profile_volume(point, 3), &
            profile_oracle(point, 2), profile_oracle(point, 3)
    end do
    close (unit)

    open (newunit=unit, file=output_directory//"/benchmark.json", &
        status="replace", action="write")
    write (unit, "(a)") "{"
    write (unit, "(a)") '  "schema": "fortfem-sheet-current-parity-v1",'
    write (unit, "(a,es24.16,a)") '  "relative_integrated_error": ', &
        relative_error, ','
    write (unit, "(a,es24.16)") '  "seconds_per_profile": ', seconds_per_profile
    write (unit, "(a)") "}"
    close (unit)
end program sheet_current_parity
