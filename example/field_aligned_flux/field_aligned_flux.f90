program field_aligned_flux
    ! Manufactured pointwise field-aligned constitutive profile.
    use fortfem_feec, only: evaluate_field_aligned_flux
    use fortfem_kinds, only: dp
    use fortplot, only: figure, legend, plot, savefig, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: point_count = 128
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(*), parameter :: output_directory = &
        "output/example/field_aligned_flux"
    real(dp), parameter :: parallel_coefficient = 100.0_dp
    real(dp), parameter :: perpendicular_coefficient = 1.0_dp
    real(dp) :: coordinate(point_count), parallel_component(point_count)
    real(dp) :: perpendicular_component(point_count), flux_norm(point_count)
    real(dp) :: direction(3), gradient(3), flux(3)
    real(dp) :: angle
    integer :: point, unit
    type(fortsparse_status_t) :: status

    call execute_command_line("mkdir -p "//output_directory)
    do point = 1, point_count
        coordinate(point) = real(point - 1, dp)/real(point_count - 1, dp)
        angle = 0.5_dp*pi*coordinate(point)
        direction = [cos(angle), sin(angle), 0.0_dp]
        gradient = [sin(2.0_dp*pi*coordinate(point)), &
            cos(2.0_dp*pi*coordinate(point)), 0.0_dp]
        call evaluate_field_aligned_flux( &
            parallel_coefficient, perpendicular_coefficient, direction, gradient, &
            flux, status)
        if (status%code /= 0) error stop "field-aligned flux evaluation failed"
        parallel_component(point) = dot_product(direction, gradient)
        perpendicular_component(point) = sqrt(max(0.0_dp, &
            dot_product(gradient, gradient) - parallel_component(point)**2))
        flux_norm(point) = sqrt(dot_product(flux, flux))
    end do

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, parallel_component, label="parallel gradient")
    call plot(coordinate, perpendicular_component, label="perpendicular gradient")
    call xlabel("normalized coordinate")
    call ylabel("gradient component")
    call title("Field-aligned anisotropy: manufactured gradient")
    call legend()
    call savefig(output_directory//"/field_aligned_flux_components_1d.png")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(coordinate, flux_norm, label="|K grad u|")
    call xlabel("normalized coordinate")
    call ylabel("flux magnitude")
    call title("Field-aligned anisotropy: k_parallel/k_perpendicular = 100")
    call legend()
    call savefig(output_directory//"/field_aligned_flux_1d.png")

    open (newunit=unit, file=output_directory//"/field_aligned_flux.csv", &
        status="replace", action="write")
    write (unit, "(a)") "coordinate,parallel_gradient,perpendicular_gradient,flux_norm"
    do point = 1, point_count
        write (unit, "(4(es24.16,1x))") coordinate(point), &
            parallel_component(point), perpendicular_component(point), &
            flux_norm(point)
    end do
    close (unit)
end program field_aligned_flux
