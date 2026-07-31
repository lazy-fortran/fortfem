module fortfem_field_aligned_flux
    !! Generated strongly anisotropic constitutive flux.
    !!
    !! For a unit field direction b and gradient g,
    !!
    !!   F = k_perpendicular g
    !!       + (k_parallel-k_perpendicular) b (b . g).
    !!
    !! FortSym owns the scalar product and its derivatives; this wrapper owns
    !! vector layout, non-negative coefficient, and unit-direction contracts.
    use fortfem_generated_field_aligned_flux, only: generated_field_aligned_flux
    use fortfem_generated_field_aligned_flux_jvp, only: &
        generated_field_aligned_flux_jvp
    use fortfem_generated_field_aligned_flux_vjp, only: &
        generated_field_aligned_flux_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    real(dp), parameter :: unit_tolerance = 1.0e-12_dp

    public :: evaluate_field_aligned_flux
    public :: evaluate_field_aligned_flux_jvp
    public :: evaluate_field_aligned_flux_vjp

contains

    pure subroutine evaluate_field_aligned_flux( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, flux, status)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), gradient(:)
        real(dp), intent(out) :: flux(:)
        type(fortsparse_status_t), intent(out) :: status

        flux = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned flux received incompatible arrays")
        if (.not. valid_inputs( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, flux)) return
        call generated_field_aligned_flux( &
            parallel_coefficient, perpendicular_coefficient, unit_direction(1), &
            unit_direction(2), unit_direction(3), gradient(1), gradient(2), &
            gradient(3), flux(1), flux(2), flux(3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_field_aligned_flux

    pure subroutine evaluate_field_aligned_flux_jvp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, parallel_coefficient_dot, perpendicular_coefficient_dot, &
            direction_dot, gradient_dot, flux_dot, status)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), gradient(:)
        real(dp), intent(in) :: parallel_coefficient_dot
        real(dp), intent(in) :: perpendicular_coefficient_dot
        real(dp), intent(in) :: direction_dot(:), gradient_dot(:)
        real(dp), intent(out) :: flux_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        flux_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned flux JVP received incompatible arrays")
        if (.not. valid_inputs( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, flux_dot)) return
        if (size(direction_dot) /= 3 .or. size(gradient_dot) /= 3) return
        call generated_field_aligned_flux_jvp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction(1), &
            unit_direction(2), unit_direction(3), gradient(1), gradient(2), &
            gradient(3), parallel_coefficient_dot, perpendicular_coefficient_dot, &
            direction_dot(1), direction_dot(2), direction_dot(3), gradient_dot(1), &
            gradient_dot(2), gradient_dot(3), flux_dot(1), flux_dot(2), flux_dot(3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_field_aligned_flux_jvp

    pure subroutine evaluate_field_aligned_flux_vjp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, flux_bar, parallel_coefficient_bar, &
            perpendicular_coefficient_bar, direction_bar, gradient_bar, status)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), gradient(:), flux_bar(:)
        real(dp), intent(out) :: parallel_coefficient_bar
        real(dp), intent(out) :: perpendicular_coefficient_bar
        real(dp), intent(out) :: direction_bar(:), gradient_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        parallel_coefficient_bar = 0.0_dp
        perpendicular_coefficient_bar = 0.0_dp
        direction_bar = 0.0_dp
        gradient_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "field-aligned flux VJP received incompatible arrays")
        if (.not. valid_inputs( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, flux_bar)) return
        if (size(direction_bar) /= 3 .or. size(gradient_bar) /= 3) return
        call generated_field_aligned_flux_vjp( &
            parallel_coefficient, perpendicular_coefficient, unit_direction(1), &
            unit_direction(2), unit_direction(3), gradient(1), gradient(2), &
            gradient(3), flux_bar(1), flux_bar(2), flux_bar(3), &
            parallel_coefficient_bar, perpendicular_coefficient_bar, &
            direction_bar(1), direction_bar(2), direction_bar(3), gradient_bar(1), &
            gradient_bar(2), gradient_bar(3))
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_field_aligned_flux_vjp

    pure logical function valid_inputs( &
            parallel_coefficient, perpendicular_coefficient, unit_direction, &
            gradient, flux)
        real(dp), intent(in) :: parallel_coefficient, perpendicular_coefficient
        real(dp), intent(in) :: unit_direction(:), gradient(:), flux(:)

        valid_inputs = parallel_coefficient >= 0.0_dp .and. &
            perpendicular_coefficient >= 0.0_dp .and. &
            size(unit_direction) == 3 .and. size(gradient) == 3 .and. &
            size(flux) == 3 .and. abs(dot_product(unit_direction, &
            unit_direction) - 1.0_dp) <= unit_tolerance
    end function valid_inputs

end module fortfem_field_aligned_flux
