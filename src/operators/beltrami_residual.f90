module fortfem_beltrami_residual
    !! Neutral multi-region Beltrami/curl-eigen residual composition.
    !!
    !! For region r, sample q, and component c, the primary residual is
    !!
    !!     R(r,q,c) = curl_B(r,q,c) - lambda(r) B(r,q,c).
    !!
    !! The caller may append independently supplied divergence, flux, and
    !! helicity rows.  Each row family is represented by a value and target
    !! array; a zero-sized pair disables that family.  This module deliberately
    !! does not compute curl, divergence, flux, helicity, constitutive laws, or
    !! equilibrium constraints.  Those operations remain in the consuming
    !! spatial/physics layer.  Exact real JVP and VJP actions make this contract
    !! usable by compatible H(curl), IGA, Fourier, and multi-region clients.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_beltrami_residual
    public :: assemble_beltrami_residual_jvp
    public :: assemble_beltrami_residual_vjp

contains

    subroutine assemble_beltrami_residual( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, &
            flux, flux_target, helicity, helicity_target, residual, &
            divergence_residual, flux_residual, helicity_residual, status)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :)
        real(dp), intent(in) :: lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(out) :: residual(:, :, :), divergence_residual(:, :)
        real(dp), intent(out) :: flux_residual(:), helicity_residual(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: region, sample, component

        residual = 0.0_dp
        divergence_residual = 0.0_dp
        flux_residual = 0.0_dp
        helicity_residual = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami residual received incompatible arrays")
        if (.not. valid_value_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual, divergence_residual, &
            flux_residual, helicity_residual)) return
        if (.not. finite_value_inputs( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target)) return

        do region = 1, size(curl_field, 1)
            do sample = 1, size(curl_field, 2)
                do component = 1, size(curl_field, 3)
                    residual(region, sample, component) = &
                        curl_field(region, sample, component) - &
                        lambda(region)*magnetic_field(region, sample, component)
                end do
            end do
        end do
        divergence_residual = divergence - divergence_target
        flux_residual = flux - flux_target
        helicity_residual = helicity - helicity_target
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_beltrami_residual

    subroutine assemble_beltrami_residual_jvp( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, curl_field_dot, &
            magnetic_field_dot, lambda_dot, divergence_dot, divergence_target_dot, &
            flux_dot, flux_target_dot, helicity_dot, helicity_target_dot, residual_dot, &
            divergence_residual_dot, flux_residual_dot, helicity_residual_dot, status)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :)
        real(dp), intent(in) :: lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: curl_field_dot(:, :, :), magnetic_field_dot(:, :, :)
        real(dp), intent(in) :: lambda_dot(:), divergence_dot(:, :)
        real(dp), intent(in) :: divergence_target_dot(:, :)
        real(dp), intent(in) :: flux_dot(:), flux_target_dot(:), helicity_dot(:)
        real(dp), intent(in) :: helicity_target_dot(:)
        real(dp), intent(out) :: residual_dot(:, :, :), divergence_residual_dot(:, :)
        real(dp), intent(out) :: flux_residual_dot(:), helicity_residual_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: region, sample, component

        residual_dot = 0.0_dp
        divergence_residual_dot = 0.0_dp
        flux_residual_dot = 0.0_dp
        helicity_residual_dot = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami residual JVP received incompatible arrays")
        if (.not. valid_value_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual_dot, &
            divergence_residual_dot, flux_residual_dot, helicity_residual_dot)) return
        if (.not. valid_direction_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, curl_field_dot, &
            magnetic_field_dot, lambda_dot, divergence_dot, divergence_target_dot, &
            flux_dot, flux_target_dot, helicity_dot, helicity_target_dot)) return
        if (.not. finite_value_inputs( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target) .or. &
            .not. finite_direction_inputs( &
            curl_field_dot, magnetic_field_dot, lambda_dot, divergence_dot, &
            divergence_target_dot, flux_dot, flux_target_dot, helicity_dot, &
            helicity_target_dot)) return

        do region = 1, size(curl_field, 1)
            do sample = 1, size(curl_field, 2)
                do component = 1, size(curl_field, 3)
                    residual_dot(region, sample, component) = &
                        curl_field_dot(region, sample, component) - &
                        lambda_dot(region)*magnetic_field(region, sample, component) - &
                        lambda(region)*magnetic_field_dot(region, sample, component)
                end do
            end do
        end do
        divergence_residual_dot = divergence_dot - divergence_target_dot
        flux_residual_dot = flux_dot - flux_target_dot
        helicity_residual_dot = helicity_dot - helicity_target_dot
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_beltrami_residual_jvp

    subroutine assemble_beltrami_residual_vjp( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual_bar, &
            divergence_residual_bar, flux_residual_bar, helicity_residual_bar, &
            curl_field_bar, magnetic_field_bar, lambda_bar, divergence_bar, &
            divergence_target_bar, flux_bar, flux_target_bar, helicity_bar, &
            helicity_target_bar, status)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :)
        real(dp), intent(in) :: lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: residual_bar(:, :, :), divergence_residual_bar(:, :)
        real(dp), intent(in) :: flux_residual_bar(:), helicity_residual_bar(:)
        real(dp), intent(out) :: curl_field_bar(:, :, :), magnetic_field_bar(:, :, :)
        real(dp), intent(out) :: lambda_bar(:), divergence_bar(:, :)
        real(dp), intent(out) :: divergence_target_bar(:, :), flux_bar(:)
        real(dp), intent(out) :: flux_target_bar(:), helicity_bar(:)
        real(dp), intent(out) :: helicity_target_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: region, sample, component

        curl_field_bar = 0.0_dp
        magnetic_field_bar = 0.0_dp
        lambda_bar = 0.0_dp
        divergence_bar = 0.0_dp
        divergence_target_bar = 0.0_dp
        flux_bar = 0.0_dp
        flux_target_bar = 0.0_dp
        helicity_bar = 0.0_dp
        helicity_target_bar = 0.0_dp
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "Beltrami residual VJP received incompatible arrays")
        if (.not. valid_value_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual_bar, divergence_residual_bar, &
            flux_residual_bar, helicity_residual_bar)) return
        if (.not. valid_vjp_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual_bar, &
            divergence_residual_bar, flux_residual_bar, helicity_residual_bar, &
            curl_field_bar, magnetic_field_bar, lambda_bar, divergence_bar, &
            divergence_target_bar, flux_bar, flux_target_bar, helicity_bar, &
            helicity_target_bar)) return
        if (.not. finite_value_inputs( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target) .or. &
            .not. all(ieee_is_finite(residual_bar)) .or. &
            .not. all(ieee_is_finite(divergence_residual_bar)) .or. &
            .not. all(ieee_is_finite(flux_residual_bar)) .or. &
            .not. all(ieee_is_finite(helicity_residual_bar))) return

        curl_field_bar = residual_bar
        do region = 1, size(curl_field, 1)
            do sample = 1, size(curl_field, 2)
                do component = 1, size(curl_field, 3)
                    magnetic_field_bar(region, sample, component) = &
                        -lambda(region)*residual_bar(region, sample, component)
                    lambda_bar(region) = lambda_bar(region) - &
                        residual_bar(region, sample, component)* &
                        magnetic_field(region, sample, component)
                end do
            end do
        end do
        divergence_bar = divergence_residual_bar
        divergence_target_bar = -divergence_residual_bar
        flux_bar = flux_residual_bar
        flux_target_bar = -flux_residual_bar
        helicity_bar = helicity_residual_bar
        helicity_target_bar = -helicity_residual_bar
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine assemble_beltrami_residual_vjp

    logical function valid_value_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual, divergence_residual, &
            flux_residual, helicity_residual) result(valid)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: residual(:, :, :), divergence_residual(:, :)
        real(dp), intent(in) :: flux_residual(:), helicity_residual(:)

        valid = size(curl_field, 1) > 0 .and. size(curl_field, 2) > 0 .and. &
            size(curl_field, 3) > 0 .and. all(shape(magnetic_field) == shape(curl_field)) .and. &
            size(lambda) == size(curl_field, 1) .and. &
            all(shape(residual) == shape(curl_field)) .and. &
            all(shape(divergence_target) == shape(divergence)) .and. &
            all(shape(divergence_residual) == shape(divergence)) .and. &
            size(flux_target) == size(flux) .and. size(flux_residual) == size(flux) .and. &
            size(helicity_target) == size(helicity) .and. &
            size(helicity_residual) == size(helicity)
    end function valid_value_shapes

    logical function valid_direction_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, curl_field_dot, magnetic_field_dot, &
            lambda_dot, divergence_dot, divergence_target_dot, flux_dot, flux_target_dot, &
            helicity_dot, helicity_target_dot) result(valid)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: curl_field_dot(:, :, :), magnetic_field_dot(:, :, :)
        real(dp), intent(in) :: lambda_dot(:), divergence_dot(:, :)
        real(dp), intent(in) :: divergence_target_dot(:, :), flux_dot(:), flux_target_dot(:)
        real(dp), intent(in) :: helicity_dot(:), helicity_target_dot(:)

        valid = all(shape(curl_field_dot) == shape(curl_field)) .and. &
            all(shape(magnetic_field_dot) == shape(magnetic_field)) .and. &
            size(lambda_dot) == size(lambda) .and. &
            all(shape(divergence_dot) == shape(divergence)) .and. &
            all(shape(divergence_target_dot) == shape(divergence_target)) .and. &
            size(flux_dot) == size(flux) .and. size(flux_target_dot) == size(flux_target) .and. &
            size(helicity_dot) == size(helicity) .and. &
            size(helicity_target_dot) == size(helicity_target)
    end function valid_direction_shapes

    logical function valid_vjp_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual_bar, &
            divergence_residual_bar, flux_residual_bar, helicity_residual_bar, &
            curl_field_bar, magnetic_field_bar, lambda_bar, divergence_bar, &
            divergence_target_bar, flux_bar, flux_target_bar, helicity_bar, &
            helicity_target_bar) result(valid)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)
        real(dp), intent(in) :: residual_bar(:, :, :), divergence_residual_bar(:, :)
        real(dp), intent(in) :: flux_residual_bar(:), helicity_residual_bar(:)
        real(dp), intent(in) :: curl_field_bar(:, :, :), magnetic_field_bar(:, :, :)
        real(dp), intent(in) :: lambda_bar(:), divergence_bar(:, :)
        real(dp), intent(in) :: divergence_target_bar(:, :), flux_bar(:)
        real(dp), intent(in) :: flux_target_bar(:), helicity_bar(:)
        real(dp), intent(in) :: helicity_target_bar(:)

        valid = valid_value_shapes( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target, residual_bar, &
            divergence_residual_bar, flux_residual_bar, helicity_residual_bar) .and. &
            all(shape(curl_field_bar) == shape(curl_field)) .and. &
            all(shape(magnetic_field_bar) == shape(magnetic_field)) .and. &
            size(lambda_bar) == size(lambda) .and. &
            all(shape(divergence_bar) == shape(divergence)) .and. &
            all(shape(divergence_target_bar) == shape(divergence_target)) .and. &
            size(flux_bar) == size(flux) .and. size(flux_target_bar) == size(flux_target) .and. &
            size(helicity_bar) == size(helicity) .and. &
            size(helicity_target_bar) == size(helicity_target)
    end function valid_vjp_shapes

    logical function finite_value_inputs( &
            curl_field, magnetic_field, lambda, divergence, divergence_target, flux, &
            flux_target, helicity, helicity_target) result(valid)
        real(dp), intent(in) :: curl_field(:, :, :), magnetic_field(:, :, :), lambda(:)
        real(dp), intent(in) :: divergence(:, :), divergence_target(:, :)
        real(dp), intent(in) :: flux(:), flux_target(:), helicity(:), helicity_target(:)

        valid = all(ieee_is_finite(curl_field)) .and. &
            all(ieee_is_finite(magnetic_field)) .and. all(ieee_is_finite(lambda)) .and. &
            all(ieee_is_finite(divergence)) .and. all(ieee_is_finite(divergence_target)) .and. &
            all(ieee_is_finite(flux)) .and. all(ieee_is_finite(flux_target)) .and. &
            all(ieee_is_finite(helicity)) .and. all(ieee_is_finite(helicity_target))
    end function finite_value_inputs

    logical function finite_direction_inputs( &
            curl_field_dot, magnetic_field_dot, lambda_dot, divergence_dot, &
            divergence_target_dot, flux_dot, flux_target_dot, helicity_dot, &
            helicity_target_dot) result(valid)
        real(dp), intent(in) :: curl_field_dot(:, :, :), magnetic_field_dot(:, :, :)
        real(dp), intent(in) :: lambda_dot(:), divergence_dot(:, :)
        real(dp), intent(in) :: divergence_target_dot(:, :), flux_dot(:), flux_target_dot(:)
        real(dp), intent(in) :: helicity_dot(:), helicity_target_dot(:)

        valid = all(ieee_is_finite(curl_field_dot)) .and. &
            all(ieee_is_finite(magnetic_field_dot)) .and. &
            all(ieee_is_finite(lambda_dot)) .and. all(ieee_is_finite(divergence_dot)) .and. &
            all(ieee_is_finite(divergence_target_dot)) .and. all(ieee_is_finite(flux_dot)) .and. &
            all(ieee_is_finite(flux_target_dot)) .and. all(ieee_is_finite(helicity_dot)) .and. &
            all(ieee_is_finite(helicity_target_dot))
    end function finite_direction_inputs

end module fortfem_beltrami_residual
