module fortfem_nested_surface_geometry
    !! A neutral fixed-topology nested-surface embedding.
    !!
    !! The caller supplies three real physical fields, represented by complex
    !! Fourier coefficients over the radial/mode registry: R, Z, and lambda.
    !! At fixed (rho, theta, zeta), the inverse coordinates are mapped to
    !! Cartesian space with phi = zeta + lambda.  No equilibrium equation,
    !! profile, flux normalization, or Boozer convention is selected here.
    !! The routine is therefore usable by external VMEC, GVEC, DESC, SPEC,
    !! or free-boundary adapters without importing their physics or formats.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_fourier_mode_registry, only: &
        evaluate_fourier_mode, fourier_mode_registry_t, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_nested_surface_geometry
    public :: evaluate_nested_surface_geometry_jvp
    public :: evaluate_nested_surface_geometry_vjp
    public :: evaluate_nested_surface_geometry_coordinate_jvp
    public :: evaluate_nested_surface_geometry_coordinate_vjp

    interface finite_real
        module procedure finite_real_1d
        module procedure finite_real_2d
        module procedure finite_real_3d
    end interface finite_real

contains

    subroutine evaluate_nested_surface_geometry( &
            registry, coefficients, rho, theta, zeta, mapped_coordinates, &
            physical_coordinates, mapped_jacobian, physical_jacobian, metric, &
            volume_jacobian, status)
        !! Evaluate inverse and Cartesian coordinates and their diagnostics.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        real(dp), intent(out) :: mapped_coordinates(:, :)
        real(dp), intent(out) :: physical_coordinates(:, :)
        real(dp), intent(out) :: mapped_jacobian(:, :, :)
        real(dp), intent(out) :: physical_jacobian(:, :, :)
        real(dp), intent(out) :: metric(:, :, :)
        real(dp), intent(out) :: volume_jacobian(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample
        real(dp) :: mapped(3), mapped_jac(3, 3)

        mapped_coordinates = 0.0_dp
        physical_coordinates = 0.0_dp
        mapped_jacobian = 0.0_dp
        physical_jacobian = 0.0_dp
        metric = 0.0_dp
        volume_jacobian = 0.0_dp
        call validate_base_inputs( &
            registry, coefficients, rho, theta, zeta, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_output_shapes( &
                rho, mapped_coordinates, physical_coordinates, &
                mapped_jacobian, physical_jacobian, metric, volume_jacobian)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface evaluator has incompatible output shapes")
            return
        end if

        do sample = 1, size(rho)
            call evaluate_modal_fields( &
                registry, coefficients, rho(sample), theta(sample), &
                zeta(sample), mapped, mapped_jac, status)
            if (status%code /= FORTSPARSE_OK) return
            mapped_coordinates(:, sample) = mapped
            call map_inverse_coordinates( &
                mapped, mapped_jac, zeta(sample), &
                physical_coordinates(:, sample), physical_jacobian(:, :, sample), &
                metric(:, :, sample), volume_jacobian(sample))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_surface_geometry

    subroutine evaluate_nested_surface_geometry_coordinate_jvp( &
            registry, coefficients, rho, theta, zeta, rho_dot, theta_dot, &
            zeta_dot, mapped_coordinates_dot, physical_coordinates_dot, status)
        !! Fixed-topology tangent for perturbations of sample coordinates.
        !!
        !! This action differentiates the returned coordinate values only.  The
        !! Jacobian, metric, and volume diagnostics have coordinate derivatives
        !! of second order and remain separate planned products.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        real(dp), intent(in) :: rho_dot(:), theta_dot(:), zeta_dot(:)
        real(dp), intent(out) :: mapped_coordinates_dot(:, :)
        real(dp), intent(out) :: physical_coordinates_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample
        real(dp) :: mapped(3), mapped_jac(3, 3)
        real(dp) :: physical(3), physical_jac(3, 3), metric(3, 3), volume
        real(dp) :: coordinate_dot(3)

        mapped_coordinates_dot = 0.0_dp
        physical_coordinates_dot = 0.0_dp
        call validate_base_inputs( &
            registry, coefficients, rho, theta, zeta, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_coordinate_jvp_shapes( &
                rho, rho_dot, theta_dot, zeta_dot, mapped_coordinates_dot, &
                physical_coordinates_dot) .or. &
                .not. finite_real(rho_dot) .or. &
                .not. finite_real(theta_dot) .or. &
                .not. finite_real(zeta_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface coordinate JVP has invalid shapes or values")
            return
        end if

        do sample = 1, size(rho)
            call evaluate_modal_fields( &
                registry, coefficients, rho(sample), theta(sample), &
                zeta(sample), mapped, mapped_jac, status)
            if (status%code /= FORTSPARSE_OK) return
            call map_inverse_coordinates( &
                mapped, mapped_jac, zeta(sample), physical, physical_jac, &
                metric, volume)
            coordinate_dot = [rho_dot(sample), theta_dot(sample), zeta_dot(sample)]
            mapped_coordinates_dot(:, sample) = matmul(mapped_jac, coordinate_dot)
            physical_coordinates_dot(:, sample) = &
                matmul(physical_jac, coordinate_dot)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_surface_geometry_coordinate_jvp

    subroutine evaluate_nested_surface_geometry_coordinate_vjp( &
            registry, coefficients, rho, theta, zeta, mapped_coordinates_bar, &
            physical_coordinates_bar, rho_bar, theta_bar, zeta_bar, status)
        !! Reverse fixed-topology action for sample-coordinate perturbations.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        real(dp), intent(in) :: mapped_coordinates_bar(:, :)
        real(dp), intent(in) :: physical_coordinates_bar(:, :)
        real(dp), intent(out) :: rho_bar(:), theta_bar(:), zeta_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample
        real(dp) :: mapped(3), mapped_jac(3, 3)
        real(dp) :: physical(3), physical_jac(3, 3), metric(3, 3), volume
        real(dp) :: coordinate_bar(3)

        rho_bar = 0.0_dp
        theta_bar = 0.0_dp
        zeta_bar = 0.0_dp
        call validate_base_inputs( &
            registry, coefficients, rho, theta, zeta, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_coordinate_vjp_shapes( &
                rho, mapped_coordinates_bar, physical_coordinates_bar, &
                rho_bar, theta_bar, zeta_bar) .or. &
                .not. finite_real(mapped_coordinates_bar) .or. &
                .not. finite_real(physical_coordinates_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface coordinate VJP has invalid shapes or values")
            return
        end if

        do sample = 1, size(rho)
            call evaluate_modal_fields( &
                registry, coefficients, rho(sample), theta(sample), &
                zeta(sample), mapped, mapped_jac, status)
            if (status%code /= FORTSPARSE_OK) return
            call map_inverse_coordinates( &
                mapped, mapped_jac, zeta(sample), physical, physical_jac, &
                metric, volume)
            coordinate_bar = matmul(transpose(mapped_jac), &
                mapped_coordinates_bar(:, sample)) + &
                matmul(transpose(physical_jac), physical_coordinates_bar(:, sample))
            rho_bar(sample) = coordinate_bar(1)
            theta_bar(sample) = coordinate_bar(2)
            zeta_bar(sample) = coordinate_bar(3)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_surface_geometry_coordinate_vjp

    subroutine evaluate_nested_surface_geometry_jvp( &
            registry, coefficients, rho, theta, zeta, coefficients_dot, &
            mapped_coordinates_dot, physical_coordinates_dot, &
            mapped_jacobian_dot, physical_jacobian_dot, metric_dot, &
            volume_jacobian_dot, status)
        !! Coefficient-only fixed-topology tangent of the full geometry map.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        complex(dp), intent(in) :: coefficients_dot(:, :)
        real(dp), intent(out) :: mapped_coordinates_dot(:, :)
        real(dp), intent(out) :: physical_coordinates_dot(:, :)
        real(dp), intent(out) :: mapped_jacobian_dot(:, :, :)
        real(dp), intent(out) :: physical_jacobian_dot(:, :, :)
        real(dp), intent(out) :: metric_dot(:, :, :)
        real(dp), intent(out) :: volume_jacobian_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample
        real(dp) :: mapped(3), mapped_jac(3, 3)
        real(dp) :: mapped_dot(3), mapped_jac_dot(3, 3)
        real(dp) :: physical(3), physical_jac(3, 3), metric(3, 3)
        real(dp) :: volume

        mapped_coordinates_dot = 0.0_dp
        physical_coordinates_dot = 0.0_dp
        mapped_jacobian_dot = 0.0_dp
        physical_jacobian_dot = 0.0_dp
        metric_dot = 0.0_dp
        volume_jacobian_dot = 0.0_dp
        call validate_base_inputs( &
            registry, coefficients, rho, theta, zeta, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_output_shapes( &
                rho, mapped_coordinates_dot, physical_coordinates_dot, &
                mapped_jacobian_dot, physical_jacobian_dot, metric_dot, &
                volume_jacobian_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface JVP has incompatible output shapes")
            return
        end if
        if (size(coefficients_dot, 1) /= 3 .or. &
                size(coefficients_dot, 2) /= size(coefficients, 2) .or. &
                .not. finite_complex(coefficients_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface JVP has invalid coefficient increments")
            return
        end if

        do sample = 1, size(rho)
            call evaluate_modal_fields( &
                registry, coefficients, rho(sample), theta(sample), &
                zeta(sample), mapped, mapped_jac, status)
            if (status%code /= FORTSPARSE_OK) return
            call evaluate_modal_fields( &
                registry, coefficients_dot, rho(sample), theta(sample), &
                zeta(sample), mapped_dot, mapped_jac_dot, status)
            if (status%code /= FORTSPARSE_OK) return
            call map_inverse_coordinates( &
                mapped, mapped_jac, zeta(sample), physical, physical_jac, &
                metric, volume)
            call map_inverse_coordinates_tangent( &
                mapped, mapped_jac, mapped_dot, mapped_jac_dot, zeta(sample), &
                physical_coordinates_dot(:, sample), &
                physical_jacobian_dot(:, :, sample), &
                metric_dot(:, :, sample), volume_jacobian_dot(sample))
            mapped_coordinates_dot(:, sample) = mapped_dot
            mapped_jacobian_dot(:, :, sample) = mapped_jac_dot
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_surface_geometry_jvp

    subroutine evaluate_nested_surface_geometry_vjp( &
            registry, coefficients, rho, theta, zeta, mapped_coordinates_bar, &
            physical_coordinates_bar, mapped_jacobian_bar, physical_jacobian_bar, &
            metric_bar, volume_jacobian_bar, coefficients_bar, status)
        !! Reverse coefficient action for all returned geometry quantities.
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        real(dp), intent(in) :: mapped_coordinates_bar(:, :)
        real(dp), intent(in) :: physical_coordinates_bar(:, :)
        real(dp), intent(in) :: mapped_jacobian_bar(:, :, :)
        real(dp), intent(in) :: physical_jacobian_bar(:, :, :)
        real(dp), intent(in) :: metric_bar(:, :, :)
        real(dp), intent(in) :: volume_jacobian_bar(:)
        complex(dp), intent(out) :: coefficients_bar(:, :)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, mode, component, coordinate
        complex(dp) :: value, radial_derivative, theta_derivative
        complex(dp) :: phi_derivative
        real(dp) :: mapped(3), mapped_jac(3, 3), physical(3)
        real(dp) :: physical_jac(3, 3), metric(3, 3), volume
        real(dp) :: mapped_bar(3), mapped_jac_bar_local(3, 3)

        coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        call validate_base_inputs( &
            registry, coefficients, rho, theta, zeta, status)
        if (status%code /= FORTSPARSE_OK) return
        if (.not. valid_vjp_shapes( &
                rho, mapped_coordinates_bar, physical_coordinates_bar, &
                mapped_jacobian_bar, physical_jacobian_bar, metric_bar, &
                volume_jacobian_bar, coefficients_bar, &
                size(registry%poloidal_modes)) .or. &
                .not. finite_real(mapped_coordinates_bar) .or. &
                .not. finite_real(physical_coordinates_bar) .or. &
                .not. finite_real(mapped_jacobian_bar) .or. &
                .not. finite_real(physical_jacobian_bar) .or. &
                .not. finite_real(metric_bar) .or. &
                .not. finite_real(volume_jacobian_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface VJP has invalid cotangent shapes or values")
            return
        end if

        do sample = 1, size(rho)
            call evaluate_modal_fields( &
                registry, coefficients, rho(sample), theta(sample), &
                zeta(sample), mapped, mapped_jac, status)
            if (status%code /= FORTSPARSE_OK) return
            call map_inverse_coordinates( &
                mapped, mapped_jac, zeta(sample), physical, physical_jac, &
                metric, volume)
            mapped_bar = mapped_coordinates_bar(:, sample)
            mapped_jac_bar_local = mapped_jacobian_bar(:, :, sample)
            call reverse_inverse_map( &
                mapped, mapped_jac, zeta(sample), physical_coordinates_bar(:, sample), &
                physical_jacobian_bar(:, :, sample), metric_bar(:, :, sample), &
                volume_jacobian_bar(sample), mapped_bar, mapped_jac_bar_local)

            do mode = 1, size(coefficients, 2)
                call evaluate_fourier_mode( &
                    registry, mode, rho(sample), theta(sample), zeta(sample), &
                    value, radial_derivative, theta_derivative, phi_derivative, status)
                if (status%code /= FORTSPARSE_OK) return
                do component = 1, 3
                    coefficients_bar(component, mode) = &
                        coefficients_bar(component, mode) + &
                        conjg(value)*mapped_bar(component)
                    coefficients_bar(component, mode) = &
                        coefficients_bar(component, mode) + &
                        conjg(radial_derivative)*mapped_jac_bar_local(component, 1)
                    coefficients_bar(component, mode) = &
                        coefficients_bar(component, mode) + &
                        conjg(theta_derivative)*mapped_jac_bar_local(component, 2)
                    coefficients_bar(component, mode) = &
                        coefficients_bar(component, mode) + &
                        conjg(phi_derivative)*mapped_jac_bar_local(component, 3)
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_nested_surface_geometry_vjp

    subroutine validate_base_inputs( &
            registry, coefficients, rho, theta, zeta, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho(:), theta(:), zeta(:)
        type(fortsparse_status_t), intent(out) :: status

        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (size(coefficients, 1) /= 3 .or. &
                size(coefficients, 2) /= size(registry%poloidal_modes) .or. &
                size(rho) < 1 .or. size(theta) /= size(rho) .or. &
                size(zeta) /= size(rho)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface inputs have incompatible dimensions")
            return
        end if
        if (.not. finite_complex(coefficients) .or. &
                .not. finite_real(rho) .or. .not. finite_real(theta) .or. &
                .not. finite_real(zeta) .or. any(rho < 0.0_dp)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "nested-surface inputs contain invalid values")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_base_inputs

    logical function valid_output_shapes( &
            rho, mapped, physical, mapped_jac, physical_jac, metric, volume) &
            result(valid)
        real(dp), intent(in) :: rho(:), mapped(:, :), physical(:, :)
        real(dp), intent(in) :: mapped_jac(:, :, :), physical_jac(:, :, :)
        real(dp), intent(in) :: metric(:, :, :), volume(:)
        integer :: sample_count

        sample_count = size(rho)
        valid = size(mapped, 1) == 3 .and. size(mapped, 2) == sample_count .and. &
            size(physical, 1) == 3 .and. size(physical, 2) == sample_count .and. &
            size(mapped_jac, 1) == 3 .and. size(mapped_jac, 2) == 3 .and. &
            size(mapped_jac, 3) == sample_count .and. &
            size(physical_jac, 1) == 3 .and. size(physical_jac, 2) == 3 .and. &
            size(physical_jac, 3) == sample_count .and. &
            size(metric, 1) == 3 .and. size(metric, 2) == 3 .and. &
            size(metric, 3) == sample_count .and. size(volume) == sample_count
    end function valid_output_shapes

    logical function valid_vjp_shapes( &
            rho, mapped, physical, mapped_jac, physical_jac, metric, volume, &
            coefficients, mode_count) result(valid)
        real(dp), intent(in) :: rho(:), mapped(:, :), physical(:, :)
        real(dp), intent(in) :: mapped_jac(:, :, :), physical_jac(:, :, :)
        real(dp), intent(in) :: metric(:, :, :), volume(:)
        complex(dp), intent(in) :: coefficients(:, :)
        integer, intent(in) :: mode_count

        valid = valid_output_shapes( &
            rho, mapped, physical, mapped_jac, physical_jac, metric, volume) .and. &
            size(coefficients, 1) == 3 .and. size(coefficients, 2) == mode_count
    end function valid_vjp_shapes

    logical function valid_coordinate_jvp_shapes( &
            rho, rho_dot, theta_dot, zeta_dot, mapped, physical) result(valid)
        real(dp), intent(in) :: rho(:), rho_dot(:), theta_dot(:), zeta_dot(:)
        real(dp), intent(in) :: mapped(:, :), physical(:, :)

        valid = size(rho_dot) == size(rho) .and. &
            size(theta_dot) == size(rho) .and. size(zeta_dot) == size(rho) .and. &
            size(mapped, 1) == 3 .and. size(mapped, 2) == size(rho) .and. &
            size(physical, 1) == 3 .and. size(physical, 2) == size(rho)
    end function valid_coordinate_jvp_shapes

    logical function valid_coordinate_vjp_shapes( &
            rho, mapped, physical, rho_bar, theta_bar, zeta_bar) result(valid)
        real(dp), intent(in) :: rho(:), mapped(:, :), physical(:, :)
        real(dp), intent(in) :: rho_bar(:), theta_bar(:), zeta_bar(:)

        valid = size(mapped, 1) == 3 .and. size(mapped, 2) == size(rho) .and. &
            size(physical, 1) == 3 .and. size(physical, 2) == size(rho) .and. &
            size(rho_bar) == size(rho) .and. size(theta_bar) == size(rho) .and. &
            size(zeta_bar) == size(rho)
    end function valid_coordinate_vjp_shapes

    subroutine evaluate_modal_fields( &
            registry, coefficients, rho, theta, zeta, fields, field_jac, status)
        type(fourier_mode_registry_t), intent(in) :: registry
        complex(dp), intent(in) :: coefficients(:, :)
        real(dp), intent(in) :: rho, theta, zeta
        real(dp), intent(out) :: fields(3), field_jac(3, 3)
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, component
        complex(dp) :: value, radial_derivative, theta_derivative
        complex(dp) :: phi_derivative

        fields = 0.0_dp
        field_jac = 0.0_dp
        do mode = 1, size(coefficients, 2)
            call evaluate_fourier_mode( &
                registry, mode, rho, theta, zeta, value, radial_derivative, &
                theta_derivative, phi_derivative, status)
            if (status%code /= FORTSPARSE_OK) return
            do component = 1, 3
                fields(component) = fields(component) + &
                    real(coefficients(component, mode)*value, dp)
                field_jac(component, 1) = field_jac(component, 1) + &
                    real(coefficients(component, mode)*radial_derivative, dp)
                field_jac(component, 2) = field_jac(component, 2) + &
                    real(coefficients(component, mode)*theta_derivative, dp)
                field_jac(component, 3) = field_jac(component, 3) + &
                    real(coefficients(component, mode)*phi_derivative, dp)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_modal_fields

    subroutine map_inverse_coordinates( &
            mapped, mapped_jac, zeta, physical, physical_jac, metric, volume)
        real(dp), intent(in) :: mapped(3), mapped_jac(3, 3), zeta
        real(dp), intent(out) :: physical(3), physical_jac(3, 3)
        real(dp), intent(out) :: metric(3, 3), volume

        real(dp) :: radius, height, lambda, phi, cosine, sine, dphi
        integer :: coordinate, first, second

        radius = mapped(1)
        height = mapped(2)
        lambda = mapped(3)
        phi = zeta + lambda
        cosine = cos(phi)
        sine = sin(phi)
        physical = [radius*cosine, radius*sine, height]
        do coordinate = 1, 3
            dphi = mapped_jac(3, coordinate)
            if (coordinate == 3) dphi = dphi + 1.0_dp
            physical_jac(:, coordinate) = [ &
                cosine*mapped_jac(1, coordinate) - radius*sine*dphi, &
                sine*mapped_jac(1, coordinate) + radius*cosine*dphi, &
                mapped_jac(2, coordinate)]
        end do
        metric = 0.0_dp
        do first = 1, 3
            do second = 1, 3
                metric(first, second) = sum( &
                    physical_jac(:, first)*physical_jac(:, second))
            end do
        end do
        volume = determinant3(physical_jac)
    end subroutine map_inverse_coordinates

    subroutine map_inverse_coordinates_tangent( &
            mapped, mapped_jac, mapped_dot, mapped_jac_dot, zeta, physical_dot, &
            physical_jac_dot, metric_dot, volume_dot)
        real(dp), intent(in) :: mapped(3), mapped_jac(3, 3)
        real(dp), intent(in) :: mapped_dot(3), mapped_jac_dot(3, 3), zeta
        real(dp), intent(out) :: physical_dot(3), physical_jac_dot(3, 3)
        real(dp), intent(out) :: metric_dot(3, 3), volume_dot

        real(dp) :: radius, lambda, phi, cosine, sine, radius_dot
        real(dp) :: lambda_dot, dphi, dphi_dot
        real(dp) :: physical_jac(3, 3), metric(3, 3), physical(3)
        integer :: first, second, index

        call map_inverse_coordinates( &
            mapped, mapped_jac, zeta, physical, physical_jac, &
            metric, volume_dot)
        radius = mapped(1)
        lambda = mapped(3)
        phi = zeta + lambda
        cosine = cos(phi)
        sine = sin(phi)
        radius_dot = mapped_dot(1)
        lambda_dot = mapped_dot(3)
        physical_dot = [ &
            cosine*radius_dot - radius*sine*lambda_dot, &
            sine*radius_dot + radius*cosine*lambda_dot, mapped_dot(2)]
        do index = 1, 3
            dphi = mapped_jac(3, index)
            if (index == 3) dphi = dphi + 1.0_dp
            dphi_dot = mapped_jac_dot(3, index)
            physical_jac_dot(1, index) = &
                cosine*mapped_jac_dot(1, index) - &
                sine*lambda_dot*mapped_jac(1, index) - &
                radius_dot*sine*dphi - radius*cosine*lambda_dot*dphi - &
                radius*sine*dphi_dot
            physical_jac_dot(2, index) = &
                sine*mapped_jac_dot(1, index) + &
                cosine*lambda_dot*mapped_jac(1, index) + &
                radius_dot*cosine*dphi - radius*sine*lambda_dot*dphi + &
                radius*cosine*dphi_dot
            physical_jac_dot(3, index) = mapped_jac_dot(2, index)
        end do
        metric_dot = 0.0_dp
        do first = 1, 3
            do second = 1, 3
                metric_dot(first, second) = sum( &
                    physical_jac_dot(:, first)*physical_jac(:, second) + &
                    physical_jac(:, first)*physical_jac_dot(:, second))
            end do
        end do
        volume_dot = determinant3_directional(physical_jac, physical_jac_dot)
    end subroutine map_inverse_coordinates_tangent

    subroutine reverse_inverse_map( &
            mapped, mapped_jac, zeta, physical_bar, physical_jac_bar, metric_bar, &
            volume_bar, mapped_bar, mapped_jac_bar)
        real(dp), intent(in) :: mapped(3), mapped_jac(3, 3), zeta
        real(dp), intent(in) :: physical_bar(3), physical_jac_bar(3, 3)
        real(dp), intent(in) :: metric_bar(3, 3), volume_bar
        real(dp), intent(inout) :: mapped_bar(3), mapped_jac_bar(3, 3)

        real(dp) :: radius, lambda, phi, cosine, sine, dphi
        real(dp) :: physical_jac(3, 3), metric(3, 3), physical(3)
        real(dp) :: cofactor(3, 3), volume
        real(dp) :: physical_jac_bar_total(3, 3)
        integer :: coordinate, first, second, component

        call map_inverse_coordinates( &
            mapped, mapped_jac, zeta, physical, physical_jac, &
            metric, volume)
        radius = mapped(1)
        lambda = mapped(3)
        phi = zeta + lambda
        cosine = cos(phi)
        sine = sin(phi)
        physical_jac_bar_total = physical_jac_bar
        do first = 1, 3
            do second = 1, 3
                do component = 1, 3
                    physical_jac_bar_total(component, first) = &
                        physical_jac_bar_total(component, first) + &
                        physical_jac(component, second)* &
                        (metric_bar(first, second) + metric_bar(second, first))
                end do
            end do
        end do
        call determinant3_gradient(physical_jac, cofactor)
        physical_jac_bar_total = physical_jac_bar_total + volume_bar*cofactor
        mapped_bar(1) = mapped_bar(1) + cosine*physical_bar(1) + &
            sine*physical_bar(2)
        mapped_bar(2) = mapped_bar(2) + physical_bar(3)
        mapped_bar(3) = mapped_bar(3) - radius*sine*physical_bar(1) + &
            radius*cosine*physical_bar(2)
        do coordinate = 1, 3
            dphi = mapped_jac(3, coordinate)
            if (coordinate == 3) dphi = dphi + 1.0_dp
            mapped_jac_bar(1, coordinate) = mapped_jac_bar(1, coordinate) + &
                cosine*physical_jac_bar_total(1, coordinate) + &
                sine*physical_jac_bar_total(2, coordinate)
            mapped_jac_bar(2, coordinate) = mapped_jac_bar(2, coordinate) + &
                physical_jac_bar_total(3, coordinate)
            mapped_jac_bar(3, coordinate) = mapped_jac_bar(3, coordinate) - &
                radius*sine*physical_jac_bar_total(1, coordinate) + &
                radius*cosine*physical_jac_bar_total(2, coordinate)
            mapped_bar(1) = mapped_bar(1) - sine*dphi* &
                physical_jac_bar_total(1, coordinate) + &
                cosine*dphi*physical_jac_bar_total(2, coordinate)
            mapped_bar(3) = mapped_bar(3) + &
                (-sine*mapped_jac(1, coordinate) - &
                radius*cosine*dphi)*physical_jac_bar_total(1, coordinate) + &
                (cosine*mapped_jac(1, coordinate) - &
                radius*sine*dphi)*physical_jac_bar_total(2, coordinate)
        end do
    end subroutine reverse_inverse_map

    pure real(dp) function determinant3(matrix) result(value)
        real(dp), intent(in) :: matrix(3, 3)

        value = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2)) - &
            matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))
    end function determinant3

    pure real(dp) function determinant3_directional( &
            matrix, matrix_dot) result(value)
        real(dp), intent(in) :: matrix(3, 3), matrix_dot(3, 3)
        real(dp) :: gradient(3, 3)

        call determinant3_gradient(matrix, gradient)
        value = sum(gradient*matrix_dot)
    end function determinant3_directional

    pure subroutine determinant3_gradient(matrix, gradient)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: gradient(3, 3)

        gradient(1, 1) = matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)
        gradient(1, 2) = matrix(2, 3)*matrix(3, 1) - matrix(2, 1)*matrix(3, 3)
        gradient(1, 3) = matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1)
        gradient(2, 1) = matrix(1, 3)*matrix(3, 2) - matrix(1, 2)*matrix(3, 3)
        gradient(2, 2) = matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1)
        gradient(2, 3) = matrix(1, 2)*matrix(3, 1) - matrix(1, 1)*matrix(3, 2)
        gradient(3, 1) = matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2)
        gradient(3, 2) = matrix(1, 3)*matrix(2, 1) - matrix(1, 1)*matrix(2, 3)
        gradient(3, 3) = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
    end subroutine determinant3_gradient

    logical function finite_real_1d(values) result(valid)
        real(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(values))
    end function finite_real_1d

    logical function finite_real_2d(values) result(valid)
        real(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(values))
    end function finite_real_2d

    logical function finite_real_3d(values) result(valid)
        real(dp), intent(in) :: values(:, :, :)

        valid = all(ieee_is_finite(values))
    end function finite_real_3d

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:, :)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

end module fortfem_nested_surface_geometry
