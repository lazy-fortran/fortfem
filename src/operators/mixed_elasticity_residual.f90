module fortfem_mixed_elasticity_residual
    !! Neutral first-order mixed elasticity residual and derivative actions.
    !!
    !! A caller supplies the compliance-like constitutive map C, the
    !! displacement-to-strain map E, and the stress-to-equilibrium map D.  The
    !! residual is
    !!
    !!     r_c = C sigma - E u,       r_e = D sigma - f.
    !!
    !! This is the algebraic contract behind Hellinger--Reissner, weak-
    !! symmetry, TDNNS, and anisotropic mixed elasticity clients.  FortFEM
    !! does not select a constitutive law, element family, or boundary model.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_mixed_elasticity_residual
    public :: assemble_mixed_elasticity_residual_jvp
    public :: assemble_mixed_elasticity_residual_vjp

contains

    subroutine assemble_mixed_elasticity_residual( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            constitutive_residual, equilibrium_residual, status)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:)
        real(dp), intent(out) :: constitutive_residual(:), equilibrium_residual(:)
        integer, intent(out) :: status

        constitutive_residual = 0.0_dp
        equilibrium_residual = 0.0_dp
        status = 1
        if (.not. valid_value_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            constitutive_residual, equilibrium_residual) .or. &
            .not. finite_values(compliance, strain_map, divergence_map, stress, &
            displacement, load)) return

        constitutive_residual = matmul(compliance, stress) - &
            matmul(strain_map, displacement)
        equilibrium_residual = matmul(divergence_map, stress) - load
        status = 0
    end subroutine assemble_mixed_elasticity_residual

    subroutine assemble_mixed_elasticity_residual_jvp( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            compliance_dot, strain_map_dot, divergence_map_dot, stress_dot, &
            displacement_dot, load_dot, constitutive_residual_dot, &
            equilibrium_residual_dot, status)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:)
        real(dp), intent(in) :: compliance_dot(:, :), strain_map_dot(:, :)
        real(dp), intent(in) :: divergence_map_dot(:, :), stress_dot(:)
        real(dp), intent(in) :: displacement_dot(:), load_dot(:)
        real(dp), intent(out) :: constitutive_residual_dot(:)
        real(dp), intent(out) :: equilibrium_residual_dot(:)
        integer, intent(out) :: status

        constitutive_residual_dot = 0.0_dp
        equilibrium_residual_dot = 0.0_dp
        status = 1
        if (.not. valid_value_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            constitutive_residual_dot, equilibrium_residual_dot) .or. &
            .not. valid_direction_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            compliance_dot, strain_map_dot, divergence_map_dot, stress_dot, &
            displacement_dot, load_dot) .or. &
            .not. finite_values(compliance, strain_map, divergence_map, stress, &
            displacement, load) .or. &
            .not. finite_direction(compliance_dot, strain_map_dot, &
            divergence_map_dot, stress_dot, displacement_dot, load_dot)) return

        constitutive_residual_dot = matmul(compliance_dot, stress) + &
            matmul(compliance, stress_dot) - matmul(strain_map_dot, displacement) - &
            matmul(strain_map, displacement_dot)
        equilibrium_residual_dot = matmul(divergence_map_dot, stress) + &
            matmul(divergence_map, stress_dot) - load_dot
        status = 0
    end subroutine assemble_mixed_elasticity_residual_jvp

    subroutine assemble_mixed_elasticity_residual_vjp( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            constitutive_residual_bar, equilibrium_residual_bar, compliance_bar, &
            strain_map_bar, divergence_map_bar, stress_bar, displacement_bar, &
            load_bar, status)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:), constitutive_residual_bar(:)
        real(dp), intent(in) :: equilibrium_residual_bar(:)
        real(dp), intent(out) :: compliance_bar(:, :), strain_map_bar(:, :)
        real(dp), intent(out) :: divergence_map_bar(:, :), stress_bar(:)
        real(dp), intent(out) :: displacement_bar(:), load_bar(:)
        integer, intent(out) :: status

        compliance_bar = 0.0_dp
        strain_map_bar = 0.0_dp
        divergence_map_bar = 0.0_dp
        stress_bar = 0.0_dp
        displacement_bar = 0.0_dp
        load_bar = 0.0_dp
        status = 1
        if (.not. valid_value_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            constitutive_residual_bar, equilibrium_residual_bar) .or. &
            .not. valid_cotangent_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            compliance_bar, strain_map_bar, divergence_map_bar, stress_bar, &
            displacement_bar, load_bar) .or. &
            .not. finite_values(compliance, strain_map, divergence_map, stress, &
            displacement, load) .or. &
            .not. all(ieee_is_finite(constitutive_residual_bar)) .or. &
            .not. all(ieee_is_finite(equilibrium_residual_bar))) return

        compliance_bar = outer_product(constitutive_residual_bar, stress)
        strain_map_bar = -outer_product(constitutive_residual_bar, displacement)
        divergence_map_bar = outer_product(equilibrium_residual_bar, stress)
        stress_bar = matmul(transpose(compliance), constitutive_residual_bar) + &
            matmul(transpose(divergence_map), equilibrium_residual_bar)
        displacement_bar = -matmul(transpose(strain_map), constitutive_residual_bar)
        load_bar = -equilibrium_residual_bar
        status = 0
    end subroutine assemble_mixed_elasticity_residual_vjp

    logical function valid_value_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            constitutive_residual, equilibrium_residual) result(valid)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:), constitutive_residual(:)
        real(dp), intent(in) :: equilibrium_residual(:)

        valid = size(compliance, 1) > 0 .and. &
            size(compliance, 1) == size(compliance, 2) .and. &
            size(strain_map, 1) == size(compliance, 1) .and. &
            size(strain_map, 2) == size(displacement) .and. &
            size(divergence_map, 2) == size(stress) .and. &
            size(divergence_map, 1) == size(load) .and. &
            size(stress) == size(compliance, 1) .and. &
            size(constitutive_residual) == size(stress) .and. &
            size(equilibrium_residual) == size(load)
    end function valid_value_shapes

    logical function valid_direction_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            compliance_dot, strain_map_dot, divergence_map_dot, stress_dot, &
            displacement_dot, load_dot) result(valid)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:), compliance_dot(:, :), strain_map_dot(:, :)
        real(dp), intent(in) :: divergence_map_dot(:, :), stress_dot(:)
        real(dp), intent(in) :: displacement_dot(:), load_dot(:)

        valid = all(shape(compliance_dot) == shape(compliance)) .and. &
            all(shape(strain_map_dot) == shape(strain_map)) .and. &
            all(shape(divergence_map_dot) == shape(divergence_map)) .and. &
            size(stress_dot) == size(stress) .and. &
            size(displacement_dot) == size(displacement) .and. &
            size(load_dot) == size(load)
    end function valid_direction_shapes

    logical function valid_cotangent_shapes( &
            compliance, strain_map, divergence_map, stress, displacement, load, &
            compliance_bar, strain_map_bar, divergence_map_bar, stress_bar, &
            displacement_bar, load_bar) result(valid)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:), compliance_bar(:, :), strain_map_bar(:, :)
        real(dp), intent(in) :: divergence_map_bar(:, :), stress_bar(:)
        real(dp), intent(in) :: displacement_bar(:), load_bar(:)

        valid = all(shape(compliance_bar) == shape(compliance)) .and. &
            all(shape(strain_map_bar) == shape(strain_map)) .and. &
            all(shape(divergence_map_bar) == shape(divergence_map)) .and. &
            size(stress_bar) == size(stress) .and. &
            size(displacement_bar) == size(displacement) .and. &
            size(load_bar) == size(load)
    end function valid_cotangent_shapes

    logical function finite_values( &
            compliance, strain_map, divergence_map, stress, displacement, load) &
            result(valid)
        real(dp), intent(in) :: compliance(:, :), strain_map(:, :)
        real(dp), intent(in) :: divergence_map(:, :), stress(:), displacement(:)
        real(dp), intent(in) :: load(:)

        valid = all(ieee_is_finite(compliance)) .and. &
            all(ieee_is_finite(strain_map)) .and. &
            all(ieee_is_finite(divergence_map)) .and. &
            all(ieee_is_finite(stress)) .and. &
            all(ieee_is_finite(displacement)) .and. all(ieee_is_finite(load))
    end function finite_values

    logical function finite_direction( &
            compliance_dot, strain_map_dot, divergence_map_dot, stress_dot, &
            displacement_dot, load_dot) result(valid)
        real(dp), intent(in) :: compliance_dot(:, :), strain_map_dot(:, :)
        real(dp), intent(in) :: divergence_map_dot(:, :), stress_dot(:)
        real(dp), intent(in) :: displacement_dot(:), load_dot(:)

        valid = all(ieee_is_finite(compliance_dot)) .and. &
            all(ieee_is_finite(strain_map_dot)) .and. &
            all(ieee_is_finite(divergence_map_dot)) .and. &
            all(ieee_is_finite(stress_dot)) .and. &
            all(ieee_is_finite(displacement_dot)) .and. all(ieee_is_finite(load_dot))
    end function finite_direction

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))
        integer :: row

        do row = 1, size(left)
            product(row, :) = left(row)*right
        end do
    end function outer_product

end module fortfem_mixed_elasticity_residual
