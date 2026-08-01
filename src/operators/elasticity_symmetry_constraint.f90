module fortfem_elasticity_symmetry_constraint
    !! Differentiable weak-symmetry constraint for mixed elasticity blocks.
    !!
    !! A caller-owned map W extracts skew or weak-symmetry coordinates from a
    !! stress vector.  The neutral residual is r_s = W sigma - target.  The
    !! map may represent an Arnold--Falk--Winther multiplier pairing, a TDNNS
    !! normal-normal trace constraint, or another fixed-topology constraint.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_elasticity_symmetry_constraint
    public :: assemble_elasticity_symmetry_constraint_jvp
    public :: assemble_elasticity_symmetry_constraint_vjp

contains

    subroutine assemble_elasticity_symmetry_constraint( &
            symmetry_map, stress, target, residual, status)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:)
        real(dp), intent(out) :: residual(:)
        integer, intent(out) :: status

        residual = 0.0_dp
        status = 1
        if (.not. valid_shapes(symmetry_map, stress, target, residual) .or. &
            .not. finite_values(symmetry_map, stress, target)) return
        residual = matmul(symmetry_map, stress) - target
        status = 0
    end subroutine assemble_elasticity_symmetry_constraint

    subroutine assemble_elasticity_symmetry_constraint_jvp( &
            symmetry_map, stress, target, symmetry_map_dot, stress_dot, target_dot, &
            residual_dot, status)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:)
        real(dp), intent(in) :: symmetry_map_dot(:, :), stress_dot(:), target_dot(:)
        real(dp), intent(out) :: residual_dot(:)
        integer, intent(out) :: status

        residual_dot = 0.0_dp
        status = 1
        if (.not. valid_shapes(symmetry_map, stress, target, residual_dot) .or. &
            .not. valid_direction_shapes(symmetry_map, stress, target, &
            symmetry_map_dot, stress_dot, target_dot) .or. &
            .not. finite_values(symmetry_map, stress, target) .or. &
            .not. finite_direction(symmetry_map_dot, stress_dot, target_dot)) return
        residual_dot = matmul(symmetry_map_dot, stress) + &
            matmul(symmetry_map, stress_dot) - target_dot
        status = 0
    end subroutine assemble_elasticity_symmetry_constraint_jvp

    subroutine assemble_elasticity_symmetry_constraint_vjp( &
            symmetry_map, stress, target, residual_bar, symmetry_map_bar, stress_bar, &
            target_bar, status)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:)
        real(dp), intent(in) :: residual_bar(:)
        real(dp), intent(out) :: symmetry_map_bar(:, :), stress_bar(:), target_bar(:)
        integer, intent(out) :: status

        symmetry_map_bar = 0.0_dp
        stress_bar = 0.0_dp
        target_bar = 0.0_dp
        status = 1
        if (.not. valid_shapes(symmetry_map, stress, target, residual_bar) .or. &
            .not. valid_cotangent_shapes(symmetry_map, stress, target, &
            symmetry_map_bar, stress_bar, target_bar) .or. &
            .not. finite_values(symmetry_map, stress, target) .or. &
            .not. all(ieee_is_finite(residual_bar))) return
        symmetry_map_bar = outer_product(residual_bar, stress)
        stress_bar = matmul(transpose(symmetry_map), residual_bar)
        target_bar = -residual_bar
        status = 0
    end subroutine assemble_elasticity_symmetry_constraint_vjp

    logical function valid_shapes(symmetry_map, stress, target, residual) result(valid)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:), residual(:)

        valid = size(symmetry_map, 1) > 0 .and. size(symmetry_map, 2) > 0 .and. &
            size(stress) == size(symmetry_map, 2) .and. &
            size(target) == size(symmetry_map, 1) .and. &
            size(residual) == size(target)
    end function valid_shapes

    logical function valid_direction_shapes( &
            symmetry_map, stress, target, symmetry_map_dot, stress_dot, target_dot) &
            result(valid)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:)
        real(dp), intent(in) :: symmetry_map_dot(:, :), stress_dot(:), target_dot(:)

        valid = all(shape(symmetry_map_dot) == shape(symmetry_map)) .and. &
            size(stress_dot) == size(stress) .and. size(target_dot) == size(target)
    end function valid_direction_shapes

    logical function valid_cotangent_shapes( &
            symmetry_map, stress, target, symmetry_map_bar, stress_bar, target_bar) &
            result(valid)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:)
        real(dp), intent(in) :: symmetry_map_bar(:, :), stress_bar(:), target_bar(:)

        valid = all(shape(symmetry_map_bar) == shape(symmetry_map)) .and. &
            size(stress_bar) == size(stress) .and. size(target_bar) == size(target)
    end function valid_cotangent_shapes

    logical function finite_values(symmetry_map, stress, target) result(valid)
        real(dp), intent(in) :: symmetry_map(:, :), stress(:), target(:)

        valid = all(ieee_is_finite(symmetry_map)) .and. &
            all(ieee_is_finite(stress)) .and. all(ieee_is_finite(target))
    end function finite_values

    logical function finite_direction(symmetry_map_dot, stress_dot, target_dot) &
            result(valid)
        real(dp), intent(in) :: symmetry_map_dot(:, :), stress_dot(:), target_dot(:)

        valid = all(ieee_is_finite(symmetry_map_dot)) .and. &
            all(ieee_is_finite(stress_dot)) .and. all(ieee_is_finite(target_dot))
    end function finite_direction

    pure function outer_product(left, right) result(product)
        real(dp), intent(in) :: left(:), right(:)
        real(dp) :: product(size(left), size(right))
        integer :: row

        do row = 1, size(left)
            product(row, :) = left(row)*right
        end do
    end function outer_product

end module fortfem_elasticity_symmetry_constraint
