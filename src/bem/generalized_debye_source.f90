module fortfem_generalized_debye_source
    !! Neutral generalized-Debye-source coordinate and residual layer.
    !!
    !! A tangential surface field is supplied in caller-owned coordinates as
    !!
    !!   j = G a + C b + H h,
    !!
    !! where G and C are two scalar-to-tangential lifts and H contains the
    !! harmonic representatives of the closed surface.  S and P are the
    !! caller's generalized-source and period operators.  This module forms
    !! the resulting current and residuals; it does not choose a Green kernel,
    !! Beltrami parameter, topology, or physical source normalization.
    !!
    !! The contract is the coordinate layer needed before a BIEST-like
    !! second-kind operator.  It remains useful for RWG/BC, high-order traces,
    !! IGA, and nonmatching surface spaces because all maps are explicit.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: assemble_generalized_debye_source_residual
    public :: assemble_generalized_debye_source_residual_jvp
    public :: assemble_generalized_debye_source_residual_vjp
    public :: assemble_generalized_debye_source_second_kind
    public :: assemble_generalized_debye_source_second_kind_jvp
    public :: assemble_generalized_debye_source_second_kind_vjp

contains

    subroutine assemble_generalized_debye_source_residual( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, surface_current, &
            source_residual, period_residual, status)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(out) :: surface_current(:), source_residual(:)
        complex(dp), intent(out) :: period_residual(:)
        integer, intent(out) :: status

        status = 1
        surface_current = cmplx(0.0_dp, 0.0_dp, dp)
        source_residual = cmplx(0.0_dp, 0.0_dp, dp)
        period_residual = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_inputs( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, surface_current, &
            source_residual, period_residual)) return

        surface_current = matmul(gradient_lift, gradient_coefficients) + &
            matmul(cogradient_lift, cogradient_coefficients) + &
            matmul(harmonic_basis, harmonic_coefficients)
        source_residual = matmul(source_operator, surface_current) - source_target
        period_residual = matmul(period_operator, surface_current) - period_target
        status = 0
    end subroutine assemble_generalized_debye_source_residual

    subroutine assemble_generalized_debye_source_residual_jvp( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, gradient_lift_dot, &
            cogradient_lift_dot, harmonic_basis_dot, source_operator_dot, &
            period_operator_dot, gradient_coefficients_dot, cogradient_coefficients_dot, &
            harmonic_coefficients_dot, source_target_dot, period_target_dot, &
            surface_current_dot, source_residual_dot, period_residual_dot, status)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(in) :: gradient_lift_dot(:, :), cogradient_lift_dot(:, :)
        complex(dp), intent(in) :: harmonic_basis_dot(:, :), source_operator_dot(:, :)
        complex(dp), intent(in) :: period_operator_dot(:, :)
        complex(dp), intent(in) :: gradient_coefficients_dot(:)
        complex(dp), intent(in) :: cogradient_coefficients_dot(:)
        complex(dp), intent(in) :: harmonic_coefficients_dot(:)
        complex(dp), intent(in) :: source_target_dot(:), period_target_dot(:)
        complex(dp), intent(out) :: surface_current_dot(:), source_residual_dot(:)
        complex(dp), intent(out) :: period_residual_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: current(:)

        status = 1
        surface_current_dot = cmplx(0.0_dp, 0.0_dp, dp)
        source_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        period_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_inputs( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, surface_current_dot, &
            source_residual_dot, period_residual_dot)) return
        if (.not. valid_direction( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, gradient_lift_dot, &
            cogradient_lift_dot, harmonic_basis_dot, source_operator_dot, &
            period_operator_dot, gradient_coefficients_dot, cogradient_coefficients_dot, &
            harmonic_coefficients_dot, source_target_dot, period_target_dot)) return

        allocate(current(size(surface_current_dot)))
        current = matmul(gradient_lift, gradient_coefficients) + &
            matmul(cogradient_lift, cogradient_coefficients) + &
            matmul(harmonic_basis, harmonic_coefficients)
        surface_current_dot = matmul(gradient_lift_dot, gradient_coefficients) + &
            matmul(gradient_lift, gradient_coefficients_dot) + &
            matmul(cogradient_lift_dot, cogradient_coefficients) + &
            matmul(cogradient_lift, cogradient_coefficients_dot) + &
            matmul(harmonic_basis_dot, harmonic_coefficients) + &
            matmul(harmonic_basis, harmonic_coefficients_dot)
        source_residual_dot = matmul(source_operator_dot, current) + &
            matmul(source_operator, surface_current_dot) - source_target_dot
        period_residual_dot = matmul(period_operator_dot, current) + &
            matmul(period_operator, surface_current_dot) - period_target_dot
        status = 0
    end subroutine assemble_generalized_debye_source_residual_jvp

    subroutine assemble_generalized_debye_source_residual_vjp( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, source_residual_bar, &
            period_residual_bar, surface_current_bar, gradient_lift_bar, &
            cogradient_lift_bar, harmonic_basis_bar, source_operator_bar, &
            period_operator_bar, gradient_coefficients_bar, cogradient_coefficients_bar, &
            harmonic_coefficients_bar, source_target_bar, period_target_bar, status)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(in) :: source_residual_bar(:), period_residual_bar(:)
        complex(dp), intent(in) :: surface_current_bar(:)
        complex(dp), intent(out) :: gradient_lift_bar(:, :), cogradient_lift_bar(:, :)
        complex(dp), intent(out) :: harmonic_basis_bar(:, :), source_operator_bar(:, :)
        complex(dp), intent(out) :: period_operator_bar(:, :)
        complex(dp), intent(out) :: gradient_coefficients_bar(:)
        complex(dp), intent(out) :: cogradient_coefficients_bar(:)
        complex(dp), intent(out) :: harmonic_coefficients_bar(:)
        complex(dp), intent(out) :: source_target_bar(:), period_target_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: current(:), current_bar(:)

        gradient_lift_bar = cmplx(0.0_dp, 0.0_dp, dp)
        cogradient_lift_bar = cmplx(0.0_dp, 0.0_dp, dp)
        harmonic_basis_bar = cmplx(0.0_dp, 0.0_dp, dp)
        source_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        period_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        gradient_coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        cogradient_coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        harmonic_coefficients_bar = cmplx(0.0_dp, 0.0_dp, dp)
        source_target_bar = cmplx(0.0_dp, 0.0_dp, dp)
        period_target_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        allocate(current(size(gradient_lift, 1)))
        current = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. valid_inputs( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, current, &
            source_residual_bar, period_residual_bar)) return
        if (.not. valid_adjoint_inputs( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, surface_current_bar, gradient_lift_bar, &
            cogradient_lift_bar, harmonic_basis_bar, source_operator_bar, &
            period_operator_bar, gradient_coefficients_bar, &
            cogradient_coefficients_bar, harmonic_coefficients_bar, &
            source_target_bar, period_target_bar)) return

        allocate(current_bar(size(surface_current_bar)))
        current = matmul(gradient_lift, gradient_coefficients) + &
            matmul(cogradient_lift, cogradient_coefficients) + &
            matmul(harmonic_basis, harmonic_coefficients)
        current_bar = surface_current_bar + &
            matmul(conjg(transpose(source_operator)), source_residual_bar) + &
            matmul(conjg(transpose(period_operator)), period_residual_bar)
        source_operator_bar = matmul(reshape(source_residual_bar, &
            [size(source_residual_bar), 1]), reshape(conjg(current), &
            [1, size(current)]))
        period_operator_bar = matmul(reshape(period_residual_bar, &
            [size(period_residual_bar), 1]), reshape(conjg(current), &
            [1, size(current)]))
        source_target_bar = -source_residual_bar
        period_target_bar = -period_residual_bar
        gradient_lift_bar = matmul(reshape(current_bar, [size(current_bar), 1]), &
            reshape(conjg(gradient_coefficients), [1, size(gradient_coefficients)]))
        cogradient_lift_bar = matmul(reshape(current_bar, [size(current_bar), 1]), &
            reshape(conjg(cogradient_coefficients), [1, size(cogradient_coefficients)]))
        harmonic_basis_bar = matmul(reshape(current_bar, [size(current_bar), 1]), &
            reshape(conjg(harmonic_coefficients), [1, size(harmonic_coefficients)]))
        gradient_coefficients_bar = matmul(conjg(transpose(gradient_lift)), current_bar)
        cogradient_coefficients_bar = matmul( &
            conjg(transpose(cogradient_lift)), current_bar)
        harmonic_coefficients_bar = matmul(conjg(transpose(harmonic_basis)), current_bar)
        status = 0
    end subroutine assemble_generalized_debye_source_residual_vjp

    subroutine assemble_generalized_debye_source_second_kind( &
            second_kind_operator, second_kind_target, gradient_lift, &
            cogradient_lift, harmonic_basis, source_operator, period_operator, &
            gradient_coefficients, cogradient_coefficients, harmonic_coefficients, &
            source_target, period_target, surface_current, second_kind_residual, &
            source_residual, period_residual, status)
        complex(dp), intent(in) :: second_kind_operator(:, :), second_kind_target(:)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(out) :: surface_current(:), second_kind_residual(:)
        complex(dp), intent(out) :: source_residual(:), period_residual(:)
        integer, intent(out) :: status

        call assemble_generalized_debye_source_residual( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, surface_current, &
            source_residual, period_residual, status)
        if (status /= 0) then
            second_kind_residual = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        if (.not. valid_second_kind_inputs( &
            second_kind_operator, second_kind_target, surface_current, &
            second_kind_residual)) then
            status = 1
            second_kind_residual = cmplx(0.0_dp, 0.0_dp, dp)
            return
        end if
        second_kind_residual = matmul(second_kind_operator, surface_current) - &
            second_kind_target
        status = 0
    end subroutine assemble_generalized_debye_source_second_kind

    subroutine assemble_generalized_debye_source_second_kind_jvp( &
            second_kind_operator, second_kind_target, gradient_lift, &
            cogradient_lift, harmonic_basis, source_operator, period_operator, &
            gradient_coefficients, cogradient_coefficients, harmonic_coefficients, &
            source_target, period_target, second_kind_operator_dot, &
            second_kind_target_dot, gradient_lift_dot, cogradient_lift_dot, &
            harmonic_basis_dot, source_operator_dot, period_operator_dot, &
            gradient_coefficients_dot, cogradient_coefficients_dot, &
            harmonic_coefficients_dot, source_target_dot, period_target_dot, &
            surface_current_dot, second_kind_residual_dot, source_residual_dot, &
            period_residual_dot, status)
        complex(dp), intent(in) :: second_kind_operator(:, :), second_kind_target(:)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(in) :: second_kind_operator_dot(:, :)
        complex(dp), intent(in) :: second_kind_target_dot(:)
        complex(dp), intent(in) :: gradient_lift_dot(:, :), cogradient_lift_dot(:, :)
        complex(dp), intent(in) :: harmonic_basis_dot(:, :), source_operator_dot(:, :)
        complex(dp), intent(in) :: period_operator_dot(:, :)
        complex(dp), intent(in) :: gradient_coefficients_dot(:)
        complex(dp), intent(in) :: cogradient_coefficients_dot(:)
        complex(dp), intent(in) :: harmonic_coefficients_dot(:)
        complex(dp), intent(in) :: source_target_dot(:), period_target_dot(:)
        complex(dp), intent(out) :: surface_current_dot(:), second_kind_residual_dot(:)
        complex(dp), intent(out) :: source_residual_dot(:), period_residual_dot(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: surface_current(:), source_residual(:)
        complex(dp), allocatable :: period_residual(:)

        surface_current_dot = cmplx(0.0_dp, 0.0_dp, dp)
        second_kind_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        source_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        period_residual_dot = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        allocate(surface_current(size(gradient_lift, 1)), &
            source_residual(size(source_operator, 1)), &
            period_residual(size(period_operator, 1)))
        call assemble_generalized_debye_source_second_kind( &
            second_kind_operator, second_kind_target, gradient_lift, cogradient_lift, &
            harmonic_basis, source_operator, period_operator, gradient_coefficients, &
            cogradient_coefficients, harmonic_coefficients, source_target, period_target, &
            surface_current, second_kind_residual_dot, source_residual, period_residual, &
            status)
        if (status /= 0) return
        if (.not. same_shape(second_kind_operator, second_kind_operator_dot) .or. &
            size(second_kind_target_dot) /= size(second_kind_target) .or. &
            .not. finite_complex(second_kind_operator_dot) .or. &
            .not. finite_complex(second_kind_target_dot)) then
            status = 1
            return
        end if
        call assemble_generalized_debye_source_residual_jvp( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, gradient_lift_dot, &
            cogradient_lift_dot, harmonic_basis_dot, source_operator_dot, &
            period_operator_dot, gradient_coefficients_dot, cogradient_coefficients_dot, &
            harmonic_coefficients_dot, source_target_dot, period_target_dot, &
            surface_current_dot, source_residual_dot, period_residual_dot, status)
        if (status /= 0) return
        second_kind_residual_dot = matmul(second_kind_operator_dot, surface_current) + &
            matmul(second_kind_operator, surface_current_dot) - second_kind_target_dot
        status = 0
    end subroutine assemble_generalized_debye_source_second_kind_jvp

    subroutine assemble_generalized_debye_source_second_kind_vjp( &
            second_kind_operator, second_kind_target, gradient_lift, &
            cogradient_lift, harmonic_basis, source_operator, period_operator, &
            gradient_coefficients, cogradient_coefficients, harmonic_coefficients, &
            source_target, period_target, second_kind_residual_bar, source_residual_bar, &
            period_residual_bar, surface_current_bar, second_kind_operator_bar, &
            second_kind_target_bar, gradient_lift_bar, cogradient_lift_bar, &
            harmonic_basis_bar, source_operator_bar, period_operator_bar, &
            gradient_coefficients_bar, cogradient_coefficients_bar, &
            harmonic_coefficients_bar, source_target_bar, period_target_bar, status)
        complex(dp), intent(in) :: second_kind_operator(:, :), second_kind_target(:)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(in) :: second_kind_residual_bar(:)
        complex(dp), intent(in) :: source_residual_bar(:), period_residual_bar(:)
        complex(dp), intent(in) :: surface_current_bar(:)
        complex(dp), intent(out) :: second_kind_operator_bar(:, :)
        complex(dp), intent(out) :: second_kind_target_bar(:)
        complex(dp), intent(out) :: gradient_lift_bar(:, :), cogradient_lift_bar(:, :)
        complex(dp), intent(out) :: harmonic_basis_bar(:, :), source_operator_bar(:, :)
        complex(dp), intent(out) :: period_operator_bar(:, :)
        complex(dp), intent(out) :: gradient_coefficients_bar(:)
        complex(dp), intent(out) :: cogradient_coefficients_bar(:)
        complex(dp), intent(out) :: harmonic_coefficients_bar(:)
        complex(dp), intent(out) :: source_target_bar(:), period_target_bar(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: surface_current(:), second_kind_residual(:)
        complex(dp), allocatable :: source_residual(:), period_residual(:)
        complex(dp), allocatable :: current_bar(:)

        second_kind_operator_bar = cmplx(0.0_dp, 0.0_dp, dp)
        second_kind_target_bar = cmplx(0.0_dp, 0.0_dp, dp)
        status = 1
        allocate(surface_current(size(gradient_lift, 1)), &
            second_kind_residual(size(second_kind_operator, 1)), &
            source_residual(size(source_operator, 1)), &
            period_residual(size(period_operator, 1)), &
            current_bar(size(gradient_lift, 1)))
        call assemble_generalized_debye_source_second_kind( &
            second_kind_operator, second_kind_target, gradient_lift, cogradient_lift, &
            harmonic_basis, source_operator, period_operator, gradient_coefficients, &
            cogradient_coefficients, harmonic_coefficients, source_target, period_target, &
            surface_current, second_kind_residual, source_residual, period_residual, &
            status)
        if (status /= 0 .or. size(second_kind_residual_bar) /= size(second_kind_target) .or. &
            size(surface_current_bar) /= size(surface_current) .or. &
            .not. finite_complex(second_kind_residual_bar) .or. &
            .not. finite_complex(surface_current_bar) .or. &
            .not. finite_complex(source_residual_bar) .or. &
            .not. finite_complex(period_residual_bar)) then
            status = 1
            return
        end if
        second_kind_operator_bar = matmul(reshape(second_kind_residual_bar, &
            [size(second_kind_residual_bar), 1]), reshape(conjg(surface_current), &
            [1, size(surface_current)]))
        second_kind_target_bar = -second_kind_residual_bar
        current_bar = surface_current_bar + matmul( &
            conjg(transpose(second_kind_operator)), second_kind_residual_bar)
        call assemble_generalized_debye_source_residual_vjp( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, source_residual_bar, &
            period_residual_bar, current_bar, gradient_lift_bar, cogradient_lift_bar, &
            harmonic_basis_bar, source_operator_bar, period_operator_bar, &
            gradient_coefficients_bar, cogradient_coefficients_bar, &
            harmonic_coefficients_bar, source_target_bar, period_target_bar, status)
    end subroutine assemble_generalized_debye_source_second_kind_vjp

    logical function valid_second_kind_inputs( &
            second_kind_operator, second_kind_target, surface_current, &
            second_kind_residual) result(valid)
        complex(dp), intent(in) :: second_kind_operator(:, :), second_kind_target(:)
        complex(dp), intent(in) :: surface_current(:), second_kind_residual(:)

        valid = size(second_kind_operator, 1) > 0 .and. &
            size(second_kind_operator, 1) == size(second_kind_operator, 2) .and. &
            size(second_kind_target) == size(second_kind_operator, 1) .and. &
            size(surface_current) == size(second_kind_operator, 2) .and. &
            size(second_kind_residual) == size(second_kind_operator, 1) .and. &
            finite_complex(second_kind_operator) .and. &
            finite_complex(second_kind_target) .and. &
            finite_complex(surface_current) .and. &
            finite_complex(second_kind_residual)
    end function valid_second_kind_inputs

    logical function valid_inputs( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, surface_current, &
            source_residual, period_residual) result(valid)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(in) :: surface_current(:), source_residual(:)
        complex(dp), intent(in) :: period_residual(:)
        integer :: tangent_count, source_count, period_count

        tangent_count = size(gradient_lift, 1)
        source_count = size(source_operator, 1)
        period_count = size(period_operator, 1)
        valid = tangent_count > 0 .and. size(gradient_lift, 2) > 0 .and. &
            size(cogradient_lift, 1) == tangent_count .and. &
            size(cogradient_lift, 2) == size(gradient_coefficients) .and. &
            size(gradient_lift, 2) == size(gradient_coefficients) .and. &
            size(harmonic_basis, 1) == tangent_count .and. &
            size(harmonic_basis, 2) == size(harmonic_coefficients) .and. &
            size(source_operator, 2) == tangent_count .and. &
            size(period_operator, 2) == tangent_count .and. &
            size(source_target) == source_count .and. &
            size(period_target) == period_count .and. &
            size(surface_current) == tangent_count .and. &
            size(source_residual) == source_count .and. &
            size(period_residual) == period_count
        if (.not. valid) return
        valid = finite_complex(gradient_lift) .and. finite_complex(cogradient_lift) .and. &
            finite_complex(harmonic_basis) .and. finite_complex(source_operator) .and. &
            finite_complex(period_operator) .and. finite_complex(gradient_coefficients) .and. &
            finite_complex(cogradient_coefficients) .and. &
            finite_complex(harmonic_coefficients) .and. finite_complex(source_target) .and. &
            finite_complex(period_target) .and. finite_complex(surface_current) .and. &
            finite_complex(source_residual) .and. finite_complex(period_residual)
    end function valid_inputs

    logical function valid_direction( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, gradient_coefficients, cogradient_coefficients, &
            harmonic_coefficients, source_target, period_target, gradient_lift_dot, &
            cogradient_lift_dot, harmonic_basis_dot, source_operator_dot, &
            period_operator_dot, gradient_coefficients_dot, cogradient_coefficients_dot, &
            harmonic_coefficients_dot, source_target_dot, period_target_dot) result(valid)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :)
        complex(dp), intent(in) :: gradient_coefficients(:)
        complex(dp), intent(in) :: cogradient_coefficients(:), harmonic_coefficients(:)
        complex(dp), intent(in) :: source_target(:), period_target(:)
        complex(dp), intent(in) :: gradient_lift_dot(:, :), cogradient_lift_dot(:, :)
        complex(dp), intent(in) :: harmonic_basis_dot(:, :), source_operator_dot(:, :)
        complex(dp), intent(in) :: period_operator_dot(:, :)
        complex(dp), intent(in) :: gradient_coefficients_dot(:)
        complex(dp), intent(in) :: cogradient_coefficients_dot(:)
        complex(dp), intent(in) :: harmonic_coefficients_dot(:)
        complex(dp), intent(in) :: source_target_dot(:), period_target_dot(:)

        valid = same_shape(gradient_lift, gradient_lift_dot) .and. &
            same_shape(cogradient_lift, cogradient_lift_dot) .and. &
            same_shape(harmonic_basis, harmonic_basis_dot) .and. &
            same_shape(source_operator, source_operator_dot) .and. &
            same_shape(period_operator, period_operator_dot) .and. &
            size(gradient_coefficients_dot) == size(gradient_coefficients) .and. &
            size(cogradient_coefficients_dot) == size(cogradient_coefficients) .and. &
            size(harmonic_coefficients_dot) == size(harmonic_coefficients) .and. &
            size(source_target_dot) == size(source_target) .and. &
            size(period_target_dot) == size(period_target) .and. &
            finite_complex(gradient_lift_dot) .and. finite_complex(cogradient_lift_dot) .and. &
            finite_complex(harmonic_basis_dot) .and. finite_complex(source_operator_dot) .and. &
            finite_complex(period_operator_dot) .and. finite_complex(gradient_coefficients_dot) .and. &
            finite_complex(cogradient_coefficients_dot) .and. &
            finite_complex(harmonic_coefficients_dot) .and. finite_complex(source_target_dot) .and. &
            finite_complex(period_target_dot)
    end function valid_direction

    logical function valid_adjoint_inputs( &
            gradient_lift, cogradient_lift, harmonic_basis, source_operator, &
            period_operator, surface_current_bar, gradient_lift_bar, &
            cogradient_lift_bar, harmonic_basis_bar, source_operator_bar, &
            period_operator_bar, gradient_coefficients_bar, cogradient_coefficients_bar, &
            harmonic_coefficients_bar, source_target_bar, period_target_bar) result(valid)
        complex(dp), intent(in) :: gradient_lift(:, :), cogradient_lift(:, :)
        complex(dp), intent(in) :: harmonic_basis(:, :), source_operator(:, :)
        complex(dp), intent(in) :: period_operator(:, :), surface_current_bar(:)
        complex(dp), intent(in) :: gradient_lift_bar(:, :), cogradient_lift_bar(:, :)
        complex(dp), intent(in) :: harmonic_basis_bar(:, :), source_operator_bar(:, :)
        complex(dp), intent(in) :: period_operator_bar(:, :)
        complex(dp), intent(in) :: gradient_coefficients_bar(:)
        complex(dp), intent(in) :: cogradient_coefficients_bar(:)
        complex(dp), intent(in) :: harmonic_coefficients_bar(:)
        complex(dp), intent(in) :: source_target_bar(:), period_target_bar(:)

        valid = size(surface_current_bar) == size(gradient_lift, 1) .and. &
            same_shape(gradient_lift, gradient_lift_bar) .and. &
            same_shape(cogradient_lift, cogradient_lift_bar) .and. &
            same_shape(harmonic_basis, harmonic_basis_bar) .and. &
            same_shape(source_operator, source_operator_bar) .and. &
            same_shape(period_operator, period_operator_bar) .and. &
            size(gradient_coefficients_bar) == size(gradient_lift, 2) .and. &
            size(cogradient_coefficients_bar) == size(cogradient_lift, 2) .and. &
            size(harmonic_coefficients_bar) == size(harmonic_basis, 2) .and. &
            size(source_target_bar) == size(source_operator, 1) .and. &
            size(period_target_bar) == size(period_operator, 1) .and. &
            finite_complex(surface_current_bar)
    end function valid_adjoint_inputs

    pure logical function same_shape(first, second) result(equal)
        complex(dp), intent(in) :: first(:, :), second(:, :)

        equal = all(shape(first) == shape(second))
    end function same_shape

    pure logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(..)

        select rank (values)
            rank (1)
            valid = all(ieee_is_finite(real(values, dp))) .and. &
                all(ieee_is_finite(aimag(values)))
            rank (2)
            valid = all(ieee_is_finite(real(values, dp))) .and. &
                all(ieee_is_finite(aimag(values)))
            rank default
            valid = .false.
        end select
    end function finite_complex

end module fortfem_generalized_debye_source
