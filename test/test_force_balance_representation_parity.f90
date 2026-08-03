program test_force_balance_representation_parity
    use check, only: check_condition, check_summary
    use fortfem_force_balance_representation_parity, only: &
        evaluate_force_balance_representation_parity, &
        evaluate_force_balance_representation_parity_jvp, &
        evaluate_force_balance_representation_parity_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: volume_count = 2, boundary_count = 1, sheet_count = 1
    integer, parameter :: test_count = 2, component_count = 2
    real(dp), parameter :: step = 1.0e-7_dp
    real(dp) :: av(volume_count, test_count, component_count)
    real(dp) :: af(volume_count, component_count), aw(volume_count)
    real(dp) :: ab(boundary_count, test_count, component_count)
    real(dp) :: abf(boundary_count, component_count), abw(boundary_count)
    real(dp) :: as(sheet_count, test_count, component_count)
    real(dp) :: asf(sheet_count, component_count), asw(sheet_count)
    real(dp) :: bv(volume_count, test_count, component_count)
    real(dp) :: bf(volume_count, component_count), bw(volume_count)
    real(dp) :: bb(boundary_count, test_count, component_count)
    real(dp) :: bbf(boundary_count, component_count), bbw(boundary_count)
    real(dp) :: bs(sheet_count, test_count, component_count)
    real(dp) :: bsf(sheet_count, component_count), bsw(sheet_count)
    real(dp) :: av_dot(volume_count, test_count, component_count)
    real(dp) :: af_dot(volume_count, component_count), aw_dot(volume_count)
    real(dp) :: ab_dot(boundary_count, test_count, component_count)
    real(dp) :: abf_dot(boundary_count, component_count), abw_dot(boundary_count)
    real(dp) :: as_dot(sheet_count, test_count, component_count)
    real(dp) :: asf_dot(sheet_count, component_count), asw_dot(sheet_count)
    real(dp) :: bv_dot(volume_count, test_count, component_count)
    real(dp) :: bf_dot(volume_count, component_count), bw_dot(volume_count)
    real(dp) :: bb_dot(boundary_count, test_count, component_count)
    real(dp) :: bbf_dot(boundary_count, component_count), bbw_dot(boundary_count)
    real(dp) :: bs_dot(sheet_count, test_count, component_count)
    real(dp) :: bsf_dot(sheet_count, component_count), bsw_dot(sheet_count)
    real(dp) :: difference(test_count, component_count)
    real(dp) :: difference_dot(test_count, component_count)
    real(dp) :: difference_plus(test_count, component_count)
    real(dp) :: difference_minus(test_count, component_count)
    real(dp) :: difference_bar(test_count, component_count)
    real(dp) :: expected(test_count, component_count), lhs, rhs
    real(dp) :: av_bar(volume_count, test_count, component_count)
    real(dp) :: af_bar(volume_count, component_count), aw_bar(volume_count)
    real(dp) :: ab_bar(boundary_count, test_count, component_count)
    real(dp) :: abf_bar(boundary_count, component_count), abw_bar(boundary_count)
    real(dp) :: as_bar(sheet_count, test_count, component_count)
    real(dp) :: asf_bar(sheet_count, component_count), asw_bar(sheet_count)
    real(dp) :: bv_bar(volume_count, test_count, component_count)
    real(dp) :: bf_bar(volume_count, component_count), bw_bar(volume_count)
    real(dp) :: bb_bar(boundary_count, test_count, component_count)
    real(dp) :: bbf_bar(boundary_count, component_count), bbw_bar(boundary_count)
    real(dp) :: bs_bar(sheet_count, test_count, component_count)
    real(dp) :: bsf_bar(sheet_count, component_count), bsw_bar(sheet_count)
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    av = reshape([1.0_dp, -0.4_dp, 0.7_dp, 0.2_dp, -0.3_dp, 0.8_dp, 0.5_dp, 1.2_dp], &
        shape(av))
    af = reshape([0.8_dp, -1.1_dp, 0.6_dp, 1.4_dp], shape(af))
    aw = [0.7_dp, 1.3_dp]
    ab = reshape([0.6_dp, -0.2_dp, 1.1_dp, 0.4_dp], shape(ab))
    abf = reshape([0.3_dp, -0.9_dp], shape(abf))
    abw = [0.8_dp]
    as = reshape([-0.7_dp, 0.2_dp, 0.5_dp, 1.0_dp], shape(as))
    asf = reshape([1.2_dp, -0.6_dp], shape(asf))
    asw = [0.9_dp]
    bv = 0.8_dp*av + 0.1_dp
    bf = 0.9_dp*af - 0.2_dp
    bw = [1.1_dp, 0.6_dp]
    bb = 0.7_dp*ab - 0.3_dp
    bbf = 1.2_dp*abf + 0.1_dp
    bbw = [0.5_dp]
    bs = 1.1_dp*as + 0.2_dp
    bsf = 0.8_dp*asf - 0.1_dp
    bsw = [1.4_dp]
    av_dot = 0.03_dp
    af_dot = -0.02_dp
    aw_dot = [0.04_dp, -0.01_dp]
    ab_dot = -0.01_dp
    abf_dot = 0.02_dp
    abw_dot = [0.03_dp]
    as_dot = 0.02_dp
    asf_dot = -0.03_dp
    asw_dot = [-0.02_dp]
    bv_dot = -0.02_dp
    bf_dot = 0.04_dp
    bw_dot = [-0.03_dp, 0.02_dp]
    bb_dot = 0.01_dp
    bbf_dot = -0.02_dp
    bbw_dot = [0.02_dp]
    bs_dot = -0.03_dp
    bsf_dot = 0.01_dp
    bsw_dot = [0.03_dp]

    call evaluate_force_balance_representation_parity( &
        av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, &
        bs, bsf, bsw, difference, status)
    call build_reference(expected)
    call record_condition(status%code == 0 .and. maxval(abs(difference - expected)) < 1.0e-14_dp, &
        "force-form parity matches an independent nested-loop oracle")

    call evaluate_force_balance_representation_parity_jvp( &
        av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, bs, bsf, bsw, &
        av_dot, af_dot, aw_dot, ab_dot, abf_dot, abw_dot, as_dot, asf_dot, asw_dot, &
        bv_dot, bf_dot, bw_dot, bb_dot, bbf_dot, bbw_dot, bs_dot, bsf_dot, bsw_dot, &
        difference_dot, status)
    call evaluate_force_balance_representation_parity( &
        av + step*av_dot, af + step*af_dot, aw + step*aw_dot, ab + step*ab_dot, &
        abf + step*abf_dot, abw + step*abw_dot, as + step*as_dot, asf + step*asf_dot, &
        asw + step*asw_dot, bv + step*bv_dot, bf + step*bf_dot, bw + step*bw_dot, &
        bb + step*bb_dot, bbf + step*bbf_dot, bbw + step*bbw_dot, bs + step*bs_dot, &
        bsf + step*bsf_dot, bsw + step*bsw_dot, difference_plus, status)
    call evaluate_force_balance_representation_parity( &
        av - step*av_dot, af - step*af_dot, aw - step*aw_dot, ab - step*ab_dot, &
        abf - step*abf_dot, abw - step*abw_dot, as - step*as_dot, asf - step*asf_dot, &
        asw - step*asw_dot, bv - step*bv_dot, bf - step*bf_dot, bw - step*bw_dot, &
        bb - step*bb_dot, bbf - step*bbf_dot, bbw - step*bbw_dot, bs - step*bs_dot, &
        bsf - step*bsf_dot, bsw - step*bsw_dot, difference_minus, status)
    call record_condition(status%code == 0 .and. maxval(abs(difference_dot - &
        (difference_plus - difference_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "force-form parity JVP matches an independent central difference")

    difference_bar = reshape([0.4_dp, -0.7_dp, 0.2_dp, 0.9_dp], shape(difference_bar))
    call evaluate_force_balance_representation_parity_vjp( &
        av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, bs, bsf, bsw, &
        difference_bar, av_bar, af_bar, aw_bar, ab_bar, abf_bar, abw_bar, as_bar, asf_bar, &
        asw_bar, bv_bar, bf_bar, bw_bar, bb_bar, bbf_bar, bbw_bar, bs_bar, bsf_bar, bsw_bar, status)
    lhs = sum(difference_bar*difference_dot)
    rhs = sum(av_bar*av_dot) + sum(af_bar*af_dot) + sum(aw_bar*aw_dot) + &
        sum(ab_bar*ab_dot) + sum(abf_bar*abf_dot) + sum(abw_bar*abw_dot) + &
        sum(as_bar*as_dot) + sum(asf_bar*asf_dot) + sum(asw_bar*asw_dot) + &
        sum(bv_bar*bv_dot) + sum(bf_bar*bf_dot) + sum(bw_bar*bw_dot) + &
        sum(bb_bar*bb_dot) + sum(bbf_bar*bbf_dot) + sum(bbw_bar*bbw_dot) + &
        sum(bs_bar*bs_dot) + sum(bsf_bar*bsf_dot) + sum(bsw_bar*bsw_dot)
    call record_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "force-form parity VJP satisfies the real dot-product identity")

    bw(1) = -1.0_dp
    call evaluate_force_balance_representation_parity( &
        av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, &
        bs, bsf, bsw, difference, status)
    call record_condition(status%code /= 0, &
        "force-form parity rejects a nonpositive representation measure")
    call check_summary("force-balance representation parity")
    if (.not. all_passed) error stop 1

contains

    subroutine build_reference(result)
        real(dp), intent(out) :: result(:, :)
        integer :: sample, test_function, component

        result = 0.0_dp
        do sample = 1, volume_count
            do test_function = 1, test_count
                do component = 1, component_count
                    result(test_function, component) = result(test_function, component) + &
                        aw(sample)*av(sample, test_function, component)*af(sample, component) - &
                        bw(sample)*bv(sample, test_function, component)*bf(sample, component)
                end do
            end do
        end do
        call add_term(ab, abf, abw, result, 1.0_dp)
        call add_term(bb, bbf, bbw, result, -1.0_dp)
        call add_term(as, asf, asw, result, 1.0_dp)
        call add_term(bs, bsf, bsw, result, -1.0_dp)
    end subroutine build_reference

    subroutine add_term(test, force, weights, result, sign)
        real(dp), intent(in) :: test(:, :, :), force(:, :), weights(:), sign
        real(dp), intent(inout) :: result(:, :)
        integer :: sample, test_function, component

        do sample = 1, size(weights)
            do test_function = 1, test_count
                do component = 1, component_count
                    result(test_function, component) = result(test_function, component) + sign* &
                        weights(sample)*test(sample, test_function, component)*force(sample, component)
                end do
            end do
        end do
    end subroutine add_term

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_force_balance_representation_parity
