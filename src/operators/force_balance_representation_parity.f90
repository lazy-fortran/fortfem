module fortfem_force_balance_representation_parity
    !! Fixed-topology parity actions for two force-balance representations.
    !!
    !! Each representation contains explicit volume, boundary, and sheet force
    !! samples.  The value action returns the residual of representation A
    !! minus representation B; the JVP and VJP use the same composition.  No
    !! force law or physical normalization is selected here.
    use fortfem_force_balance_residual, only: &
        assemble_force_balance_residual, assemble_force_balance_residual_jvp, &
        assemble_force_balance_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_force_balance_representation_parity
    public :: evaluate_force_balance_representation_parity_jvp
    public :: evaluate_force_balance_representation_parity_vjp

contains

    subroutine evaluate_force_balance_representation_parity( &
            av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, &
            bs, bsf, bsw, difference, status)
        real(dp), intent(in) :: av(:, :, :), af(:, :), aw(:)
        real(dp), intent(in) :: ab(:, :, :), abf(:, :), abw(:)
        real(dp), intent(in) :: as(:, :, :), asf(:, :), asw(:)
        real(dp), intent(in) :: bv(:, :, :), bf(:, :), bw(:)
        real(dp), intent(in) :: bb(:, :, :), bbf(:, :), bbw(:)
        real(dp), intent(in) :: bs(:, :, :), bsf(:, :), bsw(:)
        real(dp), intent(out) :: difference(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: residual_a(:, :), residual_b(:, :)

        difference = 0.0_dp
        if (.not. compatible_representations( &
            av, ab, as, bv, bb, bs, difference)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-form parity received incompatible row or component sizes")
            return
        end if
        allocate(residual_a(size(difference, 1), size(difference, 2)))
        allocate(residual_b(size(difference, 1), size(difference, 2)))
        call assemble_force_balance_residual( &
            av, af, aw, ab, abf, abw, as, asf, asw, residual_a, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_force_balance_residual( &
            bv, bf, bw, bb, bbf, bbw, bs, bsf, bsw, residual_b, status)
        if (status%code /= FORTSPARSE_OK) return
        difference = residual_a - residual_b
    end subroutine evaluate_force_balance_representation_parity

    subroutine evaluate_force_balance_representation_parity_jvp( &
            av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, &
            bs, bsf, bsw, av_dot, af_dot, aw_dot, ab_dot, abf_dot, abw_dot, &
            as_dot, asf_dot, asw_dot, bv_dot, bf_dot, bw_dot, bb_dot, bbf_dot, &
            bbw_dot, bs_dot, bsf_dot, bsw_dot, difference_dot, status)
        real(dp), intent(in) :: av(:, :, :), af(:, :), aw(:)
        real(dp), intent(in) :: ab(:, :, :), abf(:, :), abw(:)
        real(dp), intent(in) :: as(:, :, :), asf(:, :), asw(:)
        real(dp), intent(in) :: bv(:, :, :), bf(:, :), bw(:)
        real(dp), intent(in) :: bb(:, :, :), bbf(:, :), bbw(:)
        real(dp), intent(in) :: bs(:, :, :), bsf(:, :), bsw(:)
        real(dp), intent(in) :: av_dot(:, :, :), af_dot(:, :), aw_dot(:)
        real(dp), intent(in) :: ab_dot(:, :, :), abf_dot(:, :), abw_dot(:)
        real(dp), intent(in) :: as_dot(:, :, :), asf_dot(:, :), asw_dot(:)
        real(dp), intent(in) :: bv_dot(:, :, :), bf_dot(:, :), bw_dot(:)
        real(dp), intent(in) :: bb_dot(:, :, :), bbf_dot(:, :), bbw_dot(:)
        real(dp), intent(in) :: bs_dot(:, :, :), bsf_dot(:, :), bsw_dot(:)
        real(dp), intent(out) :: difference_dot(:, :)
        type(fortsparse_status_t), intent(out) :: status
        real(dp), allocatable :: residual_a_dot(:, :), residual_b_dot(:, :)

        difference_dot = 0.0_dp
        if (.not. compatible_representations( &
            av, ab, as, bv, bb, bs, difference_dot)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-form parity JVP received incompatible row or component sizes")
            return
        end if
        allocate(residual_a_dot(size(difference_dot, 1), size(difference_dot, 2)))
        allocate(residual_b_dot(size(difference_dot, 1), size(difference_dot, 2)))
        call assemble_force_balance_residual_jvp( &
            av, af, aw, ab, abf, abw, as, asf, asw, av_dot, af_dot, aw_dot, &
            ab_dot, abf_dot, abw_dot, as_dot, asf_dot, asw_dot, residual_a_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_force_balance_residual_jvp( &
            bv, bf, bw, bb, bbf, bbw, bs, bsf, bsw, bv_dot, bf_dot, bw_dot, &
            bb_dot, bbf_dot, bbw_dot, bs_dot, bsf_dot, bsw_dot, residual_b_dot, status)
        if (status%code /= FORTSPARSE_OK) return
        difference_dot = residual_a_dot - residual_b_dot
    end subroutine evaluate_force_balance_representation_parity_jvp

    subroutine evaluate_force_balance_representation_parity_vjp( &
            av, af, aw, ab, abf, abw, as, asf, asw, bv, bf, bw, bb, bbf, bbw, &
            bs, bsf, bsw, difference_bar, av_bar, af_bar, aw_bar, ab_bar, abf_bar, &
            abw_bar, as_bar, asf_bar, asw_bar, bv_bar, bf_bar, bw_bar, bb_bar, &
            bbf_bar, bbw_bar, bs_bar, bsf_bar, bsw_bar, status)
        real(dp), intent(in) :: av(:, :, :), af(:, :), aw(:)
        real(dp), intent(in) :: ab(:, :, :), abf(:, :), abw(:)
        real(dp), intent(in) :: as(:, :, :), asf(:, :), asw(:)
        real(dp), intent(in) :: bv(:, :, :), bf(:, :), bw(:)
        real(dp), intent(in) :: bb(:, :, :), bbf(:, :), bbw(:)
        real(dp), intent(in) :: bs(:, :, :), bsf(:, :), bsw(:)
        real(dp), intent(in) :: difference_bar(:, :)
        real(dp), intent(out) :: av_bar(:, :, :), af_bar(:, :), aw_bar(:)
        real(dp), intent(out) :: ab_bar(:, :, :), abf_bar(:, :), abw_bar(:)
        real(dp), intent(out) :: as_bar(:, :, :), asf_bar(:, :), asw_bar(:)
        real(dp), intent(out) :: bv_bar(:, :, :), bf_bar(:, :), bw_bar(:)
        real(dp), intent(out) :: bb_bar(:, :, :), bbf_bar(:, :), bbw_bar(:)
        real(dp), intent(out) :: bs_bar(:, :, :), bsf_bar(:, :), bsw_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        av_bar = 0.0_dp
        af_bar = 0.0_dp
        aw_bar = 0.0_dp
        ab_bar = 0.0_dp
        abf_bar = 0.0_dp
        abw_bar = 0.0_dp
        as_bar = 0.0_dp
        asf_bar = 0.0_dp
        asw_bar = 0.0_dp
        bv_bar = 0.0_dp
        bf_bar = 0.0_dp
        bw_bar = 0.0_dp
        bb_bar = 0.0_dp
        bbf_bar = 0.0_dp
        bbw_bar = 0.0_dp
        bs_bar = 0.0_dp
        bsf_bar = 0.0_dp
        bsw_bar = 0.0_dp
        if (.not. compatible_representations( &
            av, ab, as, bv, bb, bs, difference_bar)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "force-form parity VJP received incompatible row or component sizes")
            return
        end if
        call assemble_force_balance_residual_vjp( &
            av, af, aw, ab, abf, abw, as, asf, asw, difference_bar, av_bar, af_bar, &
            aw_bar, ab_bar, abf_bar, abw_bar, as_bar, asf_bar, asw_bar, status)
        if (status%code /= FORTSPARSE_OK) return
        call assemble_force_balance_residual_vjp( &
            bv, bf, bw, bb, bbf, bbw, bs, bsf, bsw, -difference_bar, bv_bar, bf_bar, &
            bw_bar, bb_bar, bbf_bar, bbw_bar, bs_bar, bsf_bar, bsw_bar, status)
    end subroutine evaluate_force_balance_representation_parity_vjp

    logical function compatible_representations( &
            av, ab, as, bv, bb, bs, difference) result(compatible)
        real(dp), intent(in) :: av(:, :, :), ab(:, :, :), as(:, :, :)
        real(dp), intent(in) :: bv(:, :, :), bb(:, :, :), bs(:, :, :)
        real(dp), intent(in) :: difference(:, :)

        compatible = size(av, 2) > 0 .and. size(av, 3) > 0 .and. &
            size(ab, 2) == size(av, 2) .and. size(ab, 3) == size(av, 3) .and. &
            size(as, 2) == size(av, 2) .and. size(as, 3) == size(av, 3) .and. &
            size(bv, 2) == size(av, 2) .and. size(bv, 3) == size(av, 3) .and. &
            size(bb, 2) == size(av, 2) .and. size(bb, 3) == size(av, 3) .and. &
            size(bs, 2) == size(av, 2) .and. size(bs, 3) == size(av, 3) .and. &
            size(difference, 1) == size(av, 2) .and. &
            size(difference, 2) == size(av, 3)
    end function compatible_representations

end module fortfem_force_balance_representation_parity
