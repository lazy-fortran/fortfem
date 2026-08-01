module fortfem_enrichment_support_activation
    !! Fixed-topology XFEM support activation from a CSR support map.
    !!
    !! A basis support is active when its level-set samples contain both signs.
    !! The logical activation map is discrete.  The returned extrema and
    !! signed margin are differentiable only while their owners and the margin
    !! branch remain fixed; the JVP/VJP reject ties and changed owners.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    public :: evaluate_enrichment_support_activation
    public :: evaluate_enrichment_support_activation_jvp
    public :: evaluate_enrichment_support_activation_vjp

contains

    subroutine evaluate_enrichment_support_activation( &
            level_values, support_offsets, support_nodes, active_mask, &
            support_min, support_max, activation_margin, min_owner, max_owner, &
            margin_branch, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        logical, intent(out) :: active_mask(:)
        real(dp), intent(out) :: support_min(:), support_max(:)
        real(dp), intent(out) :: activation_margin(:)
        integer, intent(out) :: min_owner(:), max_owner(:), margin_branch(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: basis, entry, first, last

        active_mask = .false.
        support_min = 0.0_dp
        support_max = 0.0_dp
        activation_margin = 0.0_dp
        min_owner = 0
        max_owner = 0
        margin_branch = 0
        call validate_activation_inputs( &
            level_values, support_offsets, support_nodes, active_mask, &
            support_min, support_max, activation_margin, min_owner, max_owner, &
            margin_branch, status)
        if (status%code /= FORTSPARSE_OK) return

        do basis = 1, size(active_mask)
            first = support_offsets(basis)
            last = support_offsets(basis + 1) - 1
            min_owner(basis) = support_nodes(first)
            max_owner(basis) = support_nodes(first)
            do entry = first + 1, last
                if (level_values(support_nodes(entry)) < &
                    level_values(min_owner(basis))) then
                    min_owner(basis) = support_nodes(entry)
                end if
                if (level_values(support_nodes(entry)) > &
                    level_values(max_owner(basis))) then
                    max_owner(basis) = support_nodes(entry)
                end if
            end do
            support_min(basis) = level_values(min_owner(basis))
            support_max(basis) = level_values(max_owner(basis))
            if (support_min(basis) < 0.0_dp) then
                if (support_max(basis) > 0.0_dp) active_mask(basis) = .true.
            end if
            call select_margin_branch( &
                support_min(basis), support_max(basis), &
                activation_margin(basis), margin_branch(basis))
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_activation

    subroutine evaluate_enrichment_support_activation_jvp( &
            level_values, support_offsets, support_nodes, min_owner, max_owner, &
            margin_branch, level_dot, support_min_dot, support_max_dot, &
            activation_margin_dot, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        integer, intent(in) :: min_owner(:), max_owner(:), margin_branch(:)
        real(dp), intent(in) :: level_dot(:)
        real(dp), intent(out) :: support_min_dot(:), support_max_dot(:)
        real(dp), intent(out) :: activation_margin_dot(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: basis

        support_min_dot = 0.0_dp
        support_max_dot = 0.0_dp
        activation_margin_dot = 0.0_dp
        call validate_derivative_inputs( &
            level_values, support_offsets, support_nodes, min_owner, max_owner, &
            margin_branch, level_dot, support_min_dot, support_max_dot, &
            activation_margin_dot, status)
        if (status%code /= FORTSPARSE_OK) return

        do basis = 1, size(min_owner)
            if (.not. fixed_extrema( &
                level_values, support_offsets, support_nodes, basis, &
                min_owner(basis), max_owner(basis), &
                margin_branch(basis))) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "support activation JVP crossed an owner or margin tie")
                return
            end if
            support_min_dot(basis) = level_dot(min_owner(basis))
            support_max_dot(basis) = level_dot(max_owner(basis))
            if (margin_branch(basis) == 1) then
                activation_margin_dot(basis) = support_max_dot(basis)
            else
                activation_margin_dot(basis) = -support_min_dot(basis)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_activation_jvp

    subroutine evaluate_enrichment_support_activation_vjp( &
            level_values, support_offsets, support_nodes, min_owner, max_owner, &
            margin_branch, support_min_bar, support_max_bar, &
            activation_margin_bar, level_bar, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        integer, intent(in) :: min_owner(:), max_owner(:), margin_branch(:)
        real(dp), intent(in) :: support_min_bar(:), support_max_bar(:)
        real(dp), intent(in) :: activation_margin_bar(:)
        real(dp), intent(out) :: level_bar(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: basis

        level_bar = 0.0_dp
        call validate_vjp_inputs( &
            level_values, support_offsets, support_nodes, min_owner, max_owner, &
            margin_branch, support_min_bar, support_max_bar, &
            activation_margin_bar, level_bar, status)
        if (status%code /= FORTSPARSE_OK) return

        do basis = 1, size(min_owner)
            if (.not. fixed_extrema( &
                level_values, support_offsets, support_nodes, basis, &
                min_owner(basis), max_owner(basis), &
                margin_branch(basis))) then
                call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                    "support activation VJP crossed an owner or margin tie")
                return
            end if
            level_bar(min_owner(basis)) = level_bar(min_owner(basis)) + &
                support_min_bar(basis)
            level_bar(max_owner(basis)) = level_bar(max_owner(basis)) + &
                support_max_bar(basis)
            if (margin_branch(basis) == 1) then
                level_bar(max_owner(basis)) = level_bar(max_owner(basis)) + &
                    activation_margin_bar(basis)
            else
                level_bar(min_owner(basis)) = level_bar(min_owner(basis)) - &
                    activation_margin_bar(basis)
            end if
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine evaluate_enrichment_support_activation_vjp

    subroutine validate_activation_inputs( &
            level_values, support_offsets, support_nodes, active_mask, &
            support_min, support_max, activation_margin, min_owner, max_owner, &
            margin_branch, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        logical, intent(in) :: active_mask(:)
        real(dp), intent(in) :: support_min(:), support_max(:)
        real(dp), intent(in) :: activation_margin(:)
        integer, intent(in) :: min_owner(:), max_owner(:), margin_branch(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_structure(level_values, support_offsets, support_nodes, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(active_mask) /= size(support_offsets) - 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation has an incompatible mask")
            return
        end if
        if (size(support_min) /= size(active_mask)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation has incompatible extrema")
            return
        end if
        if (size(support_max) /= size(active_mask) .or. &
            size(activation_margin) /= size(active_mask)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation has incompatible margin arrays")
            return
        end if
        if (size(min_owner) /= size(active_mask) .or. &
            size(max_owner) /= size(active_mask) .or. &
            size(margin_branch) /= size(active_mask)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation has incompatible owner arrays")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_activation_inputs

    subroutine validate_derivative_inputs( &
            level_values, support_offsets, support_nodes, min_owner, max_owner, &
            margin_branch, level_dot, support_min_dot, support_max_dot, &
            activation_margin_dot, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        integer, intent(in) :: min_owner(:), max_owner(:), margin_branch(:)
        real(dp), intent(in) :: level_dot(:)
        real(dp), intent(in) :: support_min_dot(:), support_max_dot(:)
        real(dp), intent(in) :: activation_margin_dot(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_structure(level_values, support_offsets, support_nodes, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(min_owner) /= size(support_offsets) - 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has incompatible owners")
            return
        end if
        if (size(max_owner) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has incompatible branches")
            return
        end if
        if (size(margin_branch) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has incompatible branches")
            return
        end if
        if (size(level_dot) /= size(level_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has an incompatible level increment")
            return
        end if
        if (size(support_min_dot) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has incompatible outputs")
            return
        end if
        if (size(support_max_dot) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has incompatible outputs")
            return
        end if
        if (size(activation_margin_dot) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has incompatible outputs")
            return
        end if
        if (any(min_owner < 1)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has an invalid minimum owner")
            return
        end if
        if (any(min_owner > size(level_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has an invalid minimum owner")
            return
        end if
        if (any(max_owner < 1)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has an invalid maximum owner")
            return
        end if
        if (any(max_owner > size(level_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP has an invalid maximum owner")
            return
        end if
        if (any(.not. ieee_is_finite(level_dot))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation JVP received non-finite data")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_derivative_inputs

    subroutine validate_vjp_inputs( &
            level_values, support_offsets, support_nodes, min_owner, max_owner, &
            margin_branch, support_min_bar, support_max_bar, &
            activation_margin_bar, level_bar, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        integer, intent(in) :: min_owner(:), max_owner(:), margin_branch(:)
        real(dp), intent(in) :: support_min_bar(:), support_max_bar(:)
        real(dp), intent(in) :: activation_margin_bar(:)
        real(dp), intent(in) :: level_bar(:)
        type(fortsparse_status_t), intent(out) :: status

        call validate_structure(level_values, support_offsets, support_nodes, status)
        if (status%code /= FORTSPARSE_OK) return
        if (size(min_owner) /= size(support_offsets) - 1) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has incompatible owners")
            return
        end if
        if (size(max_owner) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has incompatible branches")
            return
        end if
        if (size(margin_branch) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has incompatible branches")
            return
        end if
        if (size(support_min_bar) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has incompatible cotangents")
            return
        end if
        if (size(support_max_bar) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has incompatible cotangents")
            return
        end if
        if (size(activation_margin_bar) /= size(min_owner)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has incompatible cotangents")
            return
        end if
        if (size(level_bar) /= size(level_values)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has an incompatible level cotangent")
            return
        end if
        if (any(min_owner < 1)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has an invalid minimum owner")
            return
        end if
        if (any(min_owner > size(level_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has an invalid minimum owner")
            return
        end if
        if (any(max_owner < 1)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has an invalid maximum owner")
            return
        end if
        if (any(max_owner > size(level_values))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP has an invalid maximum owner")
            return
        end if
        if (any(.not. ieee_is_finite(support_min_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP received non-finite data")
            return
        end if
        if (any(.not. ieee_is_finite(support_max_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP received non-finite data")
            return
        end if
        if (any(.not. ieee_is_finite(activation_margin_bar))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "support activation VJP received non-finite data")
            return
        end if
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_vjp_inputs

    subroutine validate_structure(level_values, support_offsets, support_nodes, status)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        type(fortsparse_status_t), intent(out) :: status
        integer :: basis, entry, first, last, previous

        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "support activation has an invalid support map")
        if (size(level_values) < 1) return
        if (size(support_offsets) < 2) return
        if (size(support_nodes) < 1) return
        if (support_offsets(1) /= 1) return
        if (any(support_offsets < 1)) return
        if (support_offsets(size(support_offsets)) /= size(support_nodes) + 1) return
        if (any(support_offsets(2:) < &
            support_offsets(:size(support_offsets) - 1))) return
        if (any(support_offsets(2:) == &
            support_offsets(:size(support_offsets) - 1))) return
        if (any(support_nodes < 1)) return
        if (any(support_nodes > size(level_values))) return
        if (any(.not. ieee_is_finite(level_values))) return

        do basis = 1, size(support_offsets) - 1
            first = support_offsets(basis)
            last = support_offsets(basis + 1) - 1
            do entry = first, last
                if (level_values(support_nodes(entry)) == 0.0_dp) return
                do previous = first, entry - 1
                    if (support_nodes(previous) == support_nodes(entry)) return
                end do
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine validate_structure

    logical function fixed_extrema( &
            level_values, support_offsets, support_nodes, basis, min_node, &
            max_node, branch) result(valid)
        real(dp), intent(in) :: level_values(:)
        integer, intent(in) :: support_offsets(:), support_nodes(:)
        integer, intent(in) :: basis, min_node, max_node, branch
        integer :: entry, first, last, computed_min, computed_max
        integer :: min_ties, max_ties, expected_branch
        real(dp) :: min_value, max_value

        first = support_offsets(basis)
        last = support_offsets(basis + 1) - 1
        computed_min = support_nodes(first)
        computed_max = support_nodes(first)
        min_ties = 1
        max_ties = 1
        do entry = first + 1, last
            if (level_values(support_nodes(entry)) < &
                level_values(computed_min)) then
                computed_min = support_nodes(entry)
                min_ties = 1
            else if (level_values(support_nodes(entry)) == &
                    level_values(computed_min)) then
                min_ties = min_ties + 1
            end if
            if (level_values(support_nodes(entry)) > &
                level_values(computed_max)) then
                computed_max = support_nodes(entry)
                max_ties = 1
            else if (level_values(support_nodes(entry)) == &
                    level_values(computed_max)) then
                max_ties = max_ties + 1
            end if
        end do
        min_value = level_values(computed_min)
        max_value = level_values(computed_max)
        if (max_value < -min_value) then
            expected_branch = 1
        else if (-min_value < max_value) then
            expected_branch = -1
        else
            expected_branch = 0
        end if
        valid = computed_min == min_node
        if (computed_max /= max_node) valid = .false.
        if (min_ties /= 1 .or. max_ties /= 1) valid = .false.
        if (expected_branch == 0 .or. branch /= expected_branch) valid = .false.
    end function fixed_extrema

    subroutine select_margin_branch(min_value, max_value, margin, branch)
        real(dp), intent(in) :: min_value, max_value
        real(dp), intent(out) :: margin
        integer, intent(out) :: branch

        if (max_value < -min_value) then
            margin = max_value
            branch = 1
        else if (-min_value < max_value) then
            margin = -min_value
            branch = -1
        else
            margin = max_value
            branch = 0
        end if
    end subroutine select_margin_branch

end module fortfem_enrichment_support_activation
