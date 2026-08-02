program test_equation_objective_registry
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        equation_objective_block_t, equation_objective_registry_t, &
        initialize_equation_objective_registry, &
        validate_equation_objective_registry, &
        equation_objective_registry_block, &
        equation_objective_registry_block_count, &
        equation_objective_registry_total_rows, &
        pack_equation_objective_values, &
        pack_equation_objective_values_jvp, &
        pack_equation_objective_values_vjp, &
        unpack_equation_objective_values, &
        REGISTRY_KIND_EQUATION, REGISTRY_KIND_OBJECTIVE, &
        REGISTRY_KIND_CONSTRAINT
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    character(len=18), parameter :: names(3) = [character(len=18) :: &
        "force_balance", "energy", "flux_constraint"]
    character(len=12), parameter :: units(3) = [character(len=12) :: &
        "N", "J", "Wb"]
    character(len=12), parameter :: normalizations(3) = [character(len=12) :: &
        "L2", "reference", "flux_scale"]
    character(len=20), parameter :: provenance(3) = [character(len=20) :: &
        "manufactured", "external-client", "external-client"]
    integer, parameter :: kinds(3) = [REGISTRY_KIND_EQUATION, &
        REGISTRY_KIND_OBJECTIVE, REGISTRY_KIND_CONSTRAINT]
    integer, parameter :: row_counts(3) = [2, 1, 3]
    logical, parameter :: active(3) = [.true., .true., .false.]
    logical, parameter :: fixed(3) = [.false., .false., .true.]
    real(dp), parameter :: eps = 1.0e-13_dp

    type(equation_objective_registry_t) :: registry, copied
    type(equation_objective_block_t) :: block
    type(fortsparse_status_t) :: status
    real(dp) :: values(3, 3), values_dot(3, 3)
    real(dp), allocatable :: values_bar(:, :), values_unpacked(:, :)
    real(dp), allocatable :: packed(:), packed_dot(:)
    real(dp) :: packed_bar(6)
    real(dp) :: packed_expected(6), block_bar_expected(3, 3)
    integer :: index
    logical :: all_passed

    all_passed = .true.

    call initialize_equation_objective_registry( &
        registry, names, kinds, row_counts, units, normalizations, &
        provenance, active, fixed, status)
    call record_condition(status%code == 0, &
        "registry accepts deterministic equation/objective/constraint metadata")
    call record_condition(validate_equation_objective_registry(registry, status) .and. &
        status%code == 0, "registry validates initialized metadata")
    call record_condition(equation_objective_registry_block_count(registry) == 3 .and. &
        equation_objective_registry_total_rows(registry) == 6, &
        "registry exposes block and row counts")

    call equation_objective_registry_block(registry, 1, block, status)
    call record_condition(status%code == 0 .and. block%name == names(1) .and. &
        block%kind == kinds(1) .and. block%row_offset == 1 .and. &
        block%row_count == 2 .and. block%active .and. .not. block%fixed .and. &
        block%units == units(1) .and. block%normalization == normalizations(1) .and. &
        block%provenance == provenance(1), "registry stores the first block")
    call equation_objective_registry_block(registry, 2, block, status)
    call record_condition(status%code == 0 .and. block%row_offset == 3, &
        "registry assigns the second block offset after the first count")
    call equation_objective_registry_block(registry, 3, block, status)
    call record_condition(status%code == 0 .and. block%row_offset == 4 .and. &
        .not. block%active .and. block%fixed, &
        "registry retains inactive and fixed flags without changing ordering")

    values = 0.0_dp
    values(:, 1) = [1.0_dp, 2.0_dp, 99.0_dp]
    values(:, 2) = [3.0_dp, 98.0_dp, 97.0_dp]
    values(:, 3) = [4.0_dp, 5.0_dp, 6.0_dp]
    packed_expected = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    call pack_equation_objective_values(registry, values, packed, status)
    call record_condition(status%code == 0 .and. all(abs(packed - &
        packed_expected) < eps), "registry packs rows in declared block order")
    call unpack_equation_objective_values( &
        registry, packed, values_unpacked, status)
    call record_condition(status%code == 0 .and. &
        all(abs(values_unpacked - reshape([ &
        1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
        4.0_dp, 5.0_dp, 6.0_dp], shape(values_unpacked))) < eps), &
        "registry unpacks rows and clears unused block padding")

    values_dot = reshape([ &
        0.1_dp, 0.2_dp, 8.0_dp, 0.3_dp, 7.0_dp, 6.0_dp, &
        0.4_dp, 0.5_dp, 0.6_dp], shape(values_dot))
    call pack_equation_objective_values_jvp( &
        registry, values_dot, packed_dot, status)
    call record_condition(status%code == 0 .and. all(abs(packed_dot - [ &
        0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp]) < eps), &
        "registry packing JVP is the deterministic linear action")

    packed_bar = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp, -6.0_dp]
    call pack_equation_objective_values_vjp( &
        registry, packed_bar, values_bar, status)
    block_bar_expected = 0.0_dp
    block_bar_expected(:, 1) = [1.0_dp, -2.0_dp, 0.0_dp]
    block_bar_expected(:, 2) = [3.0_dp, 0.0_dp, 0.0_dp]
    block_bar_expected(:, 3) = [-4.0_dp, 5.0_dp, -6.0_dp]
    call record_condition(status%code == 0 .and. all(abs(values_bar - &
        block_bar_expected) < eps), &
        "registry packing VJP scatters with zero padding")

    call pack_equation_objective_values( &
        registry, values, packed, status)
    call record_condition(status%code == 0 .and. &
        abs(dot_product(packed_bar, packed_dot) - &
        sum(values_bar * values_dot)) < eps, &
        "registry packing satisfies the independent transpose oracle")

    copied = registry
    copied%blocks(1)%name = "changed"
    copied%blocks(1)%row_count = 9
    call record_condition(registry%blocks(1)%name == names(1) .and. &
        registry%blocks(1)%row_count == row_counts(1), &
        "registry assignment performs a deep copy")

    call initialize_equation_objective_registry( &
        registry, [character(len=18) :: "duplicate", "duplicate"], &
        [REGISTRY_KIND_EQUATION, REGISTRY_KIND_OBJECTIVE], [1, 1], &
        [character(len=12) :: "1", "1"], [character(len=12) :: "a", "b"], &
        [character(len=20) :: "test", "test"], [.true., .true.], &
        [.false., .false.], status)
    call record_condition(status%code /= 0, "registry rejects duplicate block names")

    call initialize_equation_objective_registry( &
        registry, [character(len=18) :: "bad-kind"], [99], [1], &
        [character(len=12) :: "1"], [character(len=12) :: "a"], &
        [character(len=20) :: "test"], [.true.], [.false.], status)
    call record_condition(status%code /= 0, "registry rejects unknown block kinds")

    call initialize_equation_objective_registry( &
        registry, [character(len=18) :: "bad-count"], [REGISTRY_KIND_EQUATION], &
        [-1], [character(len=12) :: "1"], [character(len=12) :: "a"], &
        [character(len=20) :: "test"], [.true.], [.false.], status)
    call record_condition(status%code /= 0, "registry rejects negative row counts")

    call initialize_equation_objective_registry( &
        registry, [character(len=18) :: "bad-tags"], [REGISTRY_KIND_EQUATION], &
        [1], [character(len=12) :: ""], [character(len=12) :: "a"], &
        [character(len=20) :: "test"], [.true.], [.false.], status)
    call record_condition(status%code /= 0, "registry rejects missing units tags")

    call pack_equation_objective_values( &
        copied, values, packed, status)
    call record_condition(status%code /= 0 .and. &
        ieee_is_finite(packed(1)), &
        "registry rejects packing against incompatible copied metadata")

    call equation_objective_registry_block(registry, 0, block, status)
    call record_condition(status%code /= 0, "registry rejects an invalid block index")

    call check_summary("Equation/objective registry")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_equation_objective_registry
