program test_shifted_vector_enriched_space
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_feec, only: &
        evaluate_shifted_vector_enriched_space, &
        evaluate_shifted_vector_enriched_space_jvp, &
        evaluate_shifted_vector_enriched_space_vjp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: component_count = 2, basis_count = 3, point_count = 4
    real(dp) :: base_values(component_count, basis_count, point_count)
    real(dp) :: level_values(point_count), anchor_values(basis_count)
    real(dp) :: enriched_values(component_count, basis_count, point_count)
    real(dp) :: base_dot(component_count, basis_count, point_count)
    real(dp) :: level_dot(point_count), anchor_dot(basis_count)
    real(dp) :: enriched_dot(component_count, basis_count, point_count)
    real(dp) :: enriched_plus(component_count, basis_count, point_count)
    real(dp) :: enriched_bar(component_count, basis_count, point_count)
    real(dp) :: base_bar(component_count, basis_count, point_count)
    real(dp) :: level_bar(point_count), anchor_bar(basis_count)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    type(fortsparse_status_t) :: status
    logical :: all_passed

    all_passed = .true.
    call build_system(base_values, level_values, anchor_values)
    call evaluate_shifted_vector_enriched_space( &
        base_values, level_values, anchor_values, enriched_values, status)
    call record_condition(status%code == 0, &
        "shifted vector enriched space assembles")
    call record_condition(maxval(abs(enriched_values - expected_values( &
        base_values, level_values, anchor_values))) < 1.0e-14_dp, &
        "shifted vector enriched space matches its sign matrix oracle")

    base_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.03_dp, 0.04_dp, -0.03_dp, 0.02_dp, &
        0.05_dp, -0.01_dp, 0.02_dp, -0.04_dp, 0.03_dp, 0.01_dp, &
        -0.02_dp, 0.03_dp, -0.01_dp, 0.05_dp, 0.04_dp, -0.03_dp, &
        0.02_dp, 0.01_dp, -0.04_dp, 0.03_dp, -0.05_dp, 0.02_dp], &
        shape(base_dot))
    level_dot = [0.02_dp, -0.03_dp, 0.01_dp, 0.04_dp]
    anchor_dot = [-0.05_dp, 0.02_dp, 0.03_dp]
    call evaluate_shifted_vector_enriched_space_jvp( &
        base_values, level_values, anchor_values, base_dot, level_dot, &
        anchor_dot, enriched_dot, status)
    call record_condition(status%code == 0, &
        "shifted vector enriched space JVP assembles")
    epsilon = 1.0e-7_dp
    call evaluate_shifted_vector_enriched_space( &
        base_values + epsilon*base_dot, level_values + epsilon*level_dot, &
        anchor_values + epsilon*anchor_dot, enriched_plus, status)
    finite_difference_error = maxval(abs(enriched_dot - &
        (enriched_plus - enriched_values)/epsilon))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "shifted vector enriched space JVP matches a forward difference")

    enriched_bar = reshape([ &
        0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp, -0.3_dp, 0.5_dp, &
        -0.2_dp, 0.6_dp, -0.4_dp, 0.1_dp, 0.7_dp, -0.5_dp, &
        -0.6_dp, 0.2_dp, 0.4_dp, -0.7_dp, 0.3_dp, 0.5_dp, &
        0.1_dp, -0.4_dp, 0.6_dp, -0.2_dp, 0.8_dp, -0.3_dp], &
        shape(enriched_bar))
    call evaluate_shifted_vector_enriched_space_vjp( &
        base_values, level_values, anchor_values, enriched_bar, base_bar, &
        level_bar, anchor_bar, status)
    call record_condition(status%code == 0, &
        "shifted vector enriched space VJP assembles")
    lhs = sum(enriched_bar*enriched_dot)
    rhs = sum(base_bar*base_dot) + dot_product(level_bar, level_dot) + &
        dot_product(anchor_bar, anchor_dot)
    call record_condition(abs(lhs - rhs) < 2.0e-13_dp, &
        "shifted vector enriched space VJP satisfies the adjoint identity")

    anchor_values(1) = 0.0_dp
    call evaluate_shifted_vector_enriched_space( &
        base_values, level_values, anchor_values, enriched_values, status)
    call record_condition(status%code /= 0, &
        "topology-event vector anchors are rejected")

    call check_summary("shifted vector enriched space")
    if (.not. all_passed) error stop 1

contains

    subroutine build_system(base, levels, anchors)
        real(dp), intent(out) :: base(:, :, :), levels(:), anchors(:)

        base = reshape([ &
            0.2_dp, 0.4_dp, 0.1_dp, -0.3_dp, 0.5_dp, -0.2_dp, &
            0.7_dp, 0.6_dp, -0.4_dp, 0.8_dp, -0.1_dp, 0.3_dp, &
            -0.2_dp, 0.3_dp, 0.6_dp, 0.1_dp, -0.5_dp, 0.4_dp, &
            0.9_dp, -0.7_dp, 0.2_dp, 0.5_dp, 0.8_dp, -0.6_dp], shape(base))
        levels = [-0.7_dp, 0.4_dp, 1.1_dp, -0.2_dp]
        anchors = [-0.3_dp, 0.6_dp, -0.9_dp]
    end subroutine build_system

    function expected_values(base, levels, anchors) result(expected)
        real(dp), intent(in) :: base(:, :, :), levels(:), anchors(:)
        real(dp) :: expected(size(base, 1), size(base, 2), size(base, 3))
        integer :: component, basis, point

        do component = 1, size(base, 1)
            do basis = 1, size(base, 2)
                do point = 1, size(base, 3)
                    expected(component, basis, point) = base(component, basis, point)* &
                        (heaviside(levels(point)) - heaviside(anchors(basis)))
                end do
            end do
        end do
    end function expected_values

    pure real(dp) function heaviside(value)
        real(dp), intent(in) :: value

        if (value > 0.0_dp) then
            heaviside = 1.0_dp
        else
            heaviside = 0.0_dp
        end if
    end function heaviside

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_shifted_vector_enriched_space
