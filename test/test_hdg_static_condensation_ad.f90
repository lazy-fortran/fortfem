program test_hdg_static_condensation_ad
    use check, only: check_condition, check_summary
    use fortfem_feec, only: assemble_hdg_static_condensation, &
        assemble_hdg_static_condensation_jvp, assemble_hdg_static_condensation_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    real(dp), parameter :: interior_matrix(2, 2) = reshape([3.0_dp, 0.5_dp, 1.0_dp, 2.0_dp], [2, 2])
    real(dp), parameter :: interior_to_trace(2, 1) = reshape([0.4_dp, -0.7_dp], [2, 1])
    real(dp), parameter :: trace_to_interior(1, 2) = reshape([1.2_dp, -0.3_dp], [1, 2])
    real(dp), parameter :: trace_matrix(1, 1) = reshape([2.5_dp], [1, 1])
    real(dp), parameter :: interior_rhs(2) = [1.1_dp, -0.8_dp]
    real(dp), parameter :: trace_rhs(1) = [0.6_dp]
    real(dp), parameter :: interior_matrix_dot(2, 2) = reshape([ &
        0.1_dp, -0.2_dp, 0.3_dp, 0.4_dp], [2, 2])
    real(dp), parameter :: interior_to_trace_dot(2, 1) = reshape([0.05_dp, 0.2_dp], [2, 1])
    real(dp), parameter :: trace_to_interior_dot(1, 2) = reshape([-0.1_dp, 0.3_dp], [1, 2])
    real(dp), parameter :: trace_matrix_dot(1, 1) = reshape([0.2_dp], [1, 1])
    real(dp), parameter :: interior_rhs_dot(2) = [0.2_dp, -0.1_dp]
    real(dp), parameter :: trace_rhs_dot(1) = [-0.15_dp]
    real(dp), parameter :: condensed_matrix_bar(1, 1) = reshape([0.7_dp], [1, 1])
    real(dp), parameter :: condensed_rhs_bar(1) = [-0.4_dp]
    real(dp) :: condensed_matrix(1, 1), condensed_rhs(1)
    real(dp) :: condensed_matrix_dot(1, 1), condensed_rhs_dot(1)
    real(dp) :: interior_matrix_bar(2, 2), interior_to_trace_bar(2, 1)
    real(dp) :: trace_to_interior_bar(1, 2), trace_matrix_bar(1, 1)
    real(dp) :: interior_rhs_bar(2), trace_rhs_bar(1), lhs, rhs
    type(fortsparse_status_t) :: status

    call assemble_hdg_static_condensation( &
        interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
        interior_rhs, trace_rhs, condensed_matrix, condensed_rhs, status)
    call check_condition(status%code == 0, &
        "HDG static condensation accepts a nonsingular local block")
    call check_condition(abs(condensed_matrix(1, 1) - 2.0472727272727273_dp) < 1.0e-13_dp .and. &
        abs(condensed_rhs(1) + 0.2154545454545455_dp) < 1.0e-13_dp, &
        "HDG static condensation matches the Schur-complement oracle")

    call assemble_hdg_static_condensation_jvp( &
        interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
        interior_rhs, trace_rhs, interior_matrix_dot, interior_to_trace_dot, &
        trace_to_interior_dot, trace_matrix_dot, interior_rhs_dot, trace_rhs_dot, &
        condensed_matrix_dot, condensed_rhs_dot, status)
    call check_condition(status%code == 0, &
        "HDG static condensation JVP accepts local block directions")

    call assemble_hdg_static_condensation_vjp( &
        interior_matrix, interior_to_trace, trace_to_interior, trace_matrix, &
        interior_rhs, trace_rhs, condensed_matrix_bar, condensed_rhs_bar, &
        interior_matrix_bar, interior_to_trace_bar, trace_to_interior_bar, &
        trace_matrix_bar, interior_rhs_bar, trace_rhs_bar, status)
    lhs = sum(condensed_matrix_bar*condensed_matrix_dot) + &
        dot_product(condensed_rhs_bar, condensed_rhs_dot)
    rhs = sum(interior_matrix_bar*interior_matrix_dot) + &
        sum(interior_to_trace_bar*interior_to_trace_dot) + &
        sum(trace_to_interior_bar*trace_to_interior_dot) + &
        sum(trace_matrix_bar*trace_matrix_dot) + &
        dot_product(interior_rhs_bar, interior_rhs_dot) + &
        dot_product(trace_rhs_bar, trace_rhs_dot)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "HDG static condensation VJP satisfies the real dot-product identity")

    call assemble_hdg_static_condensation( &
        reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2]), interior_to_trace, &
        trace_to_interior, trace_matrix, interior_rhs, trace_rhs, condensed_matrix, &
        condensed_rhs, status)
    call check_condition(status%code /= 0, &
        "HDG static condensation rejects a singular interior block")
    call check_summary("HDG static condensation AD")
end program test_hdg_static_condensation_ad
