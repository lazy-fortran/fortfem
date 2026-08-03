program test_fci_anisotropic_diffusion_field_vjp
    use check, only: check_condition, check_summary
    use fortfem_feec, only: apply_fci_anisotropic_diffusion, &
        apply_fci_anisotropic_diffusion_field_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(2, 2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2, 2])
    real(dp), parameter :: backward_map(2, 2, 2) = forward_map
    real(dp), parameter :: line_lengths(2, 2) = 1.0_dp
    real(dp), parameter :: parallel_coefficient(4) = 1.0_dp
    real(dp), parameter :: canonical_volumes(6) = 1.0_dp
    real(dp), parameter :: staggered_volumes(4) = 1.0_dp
    real(dp), parameter :: field(6) = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, &
        2.0_dp, 3.0_dp]
    real(dp), parameter :: diffusion_field_bar(6) = [ &
        0.3_dp, -0.5_dp, 0.7_dp, -0.9_dp, 1.1_dp, -1.3_dp]
    real(dp), parameter :: field_dot(6) = [ &
        -0.2_dp, 0.4_dp, -0.6_dp, 0.8_dp, -1.0_dp, 1.2_dp]
    real(dp), parameter :: expected_field_bar(6) = [ &
        -1.2_dp, 0.75_dp, -0.7_dp, 2.7_dp, -3.7_dp, 5.6_dp]
    real(dp) :: field_bar(6), diffusion_field(6), diffusion_field_dot(6)
    real(dp) :: lhs, rhs
    type(csc_t) :: perpendicular_operators(3)
    type(fortsparse_status_t) :: status

    call make_matrix( &
        [1, 2, 1, 2], [1, 1, 2, 2], [-2.0_dp, 2.0_dp, 0.5_dp, -2.0_dp], &
        perpendicular_operators(1), status)
    call make_matrix( &
        [1, 2], [1, 2], [-1.0_dp, -3.0_dp], perpendicular_operators(2), status)
    call make_matrix( &
        [1, 2], [1, 2], [-3.0_dp, -4.0_dp], perpendicular_operators(3), status)
    call apply_fci_anisotropic_diffusion_field_vjp( &
        perpendicular_operators, forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, &
        diffusion_field_bar, field_bar, status)
    call check_condition(status%code == 0, &
        "FCI anisotropic field VJP accepts square plane blocks")
    call check_condition(maxval(abs(field_bar - expected_field_bar)) < &
        1.0e-14_dp, "FCI anisotropic field VJP matches the explicit transpose oracle")

    call apply_fci_anisotropic_diffusion( &
        perpendicular_operators, forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, field_dot, &
        diffusion_field_dot, status)
    lhs = dot_product(diffusion_field_bar, diffusion_field_dot)
    rhs = dot_product(field_bar, field_dot)
    call check_condition(abs(lhs - rhs) < &
        1.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "FCI anisotropic field VJP satisfies the real dot-product identity")

    call apply_fci_anisotropic_diffusion_field_vjp( &
        perpendicular_operators(1:2), forward_map, backward_map, line_lengths, &
        parallel_coefficient, canonical_volumes, staggered_volumes, &
        diffusion_field_bar, field_bar, status)
    call check_condition(status%code /= 0, &
        "FCI anisotropic field VJP rejects a missing plane operator")
    call check_summary("FCI anisotropic diffusion field VJP")

contains

    subroutine make_matrix(rows, columns, values, matrix, status)
        integer, intent(in) :: rows(:), columns(:)
        real(dp), intent(in) :: values(:)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        call csc_from_triplet(2, 2, rows, columns, values, matrix, status)
    end subroutine make_matrix

end program test_fci_anisotropic_diffusion_field_vjp
