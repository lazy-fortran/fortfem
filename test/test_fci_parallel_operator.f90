program test_fci_parallel_operator
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_fci_parallel_gradient_csc, &
        assemble_fci_parallel_support_divergence_csc
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: forward_map(2, 2, 2) = reshape([ &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2, 2])
    real(dp), parameter :: backward_map(2, 2, 2) = forward_map
    real(dp), parameter :: line_lengths(2, 2) = reshape([ &
        1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [2, 2])
    real(dp), parameter :: canonical_volumes(6) = [ &
        1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp]
    real(dp), parameter :: staggered_volumes(4) = [ &
        1.5_dp, 0.5_dp, 2.0_dp, 0.75_dp]
    real(dp), parameter :: field(6) = [ &
        0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: constant_field(6) = 7.0_dp
    real(dp), parameter :: flux(4) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    real(dp), parameter :: expected_divergence(6) = [ &
        1.5_dp, 0.5_dp, 4.5_dp, 2.0_dp, -3.0_dp, -3.0_dp]
    real(dp), allocatable :: gradient_field(:), constant_gradient(:)
    real(dp), allocatable :: divergence_flux(:)
    real(dp) :: weighted_pairing, adjoint_pairing
    type(csc_t) :: gradient, divergence
    type(fortsparse_status_t) :: status

    call assemble_fci_parallel_gradient_csc( &
        forward_map, backward_map, line_lengths, gradient, status)
    call check_condition(status%code == 0, &
        "FCI gradient assembly accepts two mapped plane interfaces")
    call check_condition(gradient%nrow == 4 .and. gradient%ncol == 6, &
        "FCI gradient has staggered rows and canonical columns")

    allocate(gradient_field(gradient%nrow), constant_gradient(gradient%nrow))
    gradient_field = csc_matvec(gradient, field)
    constant_gradient = csc_matvec(gradient, constant_field)
    call check_condition(maxval(abs(gradient_field - 1.0_dp)) < 1.0e-13_dp, &
        "FCI gradient reproduces a linear field along a straight line")
    call check_condition(maxval(abs(constant_gradient)) < 1.0e-13_dp, &
        "FCI gradient annihilates a constant field")

    call assemble_fci_parallel_support_divergence_csc( &
        gradient, canonical_volumes, staggered_volumes, divergence, status)
    call check_condition(status%code == 0 .and. divergence%nrow == 6 .and. &
        divergence%ncol == 4, &
        "FCI support divergence has the weighted adjoint shape")

    allocate(divergence_flux(divergence%nrow))
    divergence_flux = csc_matvec(divergence, flux)
    call check_condition(maxval(abs(divergence_flux - expected_divergence)) < &
        1.0e-13_dp, &
        "FCI support divergence matches the independent flux-balance oracle")

    weighted_pairing = dot_product( &
        canonical_volumes*divergence_flux, field)
    adjoint_pairing = -dot_product( &
        staggered_volumes*gradient_field, flux)
    call check_condition(abs(weighted_pairing - adjoint_pairing) < 1.0e-13_dp, &
        "FCI support divergence is the negative volume-weighted adjoint")

    call check_summary("PARALLAX-aligned FCI support operator")
end program test_fci_parallel_operator
