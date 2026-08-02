program test_api07_free_boundary_consumer
    !! Downstream-style smoke for the canonical free-boundary port facade.
    !!
    !! The client imports only ``fortfem_boundary`` (never the compatibility
    !! umbrella or the implementation module).  The expected residual,
    !! Jacobian-vector product, and transpose identity are written directly
    !! from r = w (physical - exterior - sheet), providing an independent
    !! behavioral oracle for this API-07 consumer slice.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_free_boundary_port_residual, &
        assemble_free_boundary_port_residual_jvp, &
        assemble_free_boundary_port_residual_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    call check_value_and_jvp()
    call check_vjp_transpose()
    call check_summary("API-07 free-boundary canonical consumer")

contains

    subroutine check_value_and_jvp()
        real(dp), parameter :: physical(3, 2) = reshape([ &
            1.2_dp, -0.4_dp, 0.7_dp, 0.5_dp, 1.1_dp, -0.2_dp], [3, 2])
        real(dp), parameter :: exterior(3, 2) = reshape([ &
            0.2_dp, 0.1_dp, -0.3_dp, -0.5_dp, 0.4_dp, 0.8_dp], [3, 2])
        real(dp), parameter :: sheet(3, 2) = reshape([ &
            0.1_dp, -0.2_dp, 0.05_dp, 0.2_dp, 0.3_dp, -0.1_dp], [3, 2])
        real(dp), parameter :: weights(3) = [0.5_dp, 1.25_dp, 2.0_dp]
        real(dp), parameter :: physical_dot(3, 2) = reshape([ &
            -0.3_dp, 0.2_dp, 0.4_dp, 0.1_dp, -0.5_dp, 0.6_dp], [3, 2])
        real(dp), parameter :: exterior_dot(3, 2) = reshape([ &
            0.2_dp, -0.1_dp, 0.3_dp, -0.4_dp, 0.5_dp, 0.7_dp], [3, 2])
        real(dp), parameter :: sheet_dot(3, 2) = reshape([ &
            0.05_dp, 0.2_dp, -0.15_dp, 0.3_dp, 0.1_dp, -0.2_dp], [3, 2])
        real(dp), parameter :: weights_dot(3) = [0.1_dp, -0.2_dp, 0.3_dp]
        real(dp) :: residual(3, 2), residual_dot(3, 2)
        real(dp) :: expected(3, 2), expected_dot(3, 2)
        integer :: q, component
        type(fortsparse_status_t) :: status

        do q = 1, 3
            do component = 1, 2
                expected(q, component) = weights(q)*(physical(q, component) - &
                    exterior(q, component) - sheet(q, component))
                expected_dot(q, component) = weights_dot(q)*( &
                    physical(q, component) - exterior(q, component) - &
                    sheet(q, component)) + weights(q)*(physical_dot(q, component) - &
                    exterior_dot(q, component) - sheet_dot(q, component))
            end do
        end do

        call assemble_free_boundary_port_residual( &
            physical, exterior, weights, residual, status, sheet)
        call check_condition(status%code == 0 .and. &
            maxval(abs(residual - expected)) < 2.0e-14_dp, &
            "free-boundary consumer evaluates the weighted sheet-aware residual")

        call assemble_free_boundary_port_residual_jvp( &
            physical, exterior, weights, physical_dot, exterior_dot, &
            weights_dot, residual_dot, status, sheet, sheet_dot)
        call check_condition(status%code == 0 .and. &
            maxval(abs(residual_dot - expected_dot)) < 2.0e-14_dp, &
            "free-boundary consumer exposes the product-rule JVP")
    end subroutine check_value_and_jvp

    subroutine check_vjp_transpose()
        real(dp), parameter :: physical(2, 2) = reshape([ &
            0.9_dp, -0.6_dp, 0.3_dp, 0.8_dp], [2, 2])
        real(dp), parameter :: exterior(2, 2) = reshape([ &
            -0.1_dp, 0.2_dp, 0.5_dp, -0.4_dp], [2, 2])
        real(dp), parameter :: sheet(2, 2) = reshape([ &
            0.15_dp, -0.05_dp, 0.2_dp, 0.1_dp], [2, 2])
        real(dp), parameter :: weights(2) = [0.75_dp, 1.5_dp]
        real(dp), parameter :: residual_bar(2, 2) = reshape([ &
            0.4_dp, -0.7_dp, 0.2_dp, 0.9_dp], [2, 2])
        real(dp), parameter :: physical_dot(2, 2) = reshape([ &
            -0.2_dp, 0.6_dp, 0.3_dp, -0.4_dp], [2, 2])
        real(dp), parameter :: exterior_dot(2, 2) = reshape([ &
            0.1_dp, -0.5_dp, 0.7_dp, 0.2_dp], [2, 2])
        real(dp), parameter :: sheet_dot(2, 2) = reshape([ &
            -0.3_dp, 0.4_dp, 0.05_dp, -0.1_dp], [2, 2])
        real(dp), parameter :: weights_dot(2) = [-0.15_dp, 0.25_dp]
        real(dp) :: residual_dot(2, 2), physical_bar(2, 2)
        real(dp) :: exterior_bar(2, 2), sheet_bar(2, 2), weights_bar(2)
        real(dp) :: lhs, rhs
        integer :: status_code
        type(fortsparse_status_t) :: status

        call assemble_free_boundary_port_residual_jvp( &
            physical, exterior, weights, physical_dot, exterior_dot, &
            weights_dot, residual_dot, status, sheet, sheet_dot)
        call check_condition(status%code == 0, &
            "free-boundary consumer accepts a valid directional state")

        call assemble_free_boundary_port_residual_vjp( &
            physical, exterior, weights, residual_bar, physical_bar, &
            exterior_bar, weights_bar, status, sheet, sheet_bar)
        status_code = status%code
        call check_condition(status_code == 0 .and. &
            maxval(abs(physical_bar + exterior_bar)) < 2.0e-14_dp .and. &
            maxval(abs(sheet_bar - exterior_bar)) < 2.0e-14_dp, &
            "free-boundary consumer returns the signed trace cotangents")

        lhs = sum(residual_bar*residual_dot)
        rhs = sum(physical_bar*physical_dot) + sum(exterior_bar*exterior_dot) + &
            dot_product(weights_bar, weights_dot) + sum(sheet_bar*sheet_dot)
        call check_condition(abs(lhs - rhs) < 2.0e-14_dp, &
            "free-boundary consumer satisfies the independent VJP transpose oracle")
    end subroutine check_vjp_transpose

end program test_api07_free_boundary_consumer
