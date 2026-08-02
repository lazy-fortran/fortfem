program test_free_boundary_port_gallery_oracle
    !! Independent fast oracle for the neutral free-boundary-port gallery.
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_free_boundary_port_residual
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: sample_count = 4, component_count = 3
    real(dp), parameter :: tolerance = 2.0e-14_dp
    real(dp) :: physical_trace(sample_count, component_count)
    real(dp) :: external_target(sample_count, component_count)
    real(dp) :: sheet_current(sample_count, component_count)
    real(dp) :: weights(sample_count), residual(sample_count, component_count)
    real(dp) :: expected(sample_count, component_count)
    type(fortsparse_status_t) :: status
    integer :: sample, component

    physical_trace = reshape([ &
        1.0_dp, 0.5_dp, -0.2_dp, 0.8_dp, &
        -0.4_dp, 0.7_dp, 0.3_dp, -0.6_dp, &
        0.9_dp, -0.1_dp, 0.2_dp, 0.4_dp], shape(physical_trace))
    external_target = reshape([ &
        0.2_dp, 0.1_dp, -0.1_dp, 0.3_dp, &
        0.0_dp, 0.2_dp, -0.2_dp, 0.1_dp, &
        0.4_dp, 0.0_dp, -0.1_dp, 0.3_dp], shape(external_target))
    sheet_current = reshape([ &
        0.05_dp, -0.02_dp, 0.01_dp, 0.03_dp, &
        0.01_dp, 0.04_dp, -0.03_dp, 0.02_dp, &
        0.08_dp, -0.01_dp, 0.02_dp, -0.04_dp], shape(sheet_current))
    weights = [0.7_dp, 1.1_dp, 0.9_dp, 1.3_dp]
    do sample = 1, sample_count
        do component = 1, component_count
            expected(sample, component) = weights(sample)*( &
                physical_trace(sample, component) - &
                external_target(sample, component) - &
                sheet_current(sample, component))
        end do
    end do

    call assemble_free_boundary_port_residual(physical_trace, external_target, &
        weights, residual, status, sheet_current)
    call check_condition(status%code == 0, &
        "neutral gallery port accepts manufactured trace arrays")
    call check_condition(maxval(abs(residual - expected)) < tolerance, &
        "neutral gallery port matches independent weighted mismatch oracle")
    call check_condition(maxval(abs(residual(1, :) - residual(2, :))) > 1.0e-3_dp, &
        "neutral gallery oracle exercises nonuniform trace samples")
    call check_summary("free-boundary port gallery oracle")
end program test_free_boundary_port_gallery_oracle

