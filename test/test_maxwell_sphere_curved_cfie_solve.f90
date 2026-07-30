program test_maxwell_sphere_curved_cfie_solve
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d, &
        assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d, &
        generate_sphere_surface_mesh, &
        solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: cfie(:, :), density(:), efie(:, :), mfie(:, :)
    complex(dp), allocatable :: product(:, :), regularizer(:, :), rhs(:)
    complex(dp), allocatable :: scaled_density(:)
    complex(dp) :: polarization(3)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: direction(3), error, residual
    integer :: status
    logical :: all_passed

    all_passed = .true.
    direction = [0.0_dp, 0.0_dp, 1.0_dp]
    polarization = cmplx([1.0_dp, 0.0_dp, 0.0_dp], 0.0_dp, dp)
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d( &
        vertices, triangles, 1.0_dp, direction, polarization, 0.6_dp, 1.3_dp, &
        2, 2.0e-4_dp, 1, 0.18_dp, density, status)
    call record_condition(status == 0 .and. size(density) > 0, &
        "curved regularized CFIE PEC solve succeeds")
    call assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d( &
        vertices, triangles, 1.0_dp, 0.6_dp, 1.3_dp, 2, 2.0e-4_dp, 1, &
        0.18_dp, cfie, efie, mfie, regularizer, product, status)
    call assemble_maxwell_sphere_curved_regularized_cfie_rhs_rwg_3d( &
        vertices, triangles, 1.0_dp, direction, polarization, 0.6_dp, 1.3_dp, &
        2, regularizer, rhs, status)
    residual = sqrt(sum(abs(matmul(cfie, density) - rhs)**2))/ &
        sqrt(sum(abs(rhs)**2))
    call record_condition(status == 0 .and. residual < 4.0e-12_dp, &
        "curved CFIE density satisfies the independently assembled system")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call solve_maxwell_pec_sphere_curved_regularized_cfie_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, direction, polarization, 0.3_dp, &
        1.3_dp, 2, 2.0e-4_dp, 1, 0.18_dp, scaled_density, status)
    error = maxval(abs(scaled_density - density))/maxval(abs(density))
    call record_condition(status == 0 .and. error < 5.0e-8_dp, &
        "curved CFIE current coefficients obey electromagnetic similarity")

    call check_summary("Curved-sphere regularized CFIE solve")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_cfie_solve
