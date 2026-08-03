program test_maxwell_torus_curved_cfie_slow
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_torus_curved_regularized_cfie_rwg_3d
    use fortfem_core, only: generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    complex(dp), allocatable :: cfie(:, :), efie(:, :), mfie(:, :)
    complex(dp), allocatable :: product(:, :), regularizer(:, :)
    complex(dp), allocatable :: scaled_cfie(:, :), scaled_efie(:, :)
    complex(dp), allocatable :: scaled_mfie(:, :), scaled_product(:, :)
    complex(dp), allocatable :: scaled_regularizer(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), scaled_vertices(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 3, 3, vertices, triangles, parameters)
    call assemble_maxwell_torus_curved_regularized_cfie_rwg_3d( &
        vertices, triangles, parameters, major_radius, minor_radius, 0.45_dp, &
        1.5_dp, 1, 3.0e-4_dp, 1, 0.18_dp, cfie, efie, mfie, regularizer, &
        product, status)
    call record_condition(status == 0 .and. maxval(abs(cfie)) > 0.0_dp, &
        "exact-torus regularized Maxwell CFIE assembles")
    call record_condition(maxval(abs(cfie - (mfie - product))) < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(cfie))), &
        "exact-torus CFIE has the provenance-faithful product algebra")

    scaled_vertices = 2.0_dp*vertices
    call assemble_maxwell_torus_curved_regularized_cfie_rwg_3d( &
        scaled_vertices, triangles, parameters, 2.0_dp*major_radius, &
        2.0_dp*minor_radius, 0.225_dp, 1.5_dp, 1, 3.0e-4_dp, 1, 0.18_dp, &
        scaled_cfie, scaled_efie, scaled_mfie, scaled_regularizer, &
        scaled_product, status)
    error = maxval(abs(scaled_cfie - 2.0_dp*cfie))/maxval(abs(scaled_cfie))
    call record_condition(status == 0 .and. error < 4.0e-11_dp, &
        "exact-torus regularized CFIE obeys analytical BC length scaling")
    call record_condition(maxval(abs(scaled_regularizer - regularizer)) < &
        4.0e-11_dp*max(1.0_dp, maxval(abs(regularizer))), &
        "exact-torus BC regularizer remains similarity invariant")

    call check_summary("Exact-curved torus regularized Maxwell CFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_cfie_slow
