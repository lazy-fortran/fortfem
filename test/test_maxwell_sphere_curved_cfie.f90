program test_maxwell_sphere_curved_cfie
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: cfie(:, :), efie(:, :), mfie(:, :)
    complex(dp), allocatable :: product(:, :), regularizer(:, :)
    complex(dp), allocatable :: scaled_cfie(:, :), scaled_efie(:, :)
    complex(dp), allocatable :: scaled_mfie(:, :), scaled_product(:, :)
    complex(dp), allocatable :: scaled_regularizer(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: scaled_vertices(:, :), vertices(:, :)
    real(dp) :: error
    integer :: status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, vertices, triangles)
    call assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d( &
        vertices, triangles, 1.0_dp, 0.65_dp, 1.5_dp, 2, 2.0e-4_dp, 1, &
        0.18_dp, cfie, efie, mfie, regularizer, product, status)
    call record_condition(status == 0 .and. maxval(abs(cfie)) > 0.0_dp, &
        "curved regularized CFIE assembles")
    call record_condition(maxval(abs(cfie - (mfie - product))) < &
        3.0e-13_dp*max(1.0_dp, maxval(abs(cfie))), &
        "curved CFIE has the provenance-faithful product algebra")

    call generate_sphere_surface_mesh(2.0_dp, 0, scaled_vertices, triangles)
    call assemble_maxwell_sphere_curved_regularized_cfie_rwg_3d( &
        scaled_vertices, triangles, 2.0_dp, 0.325_dp, 1.5_dp, 2, 2.0e-4_dp, 1, &
        0.18_dp, scaled_cfie, scaled_efie, scaled_mfie, scaled_regularizer, &
        scaled_product, status)
    error = maxval(abs(scaled_cfie - 2.0_dp*cfie))/maxval(abs(scaled_cfie))
    call record_condition(status == 0 .and. error < 3.0e-9_dp, &
        "curved regularized CFIE obeys analytical BC length scaling")
    call record_condition(maxval(abs(scaled_regularizer - regularizer)) < &
        3.0e-9_dp*max(1.0_dp, maxval(abs(regularizer))), &
        "curved BC regularizer remains similarity invariant")

    call check_summary("Curved-sphere regularized Maxwell CFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_sphere_curved_cfie
