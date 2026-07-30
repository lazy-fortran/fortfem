program test_maxwell_efie_rwg_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_efie_rwg_3d, &
        assemble_maxwell_rwg_potential_operators_3d, &
        build_maxwell_rwg_surface_space
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: efie(:, :), scaled_efie(:, :)
    complex(dp), allocatable :: scalar_potential(:, :), vector_potential(:, :)
    real(dp), allocatable :: coefficients(:)
    real(dp) :: electric_field(3), vertices(3, 4)
    integer, allocatable :: edge_triangles(:, :), edge_vertices(:, :)
    integer :: boundary_triangles(3, 4), edge, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    boundary_triangles(:, 1) = [1, 3, 2]
    boundary_triangles(:, 2) = [1, 2, 4]
    boundary_triangles(:, 3) = [1, 4, 3]
    boundary_triangles(:, 4) = [2, 3, 4]
    call build_maxwell_rwg_surface_space( &
        vertices, boundary_triangles, edge_vertices, edge_triangles, status)
    call assemble_maxwell_rwg_potential_operators_3d( &
        vertices, boundary_triangles, 0.8_dp, 6, 1.0e-3_dp, 3, &
        vector_potential, scalar_potential, status)
    call record_condition(status == 0 .and. &
        maxval(abs(vector_potential - transpose(vector_potential))) < &
        2.0e-13_dp .and. &
        maxval(abs(scalar_potential - transpose(scalar_potential))) < &
        2.0e-13_dp, &
        "RWG Maxwell vector and scalar potential blocks are complex symmetric")

    electric_field = [0.7_dp, -0.4_dp, 0.2_dp]
    allocate(coefficients(size(edge_vertices, 2)))
    do edge = 1, size(edge_vertices, 2)
        coefficients(edge) = dot_product( &
            electric_field, &
            vertices(:, edge_vertices(2, edge)) - &
            vertices(:, edge_vertices(1, edge)))/norm2( &
            vertices(:, edge_vertices(2, edge)) - &
            vertices(:, edge_vertices(1, edge)))
    end do
    call record_condition(sqrt(sum(abs( &
        matmul(scalar_potential, coefficients))**2)) < 2.0e-13_dp, &
        "The charge block annihilates a divergence-free constant-field trace")

    call assemble_maxwell_efie_rwg_3d( &
        vertices, boundary_triangles, 0.8_dp, 1.7_dp, 6, 1.0e-3_dp, 3, &
        efie, status)
    call record_condition(sqrt(sum(abs(matmul(efie, coefficients) - &
        cmplx(0.0_dp, 0.8_dp*1.7_dp, dp)* &
        matmul(vector_potential, coefficients))**2)) < 2.0e-12_dp, &
        "EFIE has the exact vector-potential limit for solenoidal current")

    call assemble_maxwell_efie_rwg_3d( &
        2.0_dp*vertices, boundary_triangles, 0.4_dp, 1.7_dp, 6, &
        1.0e-3_dp, 3, scaled_efie, status)
    call record_condition(maxval(abs(scaled_efie - 4.0_dp*efie)) < 2.0e-8_dp, &
        "RWG EFIE obeys exact wave-scaled quadratic geometric scaling")

    call assemble_maxwell_efie_rwg_3d( &
        vertices, boundary_triangles, 0.0_dp, 1.7_dp, 6, 1.0e-3_dp, 3, &
        scaled_efie, status)
    call record_condition(status /= 0, &
        "RWG EFIE rejects its singular zero-frequency formulation")
    call check_summary("Three-dimensional Maxwell RWG EFIE")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_efie_rwg_3d
