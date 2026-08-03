program test_maxwell_torus_curved_fem_bem_trace
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_maxwell_fem_bem_torus_curved_system_3d, &
        assemble_maxwell_torus_curved_rwg_mass_matrix, &
        assemble_maxwell_torus_curved_rwg_nedelec_coupling_3d
    use fortfem_core, only: generate_solid_torus_tetra_mesh
    use fortfem_feec, only: build_maxwell_rwg_surface_space, &
        map_maxwell_rwg_to_tetra_nedelec_edges
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp
    real(dp), parameter :: minor_radius = 0.65_dp
    integer, allocatable :: boundary_triangles(:, :), edge_triangles(:, :)
    integer, allocatable :: nedelec_dofs(:), rwg_edges(:, :), tetrahedra(:, :)
    real(dp), allocatable :: coupling(:, :), mass(:, :), parameters(:, :)
    real(dp), allocatable :: scales(:), vertices(:, :)
    complex(dp), allocatable :: system(:, :)
    real(dp) :: error, reciprocity_error
    integer :: edge, status
    logical :: all_passed

    all_passed = .true.
    call generate_solid_torus_tetra_mesh( &
        major_radius, minor_radius, 1, 4, 4, vertices, tetrahedra, &
        boundary_triangles, parameters)
    call build_maxwell_rwg_surface_space( &
        vertices, boundary_triangles, rwg_edges, edge_triangles, status)
    call map_maxwell_rwg_to_tetra_nedelec_edges( &
        vertices, tetrahedra, rwg_edges, nedelec_dofs, scales, status)
    call record_condition(status == 0 .and. all(nedelec_dofs > 0), &
        "solid-torus surface triangulation conforms to the volume mesh")
    call assemble_maxwell_torus_curved_rwg_mass_matrix( &
        vertices, boundary_triangles, parameters, major_radius, minor_radius, &
        8, mass, status)
    call assemble_maxwell_torus_curved_rwg_nedelec_coupling_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, 1, 8, coupling, status)
    error = 0.0_dp
    if (status == 0) then
        do edge = 1, size(rwg_edges, 2)
            error = max(error, maxval(abs( &
                coupling(:, nedelec_dofs(edge)) - &
                scales(edge)*mass(:, edge))))
        end do
    end if
    write (*, "(a,es12.4)") "exact-torus RWG-Nedelec trace error: ", error
    call record_condition(status == 0 .and. error < 3.0e-11_dp, &
        "exact torus covariant trace preserves RWG edge duality")
    call assemble_maxwell_fem_bem_torus_curved_system_3d( &
        vertices, tetrahedra, parameters, boundary_triangles, major_radius, &
        minor_radius, 1.1_dp, -0.7_dp, 0.8_dp, 1.4_dp, 3, 1.0e-5_dp, 1, &
        system, status, order=2)
    reciprocity_error = huge(1.0_dp)
    if (status == 0) &
        reciprocity_error = maxval(abs(system - transpose(system)))
    write (*, "(a,es12.4)") &
        "exact-torus system reciprocity error: ", reciprocity_error
    call record_condition(status == 0 .and. &
        reciprocity_error < 3.0e-10_dp, &
        "order-two exact-torus Maxwell FEM-BEM system is reciprocal")
    call check_summary("Exact-torus Maxwell FEM-BEM trace")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_maxwell_torus_curved_fem_bem_trace
