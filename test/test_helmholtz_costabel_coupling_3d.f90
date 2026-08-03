program test_helmholtz_costabel_coupling_3d
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_helmholtz_fem_bem_costabel_3d, &
        solve_helmholtz_fem_bem_costabel_3d
    use fortfem_core, only: generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 0.7_dp
    real(dp), allocatable :: boundary_vertices(:, :), vertices(:, :)
    integer, allocatable :: tetrahedra(:, :), triangles(:, :)
    complex(dp), allocatable :: flux(:), load(:), matrix(:, :), potential(:)
    complex(dp) :: exact_flux, exact_potential
    real(dp) :: errors(2, 0:1), resonance_error, symmetry_error
    integer :: boundary_count, level, status
    logical :: all_passed

    all_passed = .true.
    exact_potential = exp(cmplx(0.0_dp, wave_number, dp))/ &
        (4.0_dp*acos(-1.0_dp))
    exact_flux = cmplx(-1.0_dp, wave_number, dp)*exact_potential
    do level = 0, 1
        call generate_sphere_surface_mesh( &
            1.0_dp, level, boundary_vertices, triangles)
        boundary_count = size(boundary_vertices, 2)
        allocate(vertices(3, boundary_count + 1))
        vertices(:, :boundary_count) = boundary_vertices
        vertices(:, boundary_count + 1) = 0.0_dp
        allocate(tetrahedra(4, size(triangles, 2)))
        tetrahedra(1, :) = boundary_count + 1
        tetrahedra(2:4, :) = triangles
        allocate(load(boundary_count + 1))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        load(boundary_count + 1) = cmplx(1.0_dp, 0.0_dp, dp)

        call assemble_helmholtz_fem_bem_costabel_3d( &
            vertices, tetrahedra, triangles, wave_number, wave_number, 8, &
            matrix, status)
        symmetry_error = maxval(abs(matrix - transpose(matrix)))
        call record_condition(status == 0 .and. symmetry_error < 3.0e-12_dp, &
            "3D Helmholtz Costabel matrix is complex symmetric")
        call solve_helmholtz_fem_bem_costabel_3d( &
            vertices, tetrahedra, triangles, load, wave_number, wave_number, &
            8, potential, flux, status)
        if (status /= 0) error stop "3D Helmholtz Costabel solve failed"
        errors(1, level) = maxval(abs( &
            potential(:boundary_count) - exact_potential))
        errors(2, level) = maxval(abs(flux - exact_flux))
        deallocate(boundary_vertices, vertices, tetrahedra, triangles, load)
        deallocate(matrix, potential, flux)
    end do

    call record_condition(all(errors(:, 1) < 0.8_dp*errors(:, 0)), &
        "3D Helmholtz Costabel traces converge under refinement")
    call record_condition(errors(1, 1) < 3.0e-2_dp, &
        "refined coupled solve matches the outgoing monopole potential")
    call record_condition(errors(2, 1) < 8.0e-2_dp, &
        "refined coupled solve matches the outgoing monopole flux")

    call generate_sphere_surface_mesh(1.0_dp, 1, boundary_vertices, triangles)
    boundary_count = size(boundary_vertices, 2)
    allocate(vertices(3, boundary_count + 1))
    vertices(:, :boundary_count) = boundary_vertices
    vertices(:, boundary_count + 1) = 0.0_dp
    allocate(tetrahedra(4, size(triangles, 2)))
    tetrahedra(1, :) = boundary_count + 1
    tetrahedra(2:4, :) = triangles
    allocate(load(boundary_count + 1))
    load = cmplx(0.0_dp, 0.0_dp, dp)
    load(boundary_count + 1) = cmplx(1.0_dp, 0.0_dp, dp)
    call solve_helmholtz_fem_bem_costabel_3d( &
        vertices, tetrahedra, triangles, load, acos(-1.0_dp), &
        acos(-1.0_dp), 8, potential, flux, status)
    exact_potential = -1.0_dp/(4.0_dp*acos(-1.0_dp))
    resonance_error = maxval(abs( &
        potential(:boundary_count) - exact_potential))
    call record_condition(status == 0 .and. resonance_error < 8.0e-2_dp, &
        "coupled solve remains accurate near the first spherical resonance")
    call check_summary("Three-dimensional Helmholtz Costabel FEM-BEM coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_costabel_coupling_3d
