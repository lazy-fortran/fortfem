program test_laplace_costabel_coupling_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_laplace_fem_bem_costabel_3d, generate_sphere_surface_mesh, &
        solve_laplace_fem_bem_costabel_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: boundary_vertices(:, :), flux(:), load(:)
    real(dp), allocatable :: matrix(:, :), potential(:), vertices(:, :)
    integer, allocatable :: tetrahedra(:, :), triangles(:, :)
    real(dp) :: boundary_error, flux_balance, symmetry_error
    integer :: boundary_count, status
    logical :: all_passed

    all_passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 1, boundary_vertices, triangles)
    boundary_count = size(boundary_vertices, 2)
    allocate(vertices(3, boundary_count + 1))
    vertices(:, :boundary_count) = boundary_vertices
    vertices(:, boundary_count + 1) = 0.0_dp
    allocate(tetrahedra(4, size(triangles, 2)))
    tetrahedra(1, :) = boundary_count + 1
    tetrahedra(2:4, :) = triangles
    allocate(load(boundary_count + 1))
    load = 0.0_dp
    load(boundary_count + 1) = 1.0_dp

    call assemble_laplace_fem_bem_costabel_3d( &
        vertices, tetrahedra, triangles, 8, matrix, status)
    symmetry_error = maxval(abs(matrix - transpose(matrix)))
    call record_condition(status == 0 .and. symmetry_error < 2.0e-12_dp, &
        "Costabel FEM-BEM matrix is symmetric")

    call solve_laplace_fem_bem_costabel_3d( &
        vertices, tetrahedra, triangles, load, 8, potential, flux, status)
    flux_balance = abs(sum(flux*triangle_areas(vertices, triangles)) + 1.0_dp)
    boundary_error = maxval(abs( &
        potential(:boundary_count) - 1.0_dp/(4.0_dp*acos(-1.0_dp))))
    if (flux_balance >= 3.0e-2_dp) then
        write(*, '(a,es14.5,a,es14.5)') &
            "Costabel flux integral ", &
            sum(flux*triangle_areas(vertices, triangles)), &
            " balance error ", flux_balance
    end if
    call record_condition(status == 0 .and. flux_balance < 3.0e-2_dp, &
        "Costabel coupling resolves the analytical point-source flux")
    call record_condition(boundary_error < 3.0e-2_dp, &
        "Costabel coupling matches the analytical sphere potential")
    call check_summary("Three-dimensional symmetric Costabel FEM-BEM coupling")
    if (.not. all_passed) error stop 1

contains

    pure function triangle_areas(mesh_vertices, cells) result(areas)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp) :: areas(size(cells, 2))

        real(dp) :: first(3), second(3)
        integer :: triangle

        do triangle = 1, size(cells, 2)
            first = mesh_vertices(:, cells(2, triangle)) - &
                mesh_vertices(:, cells(1, triangle))
            second = mesh_vertices(:, cells(3, triangle)) - &
                mesh_vertices(:, cells(1, triangle))
            areas(triangle) = 0.5_dp*norm2([ &
                first(2)*second(3) - first(3)*second(2), &
                first(3)*second(1) - first(1)*second(3), &
                first(1)*second(2) - first(2)*second(1)])
        end do
    end function triangle_areas

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_costabel_coupling_3d
