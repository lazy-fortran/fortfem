program test_scalar_helmholtz_planar_dtn_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_scalar_helmholtz_planar_dtn_p1, &
        solve_scalar_helmholtz_planar_dtn_p1_jvp, &
        solve_scalar_helmholtz_planar_dtn_p1_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: side_points = 4
    real(dp), parameter :: difference_step = 1.0e-6_dp
    complex(dp), allocatable :: dirichlet_bar(:), dirichlet_dot(:)
    complex(dp), allocatable :: dirichlet_values(:), load(:), load_bar(:)
    complex(dp), allocatable :: load_dot(:), solution(:), solution_bar(:)
    complex(dp), allocatable :: solution_dot(:), solution_minus(:)
    complex(dp), allocatable :: solution_plus(:)
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: dirichlet_nodes(:), dtn_nodes(:), triangles(:, :)
    real(dp) :: forward_pairing, period, period_bar, period_dot
    real(dp) :: reverse_pairing, wavenumber, wavenumber_bar, wavenumber_dot
    integer :: dirichlet_count, entry, node, status, status_minus, status_plus
    logical :: all_passed

    all_passed = .true.
    allocate ( &
        vertices(2, side_points**2), &
        triangles(3, 2*(side_points - 1)**2))
    call build_square_mesh(vertices, triangles)
    allocate (dtn_nodes(side_points))
    do node = 1, side_points
        dtn_nodes(node) = node*side_points
    end do
    dirichlet_count = 3*side_points - 2
    allocate (dirichlet_nodes(dirichlet_count))
    call build_dirichlet_nodes(dirichlet_nodes)
    allocate ( &
        dirichlet_values(dirichlet_count), dirichlet_dot(dirichlet_count), &
        dirichlet_bar(dirichlet_count))
    do entry = 1, dirichlet_count
        dirichlet_values(entry) = cmplx( &
            0.2_dp*sin(0.3_dp*real(entry, dp)), &
            -0.1_dp*cos(0.2_dp*real(entry, dp)), dp)
        dirichlet_dot(entry) = cmplx( &
            0.04_dp*cos(0.4_dp*real(entry, dp)), &
            0.03_dp*sin(0.5_dp*real(entry, dp)), dp)
    end do
    allocate ( &
        load(side_points**2), load_dot(side_points**2), &
        load_bar(side_points**2), solution(side_points**2), &
        solution_dot(side_points**2), solution_minus(side_points**2), &
        solution_plus(side_points**2), solution_bar(side_points**2))
    do node = 1, side_points**2
        load(node) = cmplx( &
            0.05_dp*sin(0.7_dp*real(node, dp)), &
            0.03_dp*cos(0.6_dp*real(node, dp)), dp)
        load_dot(node) = cmplx( &
            -0.02_dp*cos(0.2_dp*real(node, dp)), &
            0.01_dp*sin(0.3_dp*real(node, dp)), dp)
        solution_bar(node) = cmplx( &
            0.06_dp*cos(0.5_dp*real(node, dp)), &
            -0.04_dp*sin(0.4_dp*real(node, dp)), dp)
    end do
    wavenumber = 2.3_dp
    period = real(side_points, dp)/real(side_points - 1, dp)
    wavenumber_dot = 0.07_dp
    period_dot = -0.03_dp

    call solve_scalar_helmholtz_planar_dtn_p1( &
        vertices, triangles, dtn_nodes, wavenumber, period, load, &
        dirichlet_nodes, dirichlet_values, solution, status)
    call solve_scalar_helmholtz_planar_dtn_p1_jvp( &
        vertices, triangles, dtn_nodes, wavenumber, period, load, &
        dirichlet_nodes, dirichlet_values, wavenumber_dot, period_dot, &
        load_dot, dirichlet_dot, solution_dot, status)
    call solve_scalar_helmholtz_planar_dtn_p1( &
        vertices, triangles, dtn_nodes, &
        wavenumber - difference_step*wavenumber_dot, &
        period - difference_step*period_dot, &
        load - difference_step*load_dot, dirichlet_nodes, &
        dirichlet_values - difference_step*dirichlet_dot, &
        solution_minus, status_minus)
    call solve_scalar_helmholtz_planar_dtn_p1( &
        vertices, triangles, dtn_nodes, &
        wavenumber + difference_step*wavenumber_dot, &
        period + difference_step*period_dot, &
        load + difference_step*load_dot, dirichlet_nodes, &
        dirichlet_values + difference_step*dirichlet_dot, &
        solution_plus, status_plus)
    call record_condition(status == 0 .and. status_minus == 0 .and. &
        status_plus == 0, "FEM-DtN state JVP accepts valid inputs")
    call record_condition(maxval(abs(solution_dot - &
        (solution_plus - solution_minus)/(2.0_dp*difference_step))) < 2.0e-8_dp, &
        "FEM-DtN state JVP matches complete reassemble/refactor difference")

    call solve_scalar_helmholtz_planar_dtn_p1_vjp( &
        vertices, triangles, dtn_nodes, wavenumber, period, load, &
        dirichlet_nodes, dirichlet_values, solution, solution_bar, &
        load_bar, dirichlet_bar, wavenumber_bar, period_bar, status)
    forward_pairing = real(sum(conjg(solution_bar)*solution_dot), dp)
    reverse_pairing = real(sum(conjg(load_bar)*load_dot) + &
        sum(conjg(dirichlet_bar)*dirichlet_dot), dp) + &
        wavenumber_bar*wavenumber_dot + period_bar*period_dot
    call record_condition(status == 0, "FEM-DtN state VJP succeeds")
    call record_condition(abs(forward_pairing - reverse_pairing) < 2.0e-11_dp, &
        "FEM-DtN state products satisfy the complete adjoint identity")

    call check_summary("Scalar Helmholtz planar FEM-DtN state derivatives")
    if (.not. all_passed) error stop 1

contains

    subroutine build_square_mesh(mesh_vertices, mesh_triangles)
        real(dp), intent(out) :: mesh_vertices(:, :)
        integer, intent(out) :: mesh_triangles(:, :)

        real(dp) :: spacing
        integer :: column, lower_left, row, triangle, vertex

        spacing = 1.0_dp/real(side_points - 1, dp)
        vertex = 0
        do row = 0, side_points - 1
            do column = 0, side_points - 1
                vertex = vertex + 1
                mesh_vertices(:, vertex) = spacing*real([column, row], dp)
            end do
        end do
        triangle = 0
        do row = 1, side_points - 1
            do column = 1, side_points - 1
                lower_left = column + (row - 1)*side_points
                triangle = triangle + 1
                mesh_triangles(:, triangle) = [ &
                    lower_left, lower_left + 1, &
                    lower_left + side_points + 1]
                triangle = triangle + 1
                mesh_triangles(:, triangle) = [ &
                    lower_left, lower_left + side_points + 1, &
                    lower_left + side_points]
            end do
        end do
    end subroutine build_square_mesh

    subroutine build_dirichlet_nodes(nodes)
        integer, intent(out) :: nodes(:)

        integer :: column, entry, row

        entry = 0
        do row = 1, side_points
            entry = entry + 1
            nodes(entry) = (row - 1)*side_points + 1
        end do
        do column = 2, side_points
            entry = entry + 1
            nodes(entry) = column
        end do
        do column = 2, side_points
            entry = entry + 1
            nodes(entry) = (side_points - 1)*side_points + column
        end do
    end subroutine build_dirichlet_nodes

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_scalar_helmholtz_planar_dtn_ad
