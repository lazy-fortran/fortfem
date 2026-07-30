program test_laplace_bem_galerkin_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_laplace_single_layer_p0_hierarchical_3d, &
        assemble_laplace_calderon_p1_p0_3d, &
        evaluate_helmholtz_representation_triangles_3d, &
        evaluate_helmholtz_cfie_p0_3d, solve_helmholtz_cfie_p0_3d, &
        solve_helmholtz_dirichlet_p0_3d, &
        solve_laplace_fem_bem_johnson_nedelec_3d, &
        solve_laplace_dirichlet_p0_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), allocatable :: adjoint(:, :), density(:)
    real(dp), allocatable :: double_layer(:, :), hypersingular(:, :)
    real(dp), allocatable :: ones(:), single_layer(:, :), vertices(:, :)
    real(dp), allocatable :: fem_solution(:), flux(:), load(:)
    real(dp), allocatable :: dense_action(:), fast_action(:)
    real(dp), allocatable :: volume_vertices(:, :)
    complex(dp), allocatable :: complex_density(:), dirichlet(:)
    integer, allocatable :: triangles(:, :)
    integer, allocatable :: tetrahedra(:, :)
    complex(dp) :: exact_field, numerical_field
    real(dp) :: capacities(0:2), errors(0:2)
    integer :: interaction_count, level, status
    logical :: all_passed

    all_passed = .true.
    do level = 0, 2
        call build_sphere(level, vertices, triangles)
        call solve_laplace_dirichlet_p0_3d( &
            vertices, triangles, 1.0_dp, 8, density, capacities(level), status)
        errors(level) = abs(capacities(level) - 4.0_dp*acos(-1.0_dp))
        call record_condition(status == 0 .and. all(density > 0.0_dp), &
            "Three-dimensional Galerkin BEM gives positive sphere charge")
    end do
    call record_condition(errors(1) < 0.65_dp*errors(0) .and. &
        errors(2) < 0.65_dp*errors(1), &
        "Three-dimensional Galerkin BEM converges to sphere capacitance")
    if (errors(2) >= 0.8_dp) then
        write (*, '(A,3ES14.5,A,3ES14.5)') &
            "sphere capacities ", capacities, " errors ", errors
    end if
    call record_condition(errors(2) < 0.8_dp, &
        "Refined sphere capacitance matches the analytical value")
    call assemble_laplace_calderon_p1_p0_3d( &
        vertices, triangles, 8, single_layer, double_layer, adjoint, &
        hypersingular, status)
    allocate(ones(size(vertices, 2)))
    ones = 1.0_dp
    call record_condition(status == 0 .and. &
        maxval(abs(adjoint - transpose(double_layer))) < 2.0e-13_dp, &
        "Three-dimensional adjoint double layer is the Galerkin transpose")
    call record_condition(maxval(abs(hypersingular - &
        transpose(hypersingular))) < 2.0e-12_dp .and. &
        maxval(abs(matmul(hypersingular, ones))) < 2.0e-11_dp, &
        "Regularized 3D hypersingular operator is symmetric and kills constants")
    call record_condition(maxval(abs(matmul(double_layer, ones) + &
        0.5_dp*triangle_areas(vertices, triangles))) < 2.0e-2_dp, &
        "Three-dimensional double layer has the closed-surface constant trace")
    dense_action = matmul(single_layer, density)
    call apply_laplace_single_layer_p0_hierarchical_3d( &
        vertices, triangles, density, 0.6_dp, 6, fast_action, status, &
        interaction_count)
    call record_condition(status == 0 .and. &
        norm2(fast_action - dense_action)/norm2(dense_action) < 3.0e-2_dp, &
        "Hierarchical 3D single layer agrees with dense Galerkin action")
    if (interaction_count >= size(triangles, 2)**2/2) then
        write (*, '(A,I0,A,I0)') "hierarchical interactions ", &
            interaction_count, " dense ", size(triangles, 2)**2
    end if
    call record_condition(interaction_count < size(triangles, 2)**2/2, &
        "Hierarchical 3D single layer uses subquadratic interactions")

    allocate(volume_vertices(3, size(vertices, 2) + 1))
    volume_vertices(:, :size(vertices, 2)) = vertices
    volume_vertices(:, size(vertices, 2) + 1) = 0.0_dp
    allocate(tetrahedra(4, size(triangles, 2)))
    tetrahedra(1, :) = size(vertices, 2) + 1
    tetrahedra(2:4, :) = triangles
    allocate(load(size(volume_vertices, 2)))
    load = 0.0_dp
    load(size(load)) = 1.0_dp
    call solve_laplace_fem_bem_johnson_nedelec_3d( &
        volume_vertices, tetrahedra, triangles, load, 8, fem_solution, &
        flux, status)
    call record_condition(status == 0 .and. abs(sum( &
        flux*triangle_areas(vertices, triangles)) + 1.0_dp) < 2.0e-10_dp, &
        "Three-dimensional FEM-BEM conserves unit point-source flux")
    call record_condition(abs(sum(fem_solution(:size(vertices, 2)))/ &
        real(size(vertices, 2), dp) - 1.0_dp/(4.0_dp*acos(-1.0_dp))) < &
        2.5e-2_dp, &
        "Three-dimensional FEM-BEM matches the point-charge boundary potential")

    call solve_helmholtz_dirichlet_p0_3d( &
        vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), 0.7_dp, 8, &
        complex_density, status)
    allocate(dirichlet(size(vertices, 2)))
    dirichlet = cmplx(0.0_dp, 0.0_dp, dp)
    call evaluate_helmholtz_representation_triangles_3d( &
        vertices, triangles, dirichlet, -complex_density, &
        [0.0_dp, 0.0_dp, 2.0_dp], 0.7_dp, 8, numerical_field, status)
    exact_field = 0.5_dp*exp(cmplx(0.0_dp, 0.7_dp, dp))
    call record_condition(status == 0 .and. &
        abs(numerical_field - exact_field) < 6.0e-2_dp, &
        "Three-dimensional Helmholtz BEM matches an outgoing sphere mode")

    call solve_helmholtz_cfie_p0_3d( &
        vertices, triangles, cmplx(1.0_dp, 0.0_dp, dp), acos(-1.0_dp), &
        acos(-1.0_dp), 8, complex_density, status)
    call evaluate_helmholtz_cfie_p0_3d( &
        vertices, triangles, complex_density, [0.0_dp, 0.0_dp, 2.0_dp], &
        acos(-1.0_dp), acos(-1.0_dp), 8, numerical_field, status)
    call record_condition(status == 0 .and. &
        abs(numerical_field + 0.5_dp) < 1.2e-1_dp, &
        "Three-dimensional CFIE remains accurate at an interior resonance")

    call solve_laplace_dirichlet_p0_3d( &
        vertices(:, :2), triangles, 1.0_dp, 8, density, capacities(2), status)
    call record_condition(status /= 0, &
        "Three-dimensional Galerkin BEM rejects an invalid surface")

    call check_summary("Three-dimensional Laplace Galerkin BEM")
    if (.not. all_passed) error stop 1

contains

    subroutine build_sphere(level, vertices, triangles)
        integer, intent(in) :: level
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)

        integer :: refinement

        allocate(vertices(3, 6), triangles(3, 8))
        vertices(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 2) = [-1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        vertices(:, 4) = [0.0_dp, -1.0_dp, 0.0_dp]
        vertices(:, 5) = [0.0_dp, 0.0_dp, 1.0_dp]
        vertices(:, 6) = [0.0_dp, 0.0_dp, -1.0_dp]
        triangles(:, 1) = [1, 3, 5]
        triangles(:, 2) = [3, 2, 5]
        triangles(:, 3) = [2, 4, 5]
        triangles(:, 4) = [4, 1, 5]
        triangles(:, 5) = [3, 1, 6]
        triangles(:, 6) = [2, 3, 6]
        triangles(:, 7) = [4, 2, 6]
        triangles(:, 8) = [1, 4, 6]
        do refinement = 1, level
            call refine_sphere(vertices, triangles)
        end do
    end subroutine build_sphere

    subroutine refine_sphere(vertices, triangles)
        real(dp), allocatable, intent(inout) :: vertices(:, :)
        integer, allocatable, intent(inout) :: triangles(:, :)

        real(dp), allocatable :: expanded_vertices(:, :)
        integer, allocatable :: edge_midpoints(:), edge_vertices(:, :)
        integer, allocatable :: refined(:, :)
        integer :: edge_count, midpoint(3), old_vertex_count
        integer :: triangle, vertex_count

        old_vertex_count = size(vertices, 2)
        allocate(expanded_vertices(3, old_vertex_count + 3*size(triangles, 2)))
        expanded_vertices(:, :old_vertex_count) = vertices
        allocate(edge_vertices(2, 3*size(triangles, 2)))
        allocate(edge_midpoints(3*size(triangles, 2)))
        allocate(refined(3, 4*size(triangles, 2)))
        vertex_count = old_vertex_count
        edge_count = 0
        do triangle = 1, size(triangles, 2)
            midpoint(1) = midpoint_vertex( &
                triangles(1, triangle), triangles(2, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            midpoint(2) = midpoint_vertex( &
                triangles(2, triangle), triangles(3, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            midpoint(3) = midpoint_vertex( &
                triangles(3, triangle), triangles(1, triangle), vertices, &
                expanded_vertices, edge_vertices, edge_midpoints, &
                edge_count, vertex_count)
            refined(:, 4*triangle - 3) = &
                [triangles(1, triangle), midpoint(1), midpoint(3)]
            refined(:, 4*triangle - 2) = &
                [midpoint(1), triangles(2, triangle), midpoint(2)]
            refined(:, 4*triangle - 1) = &
                [midpoint(3), midpoint(2), triangles(3, triangle)]
            refined(:, 4*triangle) = midpoint
        end do
        vertices = expanded_vertices(:, :vertex_count)
        call move_alloc(refined, triangles)
    end subroutine refine_sphere

    function midpoint_vertex( &
            first, second, vertices, expanded_vertices, edge_vertices, &
            edge_midpoints, edge_count, vertex_count) result(midpoint)
        integer, intent(in) :: first, second
        real(dp), intent(in) :: vertices(:, :)
        real(dp), intent(inout) :: expanded_vertices(:, :)
        integer, intent(inout) :: edge_vertices(:, :), edge_midpoints(:)
        integer, intent(inout) :: edge_count, vertex_count
        integer :: midpoint

        integer :: edge, ordered(2)
        real(dp) :: midpoint_coordinate(3)

        ordered = [min(first, second), max(first, second)]
        do edge = 1, edge_count
            if (all(edge_vertices(:, edge) == ordered)) then
                midpoint = edge_midpoints(edge)
                return
            end if
        end do
        edge_count = edge_count + 1
        edge_vertices(:, edge_count) = ordered
        vertex_count = vertex_count + 1
        midpoint_coordinate = vertices(:, first) + vertices(:, second)
        expanded_vertices(:, vertex_count) = &
            midpoint_coordinate/norm2(midpoint_coordinate)
        edge_midpoints(edge_count) = vertex_count
        midpoint = vertex_count
    end function midpoint_vertex

    pure function triangle_areas(vertices, triangles) result(areas)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp) :: areas(size(triangles, 2))

        real(dp) :: first(3), second(3)
        integer :: triangle

        do triangle = 1, size(triangles, 2)
            first = vertices(:, triangles(2, triangle)) - &
                vertices(:, triangles(1, triangle))
            second = vertices(:, triangles(3, triangle)) - &
                vertices(:, triangles(1, triangle))
            areas(triangle) = 0.5_dp*norm2([ &
                first(2)*second(3) - first(3)*second(2), &
                first(3)*second(1) - first(1)*second(3), &
                first(1)*second(2) - first(2)*second(1)])
        end do
    end function triangle_areas

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_bem_galerkin_3d
