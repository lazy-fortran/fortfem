program test_laplace_symmetric_coupling_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_laplace_symmetric_coupling_p1_p0
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: refinement_count = 3
    integer, parameter :: points_per_side(refinement_count) = [3, 5, 9]
    real(dp) :: flux_error(refinement_count)
    real(dp) :: interior_error(refinement_count)
    integer :: refinement
    logical :: all_passed

    all_passed = .true.
    do refinement = 1, refinement_count
        call dipole_transmission_error( &
            points_per_side(refinement), interior_error(refinement), &
            flux_error(refinement))
    end do

    call record_condition(interior_error(2) < 0.7_dp * interior_error(1) .and. &
        interior_error(3) < 0.7_dp * interior_error(2), &
        "Dipole interior error decreases under uniform refinement")
    call record_condition(flux_error(2) < 0.7_dp * flux_error(1) .and. &
        flux_error(3) < 0.7_dp * flux_error(2), &
        "Dipole exterior flux error decreases under uniform refinement")
    call record_condition(interior_error(3) < 5.0e-3_dp .and. &
        flux_error(3) < 2.0e-2_dp, &
        "Refined dipole transmission solution reaches its accuracy target")

    call check_summary("Laplace symmetric coupling convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine dipole_transmission_error( &
            side_points, interior_error, flux_error)
        integer, intent(in) :: side_points
        real(dp), intent(out) :: interior_error, flux_error

        real(dp), allocatable :: dirichlet_jump(:), exact_flux(:)
        real(dp), allocatable :: exterior_flux(:), interior_solution(:)
        real(dp), allocatable :: neumann_jump(:), panel_end(:, :)
        real(dp), allocatable :: panel_start(:, :), vertices(:, :)
        real(dp), allocatable :: volume_load(:)
        integer, allocatable :: panel_nodes(:, :), triangles(:, :)
        real(dp) :: exact_flux_norm, perimeter
        integer :: endpoint, panel, panel_count, status, triangle_count
        integer :: vertex_count

        vertex_count = side_points**2
        triangle_count = 2 * (side_points - 1)**2
        panel_count = 4 * (side_points - 1)
        allocate(vertices(2, vertex_count), triangles(3, triangle_count))
        allocate(panel_nodes(2, panel_count))
        allocate(panel_start(2, panel_count), panel_end(2, panel_count))
        call build_square_mesh(side_points, vertices, triangles)
        call build_square_boundary( &
            side_points, vertices, panel_nodes, panel_start, panel_end)

        allocate(dirichlet_jump(vertex_count), volume_load(vertex_count))
        allocate(interior_solution(vertex_count))
        allocate(neumann_jump(panel_count), exact_flux(panel_count))
        allocate(exterior_flux(panel_count))
        dirichlet_jump = 0.0_dp
        do panel = 1, panel_count
            do endpoint = 1, 2
                dirichlet_jump(panel_nodes(endpoint, panel)) = &
                    -vertices(1, panel_nodes(endpoint, panel)) / &
                    sum(vertices(:, panel_nodes(endpoint, panel))**2)
            end do
        end do
        volume_load = 0.0_dp
        call project_dipole_flux( &
            panel_start, panel_end, exact_flux, perimeter)
        neumann_jump = -exact_flux

        call solve_laplace_symmetric_coupling_p1_p0( &
            vertices, triangles, panel_start, panel_end, panel_nodes, 20, &
            volume_load, dirichlet_jump, neumann_jump, interior_solution, &
            exterior_flux, status)
        call record_condition(status == 0, "Dipole transmission solve succeeds")
        interior_error = sqrt(sum(interior_solution**2) / real(vertex_count, dp))
        exact_flux_norm = sqrt(sum( &
            panel_lengths(panel_start, panel_end) * exact_flux**2) / perimeter)
        flux_error = sqrt(sum( &
            panel_lengths(panel_start, panel_end) * &
            (exterior_flux - exact_flux)**2) / perimeter) / exact_flux_norm
    end subroutine dipole_transmission_error

    subroutine project_dipole_flux( &
            panel_start, panel_end, flux, perimeter)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp), intent(out) :: flux(:), perimeter

        integer, parameter :: quadrature_order = 32
        real(dp) :: gradient(2), length, nodes(quadrature_order)
        real(dp) :: normal(2), point(2), weights(quadrature_order)
        integer :: node, panel

        call gauss_legendre_ab( &
            quadrature_order, 0.0_dp, 1.0_dp, nodes, weights)
        flux = 0.0_dp
        perimeter = 0.0_dp
        do panel = 1, size(panel_start, 2)
            length = norm2(panel_end(:, panel) - panel_start(:, panel))
            perimeter = perimeter + length
            normal = [panel_end(2, panel) - panel_start(2, panel), &
                panel_start(1, panel) - panel_end(1, panel)] / length
            do node = 1, quadrature_order
                point = panel_start(:, panel) + nodes(node) * &
                    (panel_end(:, panel) - panel_start(:, panel))
                gradient = dipole_gradient(point)
                flux(panel) = flux(panel) + &
                    weights(node) * dot_product(gradient, normal)
            end do
        end do
    end subroutine project_dipole_flux

    pure function dipole_gradient(point) result(gradient)
        real(dp), intent(in) :: point(2)
        real(dp) :: gradient(2)
        real(dp) :: radius_squared

        radius_squared = sum(point**2)
        gradient(1) = (point(2)**2 - point(1)**2) / radius_squared**2
        gradient(2) = -2.0_dp * point(1) * point(2) / radius_squared**2
    end function dipole_gradient

    pure function panel_lengths(panel_start, panel_end) result(lengths)
        real(dp), intent(in) :: panel_start(:, :), panel_end(:, :)
        real(dp) :: lengths(size(panel_start, 2))
        integer :: panel

        do panel = 1, size(panel_start, 2)
            lengths(panel) = &
                norm2(panel_end(:, panel) - panel_start(:, panel))
        end do
    end function panel_lengths

    subroutine build_square_mesh(side_points, vertices, triangles)
        integer, intent(in) :: side_points
        real(dp), intent(out) :: vertices(:, :)
        integer, intent(out) :: triangles(:, :)

        real(dp), parameter :: half_width = 0.1_dp
        real(dp) :: spacing
        integer :: column, lower_left, row, triangle, vertex

        spacing = 2.0_dp * half_width / real(side_points - 1, dp)
        vertex = 0
        do row = 0, side_points - 1
            do column = 0, side_points - 1
                vertex = vertex + 1
                vertices(:, vertex) = -half_width + &
                    spacing * real([column, row], dp)
            end do
        end do

        triangle = 0
        do row = 1, side_points - 1
            do column = 1, side_points - 1
                lower_left = column + (row - 1) * side_points
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + 1, &
                    lower_left + side_points + 1]
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + side_points + 1, &
                    lower_left + side_points]
            end do
        end do
    end subroutine build_square_mesh

    subroutine build_square_boundary( &
            side_points, vertices, panel_nodes, panel_start, panel_end)
        integer, intent(in) :: side_points
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(out) :: panel_nodes(:, :)
        real(dp), intent(out) :: panel_start(:, :), panel_end(:, :)

        integer :: boundary_nodes(size(panel_nodes, 2))
        integer :: index, panel, panel_count

        panel_count = size(panel_nodes, 2)
        index = 0
        do panel = 1, side_points
            index = index + 1
            boundary_nodes(index) = panel
        end do
        do panel = 2, side_points
            index = index + 1
            boundary_nodes(index) = panel * side_points
        end do
        do panel = side_points - 1, 1, -1
            index = index + 1
            boundary_nodes(index) = (side_points - 1) * side_points + panel
        end do
        do panel = side_points - 1, 2, -1
            index = index + 1
            boundary_nodes(index) = (panel - 1) * side_points + 1
        end do

        do panel = 1, panel_count
            panel_nodes(1, panel) = boundary_nodes(panel)
            panel_nodes(2, panel) = boundary_nodes(mod(panel, panel_count) + 1)
            panel_start(:, panel) = vertices(:, panel_nodes(1, panel))
            panel_end(:, panel) = vertices(:, panel_nodes(2, panel))
        end do
    end subroutine build_square_boundary

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_symmetric_coupling_convergence
