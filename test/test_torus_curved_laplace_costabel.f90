program test_torus_curved_laplace_costabel
    use check, only: check_condition, check_summary
    use fortfem_core, only: evaluate_torus_curved_panel, &
        generate_solid_torus_tetra_mesh
    use fortfem_boundary, only: &
        solve_laplace_fem_bem_costabel_torus_curved_3d, &
        triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    integer, parameter :: counts(2) = [5, 7]
    integer, allocatable :: boundary_triangles(:, :), tetrahedra(:, :)
    real(dp), allocatable :: eta(:), exact(:), flux(:), load(:), parameters(:, :)
    real(dp), allocatable :: potential(:), vertices(:, :), weights(:), xi(:)
    logical, allocatable :: is_boundary(:)
    real(dp) :: errors(2), flux_errors(2), flux_integral, jacobian
    real(dp) :: point(3), source(3)
    real(dp) :: tangent_eta(3), tangent_xi(3)
    integer :: level, local_node, node, panel, quadrature_status, status
    logical :: all_passed

    all_passed = .true.
    call triangle_duffy_quadrature(7, xi, eta, weights, quadrature_status)
    if (quadrature_status /= 0) error stop "torus area quadrature failed"
    do level = 1, 2
        call generate_solid_torus_tetra_mesh( &
            major_radius, minor_radius, 2, counts(level), counts(level) + 2, &
            vertices, tetrahedra, boundary_triangles, parameters)
        allocate(load(size(vertices, 2)))
        load = 0.0_dp
        load(1) = 1.0_dp
        source = vertices(:, 1)
        call solve_laplace_fem_bem_costabel_torus_curved_3d( &
            vertices, tetrahedra, parameters, boundary_triangles, &
            major_radius, minor_radius, load, 5, potential, flux, status)
        if (status /= 0) error stop "curved torus Costabel solve failed"
        allocate(is_boundary(size(vertices, 2)), exact(size(vertices, 2)))
        is_boundary = .false.
        do panel = 1, size(boundary_triangles, 2)
            do local_node = 1, 3
                is_boundary(boundary_triangles(local_node, panel)) = .true.
            end do
        end do
        exact = 0.0_dp
        do node = 1, size(vertices, 2)
            if (is_boundary(node)) then
                exact(node) = 1.0_dp/( &
                    4.0_dp*acos(-1.0_dp)*norm2(vertices(:, node) - source))
            end if
        end do
        errors(level) = norm2(pack(potential - exact, is_boundary))/ &
            norm2(pack(exact, is_boundary))
        flux_integral = 0.0_dp
        do panel = 1, size(boundary_triangles, 2)
            do node = 1, size(weights)
                call evaluate_torus_curved_panel( &
                    parameters(:, boundary_triangles(:, panel)), major_radius, &
                    minor_radius, xi(node), eta(node), point, tangent_xi, &
                    tangent_eta, jacobian, status)
                flux_integral = flux_integral + &
                    weights(node)*jacobian*flux(panel)
            end do
        end do
        write (*, '(a,es14.6)') "integrated exterior flux: ", flux_integral
        flux_errors(level) = abs(flux_integral + 1.0_dp)
        deallocate( &
            vertices, tetrahedra, boundary_triangles, parameters, load, &
            potential, flux, is_boundary, exact)
    end do

    write (*, '(a,2(es12.4,1x))') "boundary potential errors: ", errors
    call record_condition(flux_errors(2) < flux_errors(1), &
        "curved torus Costabel flux converges to point-source conservation")
    call record_condition(flux_errors(2) < 3.0e-2_dp, &
        "curved torus Costabel flux reaches the curved-geometry tolerance")
    call record_condition(errors(2) < 0.85_dp*errors(1), &
        "curved torus Costabel potential converges under refinement")
    call record_condition(errors(2) < 3.5e-1_dp, &
        "curved torus Costabel potential matches the point-source solution")
    call check_summary("Curved torus Laplace Costabel coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_costabel
