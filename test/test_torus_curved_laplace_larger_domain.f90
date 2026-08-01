program test_torus_curved_laplace_larger_domain
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        evaluate_laplace_representation_torus_curved_3d, &
        generate_torus_surface_mesh, solve_laplace_bem_dtn_torus_curved_3d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: minor_radius = 0.6_dp
    real(dp), parameter :: major_radii(2) = [2.0_dp, 2.5_dp]
    integer, parameter :: mesh_count = 8
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: parameters(:, :), trace(:)
    real(dp), allocatable :: vertices(:, :), dtn_flux(:)
    real(dp) :: dtn_values(2), exact, radius
    real(dp) :: displacement(3), point(3), source(3)
    real(dp) :: dtn_errors(2)
    integer :: surface, status, vertex
    logical :: all_passed

    all_passed = .true.
    point = [0.30_dp, 0.15_dp, 0.0_dp]
    source = [2.25_dp, 0.0_dp, 0.30_dp]
    displacement = point - source
    radius = norm2(displacement)
    exact = 1.0_dp/(4.0_dp*pi*radius)
    do surface = 1, 2
        call generate_torus_surface_mesh( &
            major_radii(surface), minor_radius, mesh_count, mesh_count + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)))
        do vertex = 1, size(vertices, 2)
            displacement = vertices(:, vertex) - source
            radius = norm2(displacement)
            trace(vertex) = 1.0_dp/(4.0_dp*pi*radius)
        end do
        call solve_laplace_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radii(surface), minor_radius, &
            trace, 5, dtn_flux, status)
        if (status /= 0) error stop "larger-domain torus DtN solve failed"
        call evaluate_laplace_representation_torus_curved_3d( &
            parameters, triangles, major_radii(surface), minor_radius, &
            trace, dtn_flux, point, 5, dtn_values(surface), status)
        if (status /= 0) error stop "larger-domain DtN reconstruction failed"
        dtn_errors(surface) = abs(dtn_values(surface) - exact)/abs(exact)
        deallocate(vertices, triangles, parameters, trace, dtn_flux)
    end do

    write (*, '(a,2(es12.4,1x))') "BEM/DtN reconstruction errors: ", dtn_errors
    write (*, '(a,2(es12.4,1x))') "BEM/DtN fields: ", dtn_values
    write (*, '(a,es12.4)') "exact field: ", exact
    call record_condition(maxval(dtn_errors) < 2.5e-1_dp, &
        "the toroidal BEM/DtN path matches the fundamental field")
    call record_condition(abs(dtn_values(2) - dtn_values(1))/abs(exact) < 1.5e-1_dp, &
        "moving the toroidal artificial boundary preserves the interior field")
    call check_summary("Torus curved Laplace larger-domain control")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_larger_domain
