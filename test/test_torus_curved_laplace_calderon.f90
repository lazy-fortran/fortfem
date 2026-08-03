program test_torus_curved_laplace_calderon
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        assemble_laplace_torus_curved_calderon_3d, &
        solve_laplace_bem_dtn_torus_curved_3d
    use fortfem_core, only: cartesian_to_toroidal, generate_torus_surface_mesh
    use fortfem_fourier, only: evaluate_toroidal_harmonic_p
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp, minor_radius = 0.6_dp
    real(dp), parameter :: scale = &
        sqrt(major_radius**2 - minor_radius**2)
    real(dp), parameter :: boundary_eta = acosh(major_radius/minor_radius)
    integer, parameter :: counts(2) = [5, 7]
    integer, allocatable :: triangles(:, :)
    real(dp), allocatable :: adjoint(:, :), double_layer(:, :)
    real(dp), allocatable :: flux(:), hypersingular(:, :), mass(:, :)
    real(dp), allocatable :: parameters(:, :), single_layer(:, :), trace(:)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: eta, first_term_norm, phi, residuals(2), second_term_norm
    real(dp) :: theta
    integer :: level, status, vertex
    logical :: all_passed

    all_passed = .true.
    do level = 1, 2
        call generate_torus_surface_mesh( &
            major_radius, minor_radius, counts(level), counts(level) + 2, &
            vertices, triangles, parameters)
        allocate(trace(size(vertices, 2)))
        do vertex = 1, size(vertices, 2)
            call cartesian_to_toroidal( &
                vertices(:, vertex), scale, eta, theta, phi)
            call evaluate_toroidal_harmonic_p( &
                2, 1, boundary_eta, theta, phi, trace(vertex), status)
            if (status /= 0) error stop "toroidal Calderon trace failed"
        end do
        call solve_laplace_bem_dtn_torus_curved_3d( &
            parameters, triangles, major_radius, minor_radius, trace, 5, &
            flux, status)
        if (status /= 0) error stop "toroidal Calderon flux solve failed"
        call assemble_laplace_torus_curved_calderon_3d( &
            parameters, triangles, major_radius, minor_radius, 5, &
            single_layer, double_layer, adjoint, hypersingular, mass, status)
        if (status /= 0) error stop "curved torus Calderon assembly failed"
        first_term_norm = norm2(matmul(hypersingular, trace))
        second_term_norm = norm2(matmul( &
            adjoint + 0.5_dp*transpose(mass), flux))
        residuals(level) = norm2( &
            matmul(hypersingular, trace) + &
            matmul(adjoint + 0.5_dp*transpose(mass), flux))/ &
            max(first_term_norm + second_term_norm, tiny(1.0_dp))
        if (level == 2) then
            call record_condition(maxval(abs(adjoint - &
                transpose(double_layer))) < 1.0e-13_dp, &
                "curved torus adjoint double layer is the Galerkin transpose")
            call record_condition(maxval(abs(hypersingular - &
                transpose(hypersingular))) < 2.0e-12_dp, &
                "curved torus hypersingular operator is symmetric")
            call record_condition(norm2(matmul( &
                hypersingular, spread(1.0_dp, 1, size(trace)))) < 2.0e-11_dp, &
                "curved torus hypersingular operator kills constants")
        end if
        deallocate( &
            vertices, triangles, parameters, trace, flux, single_layer, &
            double_layer, adjoint, hypersingular, mass)
    end do

    write (*, '(a,2(es12.4,1x))') "second Calderon residuals: ", residuals
    call record_condition(residuals(2) < 0.8_dp*residuals(1), &
        "curved torus second Calderon equation converges under refinement")
    call record_condition(residuals(2) < 3.0e-1_dp, &
        "curved torus second Calderon equation matches the harmonic trace")
    call check_summary("Curved torus Laplace Calderon")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_torus_curved_laplace_calderon
