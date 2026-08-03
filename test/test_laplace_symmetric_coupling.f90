program test_laplace_symmetric_coupling
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: assemble_laplace_double_layer_mixed_linear, &
        assemble_laplace_hypersingular_linear, &
        assemble_laplace_single_layer_constant, &
        assemble_laplace_symmetric_coupling_p1_p0, &
        solve_laplace_symmetric_coupling_p1_p0
    use fortfem_kinds, only: dp
    implicit none

    real(dp) :: boundary_load(3), coupling(6, 6), double_layer(3, 3)
    real(dp) :: expected(6), flux(3), hypersingular(3, 3), mass(3, 3)
    real(dp) :: normal(2), panel_end(2, 3), panel_start(2, 3)
    real(dp) :: single_layer(3, 3), state(6), trace(3), vertices(2, 3)
    real(dp) :: gradient(2), length, expected_energy, observed_energy
    real(dp) :: interior_solution(3), exterior_flux(3), neumann_jump(3)
    real(dp) :: volume_load(3)
    integer :: endpoint, panel, panel_nodes(2, 3), status, triangles(3, 1)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp]
    vertices(:, 2) = [0.25_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 0.25_dp]
    triangles(:, 1) = [1, 2, 3]
    panel_nodes(:, 1) = [1, 2]
    panel_nodes(:, 2) = [2, 3]
    panel_nodes(:, 3) = [3, 1]
    do panel = 1, 3
        panel_start(:, panel) = vertices(:, panel_nodes(1, panel))
        panel_end(:, panel) = vertices(:, panel_nodes(2, panel))
    end do

    call assemble_laplace_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, 24, &
        coupling, status)
    call record_condition(status == 0, "Costabel coupling assembles")
    call record_condition(maxval(abs( &
        coupling(1:3, 1:3) - transpose(coupling(1:3, 1:3)))) < 2.0e-14_dp, &
        "Dirichlet diagonal block is symmetric")
    call record_condition(maxval(abs( &
        coupling(4:6, 4:6) - transpose(coupling(4:6, 4:6)))) < 2.0e-14_dp, &
        "Neumann diagonal block is symmetric")
    call record_condition(maxval(abs( &
        coupling(1:3, 4:6) + transpose(coupling(4:6, 1:3)))) < 2.0e-14_dp, &
        "Transmission cross terms cancel in the coupled energy")

    call assemble_laplace_hypersingular_linear( &
        panel_start, panel_end, panel_nodes, 3, 24, hypersingular, status)
    call assemble_laplace_double_layer_mixed_linear( &
        panel_start, panel_end, panel_nodes, 3, 24, double_layer, status)
    call assemble_laplace_single_layer_constant( &
        panel_start, panel_end, 24, single_layer, status)
    mass = 0.0_dp
    boundary_load = 0.0_dp
    gradient = [1.0_dp, 2.0_dp]
    do panel = 1, 3
        length = norm2(panel_end(:, panel) - panel_start(:, panel))
        normal = [panel_end(2, panel) - panel_start(2, panel), &
            panel_start(1, panel) - panel_end(1, panel)] / length
        do endpoint = 1, 2
            mass(panel, panel_nodes(endpoint, panel)) = 0.5_dp * length
            boundary_load(panel_nodes(endpoint, panel)) = &
                boundary_load(panel_nodes(endpoint, panel)) + &
                0.5_dp * length * dot_product(gradient, normal)
        end do
    end do

    trace = 0.75_dp + vertices(1, :) + 2.0_dp * vertices(2, :)
    state = [trace, [(0.0_dp, panel=1, 3)]]
    expected(1:3) = boundary_load + matmul(hypersingular, trace)
    expected(4:6) = matmul(0.5_dp * mass - double_layer, trace)
    call record_condition(maxval(abs(matmul(coupling, state) - expected)) < &
        8.0e-14_dp, "Affine harmonic transmission identity is exact")

    flux = [0.3_dp, -0.2_dp, 0.1_dp]
    state = [trace, flux]
    observed_energy = dot_product(state, matmul(coupling, state))
    expected_energy = 0.5_dp * 0.25_dp**2 * sum(gradient**2) + &
        dot_product(trace, matmul(hypersingular, trace)) + &
        dot_product(flux, matmul(single_layer, flux))
    call record_condition(abs(observed_energy - expected_energy) < 8.0e-14_dp, &
        "Coupled energy equals the independent FEM and BEM energies")

    volume_load = 0.0_dp
    do panel = 1, 3
        length = norm2(panel_end(:, panel) - panel_start(:, panel))
        normal = [panel_end(2, panel) - panel_start(2, panel), &
            panel_start(1, panel) - panel_end(1, panel)] / length
        neumann_jump(panel) = dot_product(gradient, normal)
    end do
    call solve_laplace_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, 24, &
        volume_load, trace, neumann_jump, interior_solution, exterior_flux, &
        status)
    call record_condition(status == 0 .and. &
        maxval(abs(interior_solution - trace)) < 2.0e-12_dp .and. &
        maxval(abs(exterior_flux)) < 2.0e-12_dp, &
        "Transmission solve recovers an affine interior and zero exterior")

    triangles(1, 1) = 0
    call assemble_laplace_symmetric_coupling_p1_p0( &
        vertices, triangles, panel_start, panel_end, panel_nodes, 24, &
        coupling, status)
    call record_condition(status /= 0, "Invalid triangle node is rejected")

    call check_summary("Laplace symmetric FEM-BEM coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_laplace_symmetric_coupling
