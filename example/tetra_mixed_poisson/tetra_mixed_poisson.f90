program tetra_mixed_poisson
    use fortfem_api, only: &
        compile_tetra_mixed_form_csc, div, dx, &
        generate_structured_tetra_box_mesh, init_measures, inner, operator(*), &
        solve_symbolic_tetra_mixed_poisson_rt, test_function_t, &
        trial_function_t, vector_test_function_t, vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortplot, only: figure, legend, plot, savefig, set_yscale, title, &
        xlabel, ylabel
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/tetra_mixed_poisson"
    type(csc_t) :: divergence
    type(fortsparse_status_t) :: status
    type(test_function_t) :: pressure_test
    type(trial_function_t) :: pressure_trial
    type(vector_test_function_t) :: flux_test
    type(vector_trial_function_t) :: flux_trial
    integer, allocatable :: tetrahedra(:, :)
    real(dp), allocatable :: balance(:), flux(:), pressure(:), vertices(:, :)
    real(dp) :: balance_error(5), bounds(3, 2), cell_balance
    real(dp) :: dofs(5), jacobian(3, 3), orders(5), volume
    integer :: degree, dg_count, local_status, tetrahedron

    call init_measures()
    call execute_command_line("mkdir -p "//output_directory)
    bounds(:, 1) = 0.0_dp
    bounds(:, 2) = 1.0_dp
    call generate_structured_tetra_box_mesh( &
        bounds, [1, 1, 1], vertices, tetrahedra, local_status)
    if (local_status /= 0) error stop "tetrahedral box mesh failed"

    do degree = 0, 4
        orders(degree + 1) = real(degree, dp)
        call solve_symbolic_tetra_mixed_poisson_rt( &
            vertices, tetrahedra, degree, 2*degree + 4, &
            inner(flux_trial, flux_test)*dx, &
            (-1.0_dp)*inner(pressure_trial, div(flux_test))*dx, &
            inner(div(flux_trial), pressure_test)*dx, unit_source, &
            flux, pressure, status)
        if (status%code /= 0) error stop "symbolic tetra solve failed"
        call compile_tetra_mixed_form_csc( &
            inner(div(flux_trial), pressure_test)*dx, vertices, tetrahedra, &
            degree, 2*degree + 4, divergence, status)
        if (status%code /= 0) error stop "symbolic tetra block failed"
        balance = csc_matvec(divergence, flux)
        dg_count = (degree + 1)*(degree + 2)*(degree + 3)/6
        balance_error(degree + 1) = 0.0_dp
        do tetrahedron = 1, size(tetrahedra, 2)
            jacobian(:, 1) = &
                vertices(:, tetrahedra(2, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 2) = &
                vertices(:, tetrahedra(3, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            jacobian(:, 3) = &
                vertices(:, tetrahedra(4, tetrahedron)) - &
                vertices(:, tetrahedra(1, tetrahedron))
            volume = abs(det3(jacobian))/6.0_dp
            cell_balance = balance((tetrahedron - 1)*dg_count + 1)
            balance_error(degree + 1) = max( &
                balance_error(degree + 1), abs(cell_balance - volume))
        end do
        dofs(degree + 1) = real(size(flux) + size(pressure), dp)
        if (balance_error(degree + 1) >= 2.0e-9_dp) &
            error stop "tetrahedral cell balance regression"
        deallocate(balance, flux, pressure)
    end do

    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(orders, max(balance_error, epsilon(1.0_dp)), &
        label="maximum cell balance error", marker="o")
    call set_yscale("log")
    call xlabel("RT-DG polynomial degree")
    call ylabel("absolute balance error")
    call title("Symbolic tetrahedral mixed Poisson conservation")
    call legend()
    call savefig(output_directory//"/tetra_mixed_conservation_1d.png")

    call figure(figsize=[8.5_dp, 5.5_dp])
    call plot(orders, dofs, label="coupled unknowns", marker="s")
    call xlabel("RT-DG polynomial degree")
    call ylabel("degrees of freedom")
    call title("Tetrahedral mixed-system order growth")
    call legend()
    call savefig(output_directory//"/tetra_mixed_dofs_1d.png")

contains

    pure real(dp) function unit_source(x, y, z) result(value)
        real(dp), intent(in) :: x, y, z

        associate(unused => x + y + z)
        end associate
        value = 1.0_dp
    end function unit_source

end program tetra_mixed_poisson
