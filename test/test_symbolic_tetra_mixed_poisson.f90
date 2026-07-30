program test_symbolic_tetra_mixed_poisson
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        compile_tetra_mixed_form_csc, div, dx, &
        generate_structured_tetra_box_mesh, init_measures, inner, operator(*), &
        solve_symbolic_tetra_mixed_poisson_rt, test_function_t, &
        trial_function_t, vector_test_function_t, vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: divergence
    type(fortsparse_status_t) :: status
    type(test_function_t) :: pressure_test
    type(trial_function_t) :: pressure_trial
    type(vector_test_function_t) :: flux_test
    type(vector_trial_function_t) :: flux_trial
    integer, allocatable :: tetrahedra(:, :)
    real(dp), allocatable :: balance(:), cell_balance(:), flux(:), pressure(:)
    real(dp), allocatable :: vertices(:, :), volumes(:)
    real(dp) :: bounds(3, 2), jacobian(3, 3)
    integer :: cell_count, degree, dof_count, local_status, tetrahedron
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    bounds(:, 1) = 0.0_dp
    bounds(:, 2) = 1.0_dp
    call generate_structured_tetra_box_mesh( &
        bounds, [1, 1, 1], vertices, tetrahedra, local_status)
    if (local_status /= 0) error stop "tetrahedral box mesh failed"

    allocate(volumes(size(tetrahedra, 2)))
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
        volumes(tetrahedron) = abs(det3(jacobian))/6.0_dp
    end do
    do degree = 0, 6
        cell_count = merge(size(tetrahedra, 2), 1, degree <= 5)
        call solve_symbolic_tetra_mixed_poisson_rt( &
            vertices, tetrahedra(:, :cell_count), degree, 2*degree + 4, &
            inner(flux_trial, flux_test)*dx, &
            (-1.0_dp)*inner(pressure_trial, div(flux_test))*dx, &
            inner(div(flux_trial), pressure_test)*dx, unit_source, &
            flux, pressure, status)
        if (status%code /= 0) error stop "symbolic tetra mixed solve failed"
        call compile_tetra_mixed_form_csc( &
            inner(div(flux_trial), pressure_test)*dx, vertices, &
            tetrahedra(:, :cell_count), degree, 2*degree + 4, divergence, &
            status)
        if (status%code /= 0) error stop "symbolic tetra balance block failed"
        balance = csc_matvec(divergence, flux)
        dof_count = (degree + 1)*(degree + 2)*(degree + 3)/6
        allocate(cell_balance(cell_count))
        do tetrahedron = 1, cell_count
            cell_balance(tetrahedron) = &
                balance((tetrahedron - 1)*dof_count + 1)
        end do
        call record_condition( &
            size(pressure) == dof_count*cell_count .and. &
            maxval(abs(cell_balance - volumes(:cell_count))) < 2.0e-9_dp, &
            "symbolic tetra RT-DG solve conserves every cell source")
        call record_condition( &
            abs(sum(cell_balance) - sum(volumes(:cell_count))) < 2.0e-9_dp, &
            "symbolic tetra RT-DG solve has exact domain balance")
        deallocate(balance, cell_balance, flux, pressure)
    end do
    call check_summary("Public symbolic tetrahedral RT-DG Poisson")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function unit_source(x, y, z) result(value)
        real(dp), intent(in) :: x, y, z

        associate(unused => x + y + z)
        end associate
        value = 1.0_dp
    end function unit_source

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_symbolic_tetra_mixed_poisson
