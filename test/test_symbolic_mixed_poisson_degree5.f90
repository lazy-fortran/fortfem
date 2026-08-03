program test_symbolic_mixed_poisson_degree5
    use check, only: check_condition, check_summary
    use fortfem_core, only: &
        function_space_t, mesh_t, rectangle_mesh, vector_function_space_t
    use fortfem_feec, only: &
        compile_mixed_form_csc, div, dx, function_space, init_measures, inner, &
        operator(*), solve_symbolic_mixed_poisson_rt, test_function, &
        test_function_t, trial_function, trial_function_t, vector_function_space, &
        vector_test_function, vector_test_function_t, vector_trial_function, &
        vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    integer, parameter :: degree = 5
    integer, parameter :: local_pressure_count = (degree + 1)*(degree + 2)/2
    type(csc_t) :: divergence
    type(fortsparse_status_t) :: status
    type(mesh_t), target :: mesh
    type(function_space_t), target :: pressure_space
    type(vector_function_space_t), target :: flux_space
    type(test_function_t) :: pressure_test
    type(trial_function_t) :: pressure_trial
    type(vector_test_function_t) :: flux_test
    type(vector_trial_function_t) :: flux_trial
    real(dp), allocatable :: balance(:), flux(:), pressure(:)
    real(dp) :: area, cell_balance
    integer :: cell, first
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    mesh = rectangle_mesh(2, 2, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
    flux_space = vector_function_space(mesh, "RT", degree)
    pressure_space = function_space(mesh, "DG", degree)
    flux_trial = vector_trial_function(flux_space)
    flux_test = vector_test_function(flux_space)
    pressure_trial = trial_function(pressure_space)
    pressure_test = test_function(pressure_space)

    call solve_symbolic_mixed_poisson_rt( &
        mesh%data, degree, 2*degree + 4, &
        inner(flux_trial, flux_test)*dx, &
        (-1.0_dp)*inner(pressure_trial, div(flux_test))*dx, &
        inner(div(flux_trial), pressure_test)*dx, unit_source, &
        flux, pressure, status)
    call record_condition(status%code == 0, &
        "symbolic RT5-DG5 mixed Poisson solve succeeds")
    if (status%code == 0) then
        call compile_mixed_form_csc( &
            inner(div(flux_trial), pressure_test)*dx, mesh%data, "RT", &
            degree, 2*degree + 4, divergence, status)
        balance = csc_matvec(divergence, flux)
        do cell = 1, mesh%data%n_triangles
            first = (cell - 1)*local_pressure_count + 1
            cell_balance = sum(balance(first:first + local_pressure_count - 1))
            area = triangle_area( &
                mesh%data%vertices(:, mesh%data%triangles(:, cell)))
            call record_condition(abs(cell_balance - area) < 2.0e-10_dp, &
                "RT5-DG5 flux reproduces the independent cell source balance")
        end do
    end if
    call check_summary("Symbolic degree-five RT-DG mixed Poisson")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function unit_source(x, y) result(value)
        real(dp), intent(in) :: x, y

        associate(unused => x + y)
        end associate
        value = 1.0_dp
    end function unit_source

    pure real(dp) function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(2, 3)

        area = 0.5_dp*abs( &
            (vertices(1, 2) - vertices(1, 1))* &
            (vertices(2, 3) - vertices(2, 1)) - &
            (vertices(2, 2) - vertices(2, 1))* &
            (vertices(1, 3) - vertices(1, 1)))
    end function triangle_area

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_symbolic_mixed_poisson_degree5
