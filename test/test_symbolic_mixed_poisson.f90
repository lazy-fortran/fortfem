program test_symbolic_mixed_poisson
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        compile_mixed_form_csc, div, dx, init_measures, inner, &
        function_space, function_space_t, mesh_t, operator(*), rectangle_mesh, &
        solve_symbolic_mixed_poisson_rt, &
        test_function, test_function_t, trial_function, trial_function_t, &
        vector_function_space, vector_function_space_t, vector_test_function, &
        vector_test_function_t, vector_trial_function, vector_trial_function_t
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_matvec, csc_t, fortsparse_status_t
    implicit none

    type(csc_t) :: divergence
    type(fortsparse_status_t) :: status
    type(mesh_t), target :: mesh
    type(test_function_t) :: pressure_test
    type(trial_function_t) :: pressure_trial
    type(function_space_t), target :: pressure_space
    type(vector_function_space_t), target :: flux_space
    type(vector_test_function_t) :: flux_test
    type(vector_trial_function_t) :: flux_trial
    real(dp), allocatable :: balance(:), flux(:), pressure(:)
    real(dp), allocatable :: expected_load(:)
    real(dp) :: a(2), b(2), c(2)
    integer :: triangle
    logical :: all_passed

    all_passed = .true.
    call init_measures()
    mesh = rectangle_mesh(2, 2, [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
    flux_space = vector_function_space(mesh, "RT", 0)
    pressure_space = function_space(mesh, "DG", 0)
    flux_trial = vector_trial_function(flux_space)
    flux_test = vector_test_function(flux_space)
    pressure_trial = trial_function(pressure_space)
    pressure_test = test_function(pressure_space)

    call solve_symbolic_mixed_poisson_rt( &
        mesh%data, 0, 4, &
        inner(flux_trial, flux_test)*dx, &
        (-1.0_dp)*inner(pressure_trial, div(flux_test))*dx, &
        inner(div(flux_trial), pressure_test)*dx, unit_source, &
        flux, pressure, status)
    if (status%code /= 0) error stop "symbolic mixed Poisson solve failed"

    call compile_mixed_form_csc( &
        inner(div(flux_trial), pressure_test)*dx, mesh%data, "RT", 0, 4, &
        divergence, status)
    if (status%code /= 0) error stop "symbolic balance block failed"
    balance = csc_matvec(divergence, flux)
    allocate(expected_load(mesh%data%n_triangles))
    do triangle = 1, mesh%data%n_triangles
        a = mesh%data%vertices(:, mesh%data%triangles(1, triangle))
        b = mesh%data%vertices(:, mesh%data%triangles(2, triangle))
        c = mesh%data%vertices(:, mesh%data%triangles(3, triangle))
        expected_load(triangle) = 0.5_dp*abs( &
            (b(1) - a(1))*(c(2) - a(2)) - &
            (b(2) - a(2))*(c(1) - a(1)))
    end do
    call record_condition( &
        size(pressure) == mesh%data%n_triangles .and. &
        maxval(abs(balance - expected_load)) < 2.0e-12_dp, &
        "symbolic mixed solve conserves the unit source in every cell")
    call record_condition( &
        abs(sum(balance) - 1.0_dp) < 2.0e-12_dp, &
        "symbolic mixed solve has the exact global flux balance")
    call check_summary("Public symbolic RT-DG mixed Poisson")
    if (.not. all_passed) error stop 1

contains

    pure real(dp) function unit_source(x, y) result(value)
        real(dp), intent(in) :: x, y

        associate(unused => x + y)
        end associate
        value = 1.0_dp
    end function unit_source

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_symbolic_mixed_poisson
