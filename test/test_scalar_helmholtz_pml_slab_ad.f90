program test_scalar_helmholtz_pml_slab_ad
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        solve_scalar_helmholtz_pml_slab_1d, &
        solve_scalar_helmholtz_pml_slab_1d_jvp, &
        solve_scalar_helmholtz_pml_slab_1d_vjp
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: step = 2.0e-6_dp
    real(dp) :: nodes(5), nodes_dot(5), nodes_bar(5)
    real(dp) :: physical_end, physical_end_dot, physical_end_bar
    real(dp) :: wave_number, wave_number_dot, wave_number_bar
    real(dp) :: sigma_max, sigma_max_dot, sigma_max_bar
    complex(dp) :: left_value, left_value_dot, left_value_bar
    complex(dp), allocatable :: solution(:), solution_dot(:)
    complex(dp), allocatable :: solution_plus(:), solution_minus(:)
    complex(dp) :: solution_bar(5), lhs, rhs
    integer :: polynomial_degree, status, status_minus, status_plus

    nodes = [0.0_dp, 0.4_dp, 0.8_dp, 1.2_dp, 1.6_dp]
    nodes_dot = [0.006_dp, -0.004_dp, 0.008_dp, -0.005_dp, 0.007_dp]
    physical_end = 0.8_dp
    physical_end_dot = 0.012_dp
    wave_number = 2.0_dp
    wave_number_dot = 0.05_dp
    sigma_max = 3.0_dp
    sigma_max_dot = -0.08_dp
    polynomial_degree = 2
    left_value = cmplx(1.0_dp, -0.2_dp, dp)
    left_value_dot = cmplx(-0.03_dp, 0.04_dp, dp)
    solution_bar = [ &
        cmplx(0.0_dp, 0.0_dp, dp), cmplx(0.2_dp, -0.1_dp, dp), &
        cmplx(-0.1_dp, 0.3_dp, dp), cmplx(0.4_dp, 0.2_dp, dp), &
        cmplx(0.0_dp, 0.0_dp, dp)]

    call solve_scalar_helmholtz_pml_slab_1d( &
        nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
        left_value, solution, status)
    call check_condition(status == 0, "PML slab primal solve succeeds")
    call check_condition(size(solution) == size(nodes), &
        "PML slab primal solution has nodal size")

    call solve_scalar_helmholtz_pml_slab_1d_jvp( &
        nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
        left_value, nodes_dot, physical_end_dot, wave_number_dot, &
        sigma_max_dot, left_value_dot, solution_dot, status)
    call solve_scalar_helmholtz_pml_slab_1d( &
        nodes + step*nodes_dot, physical_end + step*physical_end_dot, &
        wave_number + step*wave_number_dot, sigma_max + step*sigma_max_dot, &
        polynomial_degree, left_value + step*left_value_dot, solution_plus, &
        status_plus)
    call solve_scalar_helmholtz_pml_slab_1d( &
        nodes - step*nodes_dot, physical_end - step*physical_end_dot, &
        wave_number - step*wave_number_dot, sigma_max - step*sigma_max_dot, &
        polynomial_degree, left_value - step*left_value_dot, solution_minus, &
        status_minus)
    call check_condition(status == 0, "PML slab JVP succeeds")
    call check_condition(status_plus == 0, "positive PML slab perturbation succeeds")
    call check_condition(status_minus == 0, "negative PML slab perturbation succeeds")
    call check_condition( &
        maxval(abs(solution_dot - (solution_plus - solution_minus)/(2.0_dp*step))) &
        < 2.0e-7_dp, "PML slab JVP matches central differences")

    call solve_scalar_helmholtz_pml_slab_1d_vjp( &
        nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
        left_value, solution, solution_bar, nodes_bar, physical_end_bar, &
        wave_number_bar, sigma_max_bar, left_value_bar, status)
    lhs = sum(conjg(solution_bar)*solution_dot)
    rhs = cmplx(sum(nodes_bar*nodes_dot) + physical_end_bar*physical_end_dot + &
        wave_number_bar*wave_number_dot + sigma_max_bar*sigma_max_dot, 0.0_dp, dp) + &
        conjg(left_value_bar)*left_value_dot
    call check_condition(status == 0, "PML slab VJP succeeds")
    call check_condition( &
        abs(real(lhs - rhs, dp)) < 3.0e-10_dp*max(1.0_dp, abs(real(lhs, dp))), &
        "PML slab products obey the real adjoint identity")

    call check_summary("Differentiable scalar Helmholtz PML slab")
end program test_scalar_helmholtz_pml_slab_ad
