program test_tetra_mixed_poisson_state_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: solve_tetra_mixed_poisson_state, &
        solve_tetra_mixed_poisson_state_jvp, &
        solve_tetra_mixed_poisson_state_vjp
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: degree = 0, quadrature_degree = 4
    integer, parameter :: tetrahedra(4, 2) = reshape([ &
        1, 2, 3, 4, 2, 3, 4, 5], [4, 2])
    real(dp), parameter :: step = 2.0e-7_dp
    real(dp) :: vertices(3, 5), vertices_dot(3, 5), vertices_bar(3, 5)
    real(dp) :: load(2), load_bar(2), load_dot(2)
    real(dp), allocatable :: flux(:), flux_bar(:), flux_dot(:)
    real(dp), allocatable :: flux_minus(:), flux_plus(:)
    real(dp), allocatable :: pressure(:), pressure_bar(:), pressure_dot(:)
    real(dp), allocatable :: pressure_minus(:), pressure_plus(:)
    real(dp) :: mass_coefficient, mass_coefficient_bar, mass_coefficient_dot
    real(dp) :: lhs, rhs
    integer :: dof, failures, status

    failures = 0
    vertices = reshape([ &
        0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp, 1.0_dp], [3, 5])
    vertices_dot = reshape([ &
        0.01_dp, -0.02_dp, 0.015_dp, &
        -0.03_dp, 0.01_dp, 0.02_dp, &
        0.025_dp, -0.015_dp, 0.005_dp, &
        -0.01_dp, 0.03_dp, -0.02_dp, &
        0.02_dp, -0.01_dp, 0.025_dp], [3, 5])
    load = [1.3_dp, -0.7_dp]
    load_dot = [0.08_dp, -0.11_dp]
    mass_coefficient = 1.6_dp
    mass_coefficient_dot = -0.13_dp

    call solve_tetra_mixed_poisson_state( &
        vertices, tetrahedra, degree, quadrature_degree, mass_coefficient, &
        load, flux, pressure, status)
    call check(status == 0, "tetra mixed Poisson primal state succeeds")
    allocate(flux_bar(size(flux)), pressure_bar(size(pressure)))
    flux_bar = [(0.02_dp*dof - 0.07_dp, dof = 1, size(flux))]
    pressure_bar = [0.3_dp, -0.25_dp]
    call solve_tetra_mixed_poisson_state_jvp( &
        vertices, tetrahedra, degree, quadrature_degree, mass_coefficient, &
        load, vertices_dot, mass_coefficient_dot, load_dot, flux_dot, &
        pressure_dot, status)
    call check(status == 0, "tetra mixed Poisson JVP succeeds")
    call solve_tetra_mixed_poisson_state( &
        vertices + step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        mass_coefficient + step*mass_coefficient_dot, load + step*load_dot, &
        flux_plus, pressure_plus, status)
    call solve_tetra_mixed_poisson_state( &
        vertices - step*vertices_dot, tetrahedra, degree, quadrature_degree, &
        mass_coefficient - step*mass_coefficient_dot, load - step*load_dot, &
        flux_minus, pressure_minus, status)
    call check(maxval(abs( &
        flux_dot - (flux_plus - flux_minus)/(2.0_dp*step))) < 2.0e-8_dp, &
        "tetra mixed flux JVP matches independent full reassembly")
    call check(maxval(abs( &
        pressure_dot - (pressure_plus - pressure_minus)/(2.0_dp*step))) < &
        2.0e-8_dp, &
        "tetra mixed pressure JVP matches independent full reassembly")

    call solve_tetra_mixed_poisson_state_vjp( &
        vertices, tetrahedra, degree, quadrature_degree, mass_coefficient, &
        load, flux, pressure, flux_bar, pressure_bar, vertices_bar, &
        mass_coefficient_bar, load_bar, status)
    call check(status == 0, "tetra mixed Poisson VJP succeeds")
    lhs = dot_product(flux_bar, flux_dot) + &
        dot_product(pressure_bar, pressure_dot)
    rhs = sum(vertices_bar*vertices_dot) + &
        mass_coefficient_bar*mass_coefficient_dot + &
        dot_product(load_bar, load_dot)
    call check(abs(lhs - rhs) < 3.0e-10_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "tetra mixed products accumulate the moving-mesh adjoint")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_tetra_mixed_poisson_state_ad
