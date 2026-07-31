program test_mixed_rt_system_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_api, only: solve_mixed_rt_system, solve_mixed_rt_system_jvp, &
        solve_mixed_rt_system_vjp
    use fortfem_kinds, only: dp
    use fortsparse, only: csc_from_triplet, csc_t, fortsparse_status_t
    implicit none

    real(dp), parameter :: step = 2.0e-7_dp
    type(csc_t) :: divergence, divergence_dot, flux_mass, flux_mass_dot
    type(csc_t) :: divergence_minus, divergence_plus
    type(csc_t) :: flux_mass_minus, flux_mass_plus
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: load(2), load_bar(2), load_dot(2)
    real(dp) :: flux(3), flux_bar(3), flux_dot(3)
    real(dp) :: flux_minus(3), flux_plus(3)
    real(dp) :: pressure(2), pressure_bar(2), pressure_dot(2)
    real(dp) :: pressure_minus(2), pressure_plus(2)
    real(dp), allocatable :: divergence_values_bar(:), flux_mass_values_bar(:)
    real(dp) :: lhs, rhs
    integer :: failures, status

    failures = 0
    call make_flux_mass( &
        [4.0_dp, 0.4_dp, 3.2_dp, -0.2_dp, 2.5_dp], flux_mass)
    call make_flux_mass( &
        [0.12_dp, -0.03_dp, 0.08_dp, 0.02_dp, -0.05_dp], flux_mass_dot)
    call make_divergence([1.0_dp, -0.5_dp, 0.3_dp, 0.8_dp], divergence)
    call make_divergence([0.04_dp, -0.02_dp, 0.03_dp, -0.01_dp], &
        divergence_dot)
    load = [1.3_dp, -0.7_dp]
    load_dot = [0.08_dp, -0.11_dp]
    flux_bar = [0.4_dp, -0.2_dp, 0.1_dp]
    pressure_bar = [-0.3_dp, 0.25_dp]

    call solve_mixed_rt_system( &
        flux_mass, divergence, load, flux, pressure, status)
    call check(status == 0, "mixed RT primal solve succeeds")
    call solve_mixed_rt_system_jvp( &
        flux_mass, divergence, load, flux_mass_dot, divergence_dot, load_dot, &
        flux_dot, pressure_dot, status)
    call check(status == 0, "mixed RT JVP succeeds")

    call make_flux_mass( &
        flux_mass%val + step*flux_mass_dot%val, flux_mass_plus)
    call make_flux_mass( &
        flux_mass%val - step*flux_mass_dot%val, flux_mass_minus)
    call make_divergence( &
        divergence%val + step*divergence_dot%val, divergence_plus)
    call make_divergence( &
        divergence%val - step*divergence_dot%val, divergence_minus)
    call solve_mixed_rt_system( &
        flux_mass_plus, divergence_plus, load + step*load_dot, flux_plus, &
        pressure_plus, status)
    call solve_mixed_rt_system( &
        flux_mass_minus, divergence_minus, load - step*load_dot, flux_minus, &
        pressure_minus, status)
    call check(maxval(abs( &
        flux_dot - (flux_plus - flux_minus)/(2.0_dp*step))) < 2.0e-9_dp, &
        "mixed RT flux JVP matches independent central solves")
    call check(maxval(abs( &
        pressure_dot - (pressure_plus - pressure_minus)/(2.0_dp*step))) < &
        1.0e-8_dp, &
        "mixed RT pressure JVP matches independent central solves")

    allocate(flux_mass_values_bar(flux_mass%nnz))
    allocate(divergence_values_bar(divergence%nnz))
    call solve_mixed_rt_system_vjp( &
        flux_mass, divergence, load, flux, pressure, flux_bar, pressure_bar, &
        flux_mass_values_bar, divergence_values_bar, load_bar, status)
    call check(status == 0, "mixed RT VJP succeeds")
    lhs = dot_product(flux_bar, flux_dot) + &
        dot_product(pressure_bar, pressure_dot)
    rhs = dot_product(flux_mass_values_bar, flux_mass_dot%val) + &
        dot_product(divergence_values_bar, divergence_dot%val) + &
        dot_product(load_bar, load_dot)
    call check(abs(lhs - rhs) < 2.0e-11_dp, &
        "mixed RT products satisfy the complete adjoint identity")

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine make_flux_mass(values, matrix)
        real(dp), intent(in) :: values(5)
        type(csc_t), intent(out) :: matrix

        integer, parameter :: rows(5) = [1, 2, 2, 3, 3]
        integer, parameter :: columns(5) = [1, 1, 2, 2, 3]

        call csc_from_triplet( &
            3, 3, rows, columns, values, matrix, sparse_status)
        call check(sparse_status%code == 0, "flux mass fixture is valid CSC")
    end subroutine make_flux_mass

    subroutine make_divergence(values, matrix)
        real(dp), intent(in) :: values(4)
        type(csc_t), intent(out) :: matrix

        integer, parameter :: rows(4) = [1, 2, 1, 2]
        integer, parameter :: columns(4) = [1, 2, 3, 3]

        call csc_from_triplet( &
            2, 3, rows, columns, values, matrix, sparse_status)
        call check(sparse_status%code == 0, "divergence fixture is valid CSC")
    end subroutine make_divergence

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_mixed_rt_system_ad
