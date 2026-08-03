program test_triangle_vector_basis_ad
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortfem_feec, only: &
        evaluate_triangle_bdm, evaluate_triangle_bdm_jvp, &
        evaluate_triangle_bdm_vjp, evaluate_triangle_nedelec_first_kind, &
        evaluate_triangle_nedelec_first_kind_jvp, &
        evaluate_triangle_nedelec_first_kind_vjp, &
        evaluate_triangle_nedelec_second_kind, &
        evaluate_triangle_nedelec_second_kind_jvp, &
        evaluate_triangle_nedelec_second_kind_vjp, &
        evaluate_triangle_raviart_thomas, &
        evaluate_triangle_raviart_thomas_jvp, &
        evaluate_triangle_raviart_thomas_vjp, initialize_triangle_bdm, &
        initialize_triangle_nedelec_first_kind, &
        initialize_triangle_nedelec_second_kind, &
        initialize_triangle_raviart_thomas, triangle_bdm_basis_t, &
        triangle_bdm_dof_count, triangle_nedelec_dof_count, &
        triangle_nedelec_first_kind_t, &
        triangle_nedelec_second_kind_dof_count, &
        triangle_nedelec_second_kind_t, triangle_rt_basis_t, &
        triangle_rt_dof_count
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: eta = 0.27_dp, eta_dot = -0.19_dp
    real(dp), parameter :: step = 2.0e-6_dp
    real(dp), parameter :: xi = 0.31_dp, xi_dot = 0.23_dp
    type(triangle_bdm_basis_t) :: bdm_basis
    type(triangle_nedelec_first_kind_t) :: nedelec_first_basis
    type(triangle_nedelec_second_kind_t) :: nedelec_second_basis
    type(triangle_rt_basis_t) :: rt_basis
    real(dp), allocatable :: scalars(:), scalars_bar(:), scalars_dot(:)
    real(dp), allocatable :: scalars_minus(:), scalars_plus(:)
    real(dp), allocatable :: values(:, :), values_bar(:, :), values_dot(:, :)
    real(dp), allocatable :: values_minus(:, :), values_plus(:, :)
    real(dp) :: eta_bar, xi_bar
    integer :: failures, status

    failures = 0
    call initialize_triangle_nedelec_first_kind(3, nedelec_first_basis, status)
    call initialize_arrays(triangle_nedelec_dof_count(nedelec_first_basis))
    call check_nedelec_first()
    call initialize_triangle_raviart_thomas(2, rt_basis, status)
    call initialize_arrays(triangle_rt_dof_count(rt_basis))
    call check_rt()
    call initialize_triangle_nedelec_second_kind(2, nedelec_second_basis, status)
    call initialize_arrays( &
        triangle_nedelec_second_kind_dof_count(nedelec_second_basis))
    call check_nedelec_second()
    call initialize_triangle_bdm(2, bdm_basis, status)
    call initialize_arrays(triangle_bdm_dof_count(bdm_basis))
    call check_bdm()

    if (failures > 0) then
        write(error_unit, "(i0,a)") failures, " test(s) failed"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine initialize_arrays(count)
        integer, intent(in) :: count
        integer :: dof

        if (allocated(scalars)) then
            deallocate( &
                scalars, scalars_bar, scalars_dot, scalars_minus, scalars_plus)
            deallocate( &
                values, values_bar, values_dot, values_minus, values_plus)
        end if
        allocate(scalars(count), scalars_bar(count), scalars_dot(count))
        allocate(scalars_minus(count), scalars_plus(count))
        allocate(values(2, count), values_bar(2, count), values_dot(2, count))
        allocate(values_minus(2, count), values_plus(2, count))
        do dof = 1, count
            values_bar(:, dof) = [0.02_dp*dof - 0.11_dp, 0.07_dp - 0.01_dp*dof]
            scalars_bar(dof) = 0.015_dp*dof - 0.04_dp
        end do
    end subroutine initialize_arrays

    subroutine check_nedelec_first()
        call evaluate_triangle_nedelec_first_kind_jvp( &
            nedelec_first_basis, xi, eta, xi_dot, eta_dot, values_dot, &
            scalars_dot, status)
        call evaluate_triangle_nedelec_first_kind( &
            nedelec_first_basis, xi + step*xi_dot, eta + step*eta_dot, &
            values_plus, scalars_plus, status)
        call evaluate_triangle_nedelec_first_kind( &
            nedelec_first_basis, xi - step*xi_dot, eta - step*eta_dot, &
            values_minus, scalars_minus, status)
        call evaluate_triangle_nedelec_first_kind_vjp( &
            nedelec_first_basis, xi, eta, values_bar, scalars_bar, xi_bar, &
            eta_bar, status)
        call check_products("first-kind Nedelec")
    end subroutine check_nedelec_first

    subroutine check_rt()
        call evaluate_triangle_raviart_thomas_jvp( &
            rt_basis, xi, eta, xi_dot, eta_dot, values_dot, scalars_dot, status)
        call evaluate_triangle_raviart_thomas( &
            rt_basis, xi + step*xi_dot, eta + step*eta_dot, values_plus, &
            scalars_plus, status)
        call evaluate_triangle_raviart_thomas( &
            rt_basis, xi - step*xi_dot, eta - step*eta_dot, values_minus, &
            scalars_minus, status)
        call evaluate_triangle_raviart_thomas_vjp( &
            rt_basis, xi, eta, values_bar, scalars_bar, xi_bar, eta_bar, status)
        call check_products("Raviart-Thomas")
    end subroutine check_rt

    subroutine check_nedelec_second()
        call evaluate_triangle_nedelec_second_kind_jvp( &
            nedelec_second_basis, xi, eta, xi_dot, eta_dot, values_dot, &
            scalars_dot, status)
        call evaluate_triangle_nedelec_second_kind( &
            nedelec_second_basis, xi + step*xi_dot, eta + step*eta_dot, &
            values_plus, scalars_plus, status)
        call evaluate_triangle_nedelec_second_kind( &
            nedelec_second_basis, xi - step*xi_dot, eta - step*eta_dot, &
            values_minus, scalars_minus, status)
        call evaluate_triangle_nedelec_second_kind_vjp( &
            nedelec_second_basis, xi, eta, values_bar, scalars_bar, xi_bar, &
            eta_bar, status)
        call check_products("second-kind Nedelec")
    end subroutine check_nedelec_second

    subroutine check_bdm()
        call evaluate_triangle_bdm_jvp( &
            bdm_basis, xi, eta, xi_dot, eta_dot, values_dot, scalars_dot, status)
        call evaluate_triangle_bdm( &
            bdm_basis, xi + step*xi_dot, eta + step*eta_dot, values_plus, &
            scalars_plus, status)
        call evaluate_triangle_bdm( &
            bdm_basis, xi - step*xi_dot, eta - step*eta_dot, values_minus, &
            scalars_minus, status)
        call evaluate_triangle_bdm_vjp( &
            bdm_basis, xi, eta, values_bar, scalars_bar, xi_bar, eta_bar, status)
        call check_products("BDM")
    end subroutine check_bdm

    subroutine check_products(label)
        character(*), intent(in) :: label
        real(dp) :: lhs, rhs

        call check(status == 0, label//" basis products succeed")
        call check(maxval(abs( &
            values_dot - (values_plus - values_minus)/(2.0_dp*step))) < &
            5.0e-8_dp, label//" value JVP matches re-evaluation")
        call check(maxval(abs( &
            scalars_dot - (scalars_plus - scalars_minus)/(2.0_dp*step))) < &
            5.0e-8_dp, label//" derivative JVP matches re-evaluation")
        lhs = sum(values_bar*values_dot) + dot_product(scalars_bar, scalars_dot)
        rhs = xi_bar*xi_dot + eta_bar*eta_dot
        call check(abs(lhs - rhs) < 2.0e-12_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
            label//" basis products obey adjoint identity")
    end subroutine check_products

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write(error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_triangle_vector_basis_ad
