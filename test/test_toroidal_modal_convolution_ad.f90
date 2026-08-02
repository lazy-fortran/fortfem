program test_toroidal_modal_convolution_ad
    use, intrinsic :: iso_fortran_env, only: real64
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_toroidal_modal_convolution, &
        apply_toroidal_modal_convolution_jvp, &
        apply_toroidal_modal_convolution_vjp, &
        initialize_fourier_mode_registry, fourier_mode_registry_t
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: dp = real64, mode_count = 5
    integer :: poloidal_modes(mode_count), toroidal_modes(mode_count)
    complex(dp) :: kernel(mode_count), source(mode_count), potential(mode_count)
    complex(dp) :: kernel_dot(mode_count), source_dot(mode_count)
    complex(dp) :: potential_dot(mode_count), potential_plus(mode_count)
    complex(dp) :: kernel_bar(mode_count), source_bar(mode_count)
    complex(dp) :: potential_bar(mode_count)
    complex(dp) :: expected(mode_count), expected_dot(mode_count)
    real(dp) :: epsilon, finite_difference_error, lhs, rhs
    integer :: first_mode, second_mode, output_mode, status
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: sparse_status
    logical :: all_passed

    all_passed = .true.
    poloidal_modes = [-2, -1, 0, 1, 2]
    toroidal_modes = [-1, 0, 0, 0, 1]
    call initialize_fourier_mode_registry( &
        registry, poloidal_modes, toroidal_modes, 1, 0.0_dp, 0.0_dp, .false., &
        status=sparse_status)
    call record_condition(sparse_status%code == 0, &
        "toroidal convolution registry initializes")

    kernel = [ &
        cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.2_dp, dp), &
        cmplx(0.5_dp, 0.0_dp, dp), cmplx(0.4_dp, -0.2_dp, dp), &
        cmplx(-0.1_dp, 0.3_dp, dp)]
    source = [ &
        cmplx(0.6_dp, 0.1_dp, dp), cmplx(-0.4_dp, 0.3_dp, dp), &
        cmplx(0.2_dp, -0.5_dp, dp), cmplx(0.7_dp, 0.2_dp, dp), &
        cmplx(-0.2_dp, 0.4_dp, dp)]
    call apply_toroidal_modal_convolution( &
        registry, kernel, source, potential, sparse_status)
    call record_condition(sparse_status%code == 0, &
        "toroidal convolution action assembles")
    expected = cmplx(0.0_dp, 0.0_dp, dp)
    do output_mode = 1, mode_count
        do first_mode = 1, mode_count
            second_mode = find_mode(registry, &
                poloidal_modes(output_mode) - poloidal_modes(first_mode), &
                toroidal_modes(output_mode) - toroidal_modes(first_mode))
            if (second_mode == 0) cycle
            expected(output_mode) = expected(output_mode) + &
                kernel(second_mode)*source(first_mode)
        end do
    end do
    call record_condition(maxval(abs(potential - expected)) < 1.0e-14_dp, &
        "toroidal convolution matches the independent retained-mode oracle")

    kernel_dot = [ &
        cmplx(-0.02_dp, 0.03_dp, dp), cmplx(0.04_dp, -0.01_dp, dp), &
        cmplx(0.01_dp, 0.02_dp, dp), cmplx(-0.03_dp, 0.01_dp, dp), &
        cmplx(0.02_dp, -0.04_dp, dp)]
    source_dot = [ &
        cmplx(0.03_dp, -0.02_dp, dp), cmplx(-0.01_dp, 0.04_dp, dp), &
        cmplx(0.05_dp, 0.01_dp, dp), cmplx(-0.02_dp, -0.03_dp, dp), &
        cmplx(0.01_dp, 0.05_dp, dp)]
    call apply_toroidal_modal_convolution_jvp( &
        registry, kernel, source, kernel_dot, source_dot, potential_dot, &
        sparse_status)
    call record_condition(sparse_status%code == 0, &
        "toroidal convolution JVP assembles")
    expected_dot = cmplx(0.0_dp, 0.0_dp, dp)
    do output_mode = 1, mode_count
        do first_mode = 1, mode_count
            second_mode = find_mode(registry, &
                poloidal_modes(output_mode) - poloidal_modes(first_mode), &
                toroidal_modes(output_mode) - toroidal_modes(first_mode))
            if (second_mode == 0) cycle
            expected_dot(output_mode) = expected_dot(output_mode) + &
                kernel_dot(second_mode)*source(first_mode) + &
                kernel(second_mode)*source_dot(first_mode)
        end do
    end do
    call record_condition(maxval(abs(potential_dot - expected_dot)) < 1.0e-14_dp, &
        "toroidal convolution JVP matches the independent product oracle")

    epsilon = 1.0e-7_dp
    call apply_toroidal_modal_convolution( &
        registry, kernel + epsilon*kernel_dot, source + epsilon*source_dot, &
        potential_plus, sparse_status)
    finite_difference_error = maxval(abs(potential_dot - &
        (potential_plus - potential)/epsilon))
    call record_condition(finite_difference_error < 2.0e-8_dp, &
        "toroidal convolution JVP matches a forward difference")

    potential_bar = [ &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(-0.3_dp, 0.1_dp, dp), &
        cmplx(0.4_dp, 0.2_dp, dp), cmplx(-0.1_dp, 0.5_dp, dp), &
        cmplx(0.5_dp, -0.2_dp, dp)]
    call apply_toroidal_modal_convolution_vjp( &
        registry, kernel, source, potential_bar, kernel_bar, source_bar, &
        sparse_status)
    call record_condition(sparse_status%code == 0, &
        "toroidal convolution VJP assembles")
    lhs = real(sum(conjg(potential_bar)*potential_dot), dp)
    rhs = real(sum(conjg(kernel_bar)*kernel_dot) + &
        sum(conjg(source_bar)*source_dot), dp)
    call record_condition(abs(lhs - rhs) < 2.0e-12_dp, &
        "toroidal convolution VJP satisfies the real complex adjoint identity")

    call check_summary("toroidal modal convolution")
    if (.not. all_passed) error stop 1

contains

    integer function find_mode(registry, poloidal_mode, toroidal_mode) result(index)
        type(fourier_mode_registry_t), intent(in) :: registry
        integer, intent(in) :: poloidal_mode, toroidal_mode
        integer :: mode

        index = 0
        do mode = 1, size(registry%poloidal_modes)
            if (registry%poloidal_modes(mode) /= poloidal_mode .or. &
                registry%toroidal_modes(mode) /= toroidal_mode) cycle
            index = mode
            return
        end do
    end function find_mode

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_toroidal_modal_convolution_ad
