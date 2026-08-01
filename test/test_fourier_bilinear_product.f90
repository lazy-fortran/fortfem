program test_fourier_bilinear_product
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        apply_fourier_bilinear_product, &
        apply_fourier_bilinear_product_jvp, &
        apply_fourier_bilinear_product_vjp, &
        build_fourier_mode_padded_registry, &
        find_fourier_mode, fourier_mode_registry_t, &
        initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: input_mode_count = 3, point_count = 2
    integer, parameter :: left_count = 2, right_count = 2, output_count = 2
    integer, parameter :: input_poloidal(input_mode_count) = [-1, 0, 1]
    integer, parameter :: input_toroidal(input_mode_count) = [0, 0, 0]
    real(dp), parameter :: eps = 1.0e-7_dp
    real(dp), parameter :: input_phase = 0.13_dp
    real(dp), parameter :: time_phase = -0.17_dp
    type(fourier_mode_registry_t) :: input_registry, output_registry
    type(fortsparse_status_t) :: status
    complex(dp), allocatable :: coupling(:, :, :), coupling_dot(:, :, :)
    complex(dp), allocatable :: coupling_bar(:, :, :)
    complex(dp), allocatable :: left(:, :, :), right(:, :, :)
    complex(dp), allocatable :: left_dot(:, :, :), right_dot(:, :, :)
    complex(dp), allocatable :: left_bar(:, :, :), right_bar(:, :, :)
    complex(dp), allocatable :: product(:, :, :), product_dot(:, :, :)
    complex(dp), allocatable :: product_plus(:, :, :), product_minus(:, :, :)
    complex(dp), allocatable :: product_bar(:, :, :), oracle(:, :, :)
    complex(dp), allocatable :: oracle_dot(:, :, :), bad_product(:, :, :)
    integer :: output_mode_count
    integer :: first_mode, second_mode, output_mode, point
    integer :: component, left_component, right_component
    complex(dp) :: term, term_dot
    real(dp) :: lhs, rhs

    call initialize_fourier_mode_registry( &
        input_registry, input_poloidal, input_toroidal, 2, input_phase, &
        time_phase, .false., status=status)
    call check_condition(status%code == 0, &
        "bilinear Fourier application accepts the input registry")
    call build_fourier_mode_padded_registry( &
        input_registry, output_registry, status)
    if (status%code /= 0) error stop 1
    output_mode_count = size(output_registry%poloidal_modes)
    call check_condition(output_mode_count > input_mode_count, &
        "bilinear Fourier application receives a padded output registry")

    allocate(coupling(output_count, left_count, right_count), &
        coupling_dot(output_count, left_count, right_count), &
        coupling_bar(output_count, left_count, right_count), &
        left(left_count, point_count, input_mode_count), &
        right(right_count, point_count, input_mode_count), &
        left_dot(left_count, point_count, input_mode_count), &
        right_dot(right_count, point_count, input_mode_count), &
        left_bar(left_count, point_count, input_mode_count), &
        right_bar(right_count, point_count, input_mode_count), &
        product(output_count, point_count, output_mode_count), &
        product_dot(output_count, point_count, output_mode_count), &
        product_plus(output_count, point_count, output_mode_count), &
        product_minus(output_count, point_count, output_mode_count), &
        product_bar(output_count, point_count, output_mode_count), &
        oracle(output_count, point_count, output_mode_count), &
        oracle_dot(output_count, point_count, output_mode_count), &
        bad_product(1, point_count, output_mode_count))

    do component = 1, output_count
        do left_component = 1, left_count
            do right_component = 1, right_count
                coupling(component, left_component, right_component) = cmplx( &
                    0.11_dp*real(component + left_component + right_component, dp), &
                    -0.03_dp*real(component + 2*left_component + right_component, dp), dp)
                coupling_dot(component, left_component, right_component) = cmplx( &
                    -0.02_dp*real(component + left_component, dp), &
                    0.04_dp*real(right_component, dp), dp)
            end do
        end do
    end do
    do first_mode = 1, input_mode_count
        do point = 1, point_count
            do left_component = 1, left_count
                left(left_component, point, first_mode) = cmplx( &
                    0.2_dp*real(left_component + point + first_mode, dp), &
                    -0.07_dp*real(left_component + first_mode, dp), dp)
                left_dot(left_component, point, first_mode) = cmplx( &
                    -0.03_dp*real(point + first_mode, dp), &
                    0.02_dp*real(left_component, dp), dp)
            end do
            do right_component = 1, right_count
                right(right_component, point, first_mode) = cmplx( &
                    -0.13_dp*real(right_component + point, dp), &
                    0.06_dp*real(right_component + first_mode, dp), dp)
                right_dot(right_component, point, first_mode) = cmplx( &
                    0.05_dp*real(right_component + first_mode, dp), &
                    -0.01_dp*real(point, dp), dp)
            end do
        end do
    end do

    oracle = cmplx(0.0_dp, 0.0_dp, dp)
    oracle_dot = cmplx(0.0_dp, 0.0_dp, dp)
    do first_mode = 1, input_mode_count
        do second_mode = 1, input_mode_count
            output_mode = find_fourier_mode(output_registry, &
                input_poloidal(first_mode) + input_poloidal(second_mode), &
                input_toroidal(first_mode) + input_toroidal(second_mode))
            if (output_mode == 0) error stop 1
            do point = 1, point_count
                do component = 1, output_count
                    term = cmplx(0.0_dp, 0.0_dp, dp)
                    term_dot = cmplx(0.0_dp, 0.0_dp, dp)
                    do left_component = 1, left_count
                        do right_component = 1, right_count
                            term = term + coupling(component, left_component, &
                                right_component)*left(left_component, point, first_mode)* &
                                right(right_component, point, second_mode)
                            term_dot = term_dot + &
                                coupling_dot(component, left_component, right_component)* &
                                left(left_component, point, first_mode)* &
                                right(right_component, point, second_mode) + &
                                coupling(component, left_component, right_component)* &
                                left_dot(left_component, point, first_mode)* &
                                right(right_component, point, second_mode) + &
                                coupling(component, left_component, right_component)* &
                                left(left_component, point, first_mode)* &
                                right_dot(right_component, point, second_mode)
                        end do
                    end do
                    oracle(component, point, output_mode) = &
                        oracle(component, point, output_mode) + term
                    oracle_dot(component, point, output_mode) = &
                        oracle_dot(component, point, output_mode) + term_dot
                end do
            end do
        end do
    end do

    call apply_fourier_bilinear_product( &
        input_registry, output_registry, coupling, left, right, product, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(product - oracle)) < 1.0e-13_dp, &
        "bilinear Fourier application matches an independent padded oracle")
    call apply_fourier_bilinear_product_jvp( &
        input_registry, output_registry, coupling, left, right, coupling_dot, &
        left_dot, right_dot, product_dot, status)
    call apply_fourier_bilinear_product( &
        input_registry, output_registry, coupling + eps*coupling_dot, &
        left + eps*left_dot, right + eps*right_dot, product_plus, status)
    call apply_fourier_bilinear_product( &
        input_registry, output_registry, coupling - eps*coupling_dot, &
        left - eps*left_dot, right - eps*right_dot, product_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(product_dot - oracle_dot)) < 1.0e-13_dp .and. &
        maxval(abs((product_plus - product_minus)/(2.0_dp*eps) - product_dot)) &
        < 1.0e-8_dp, "bilinear Fourier application JVP matches product rule")

    do component = 1, output_count
        do point = 1, point_count
            do output_mode = 1, output_mode_count
                product_bar(component, point, output_mode) = cmplx( &
                    0.06_dp*real(component + point + output_mode, dp), &
                    -0.02_dp*real(component + output_mode, dp), dp)
            end do
        end do
    end do
    call apply_fourier_bilinear_product_vjp( &
        input_registry, output_registry, coupling, left, right, product_bar, &
        left_bar, right_bar, coupling_bar, status)
    lhs = real(sum(conjg(product_bar)*product_dot), dp)
    rhs = real(sum(conjg(left_bar)*left_dot) + sum(conjg(right_bar)*right_dot) + &
        sum(conjg(coupling_bar)*coupling_dot), dp)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "bilinear Fourier application VJP satisfies the complex dot identity")

    call apply_fourier_bilinear_product( &
        input_registry, output_registry, coupling, left, right, bad_product, status)
    call check_condition(status%code /= 0, &
        "bilinear Fourier application rejects an incompatible output shape")
    call check_summary("padded bilinear Fourier application")
end program test_fourier_bilinear_product
