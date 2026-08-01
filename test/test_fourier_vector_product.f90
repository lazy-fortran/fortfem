program test_fourier_vector_product
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_fourier_vector_product, &
        assemble_fourier_vector_product_jvp, &
        assemble_fourier_vector_product_vjp, &
        build_fourier_mode_triad_map, fourier_mode_registry_t, &
        initialize_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: mode_count = 3, point_count = 2
    integer, parameter :: left_count = 2, right_count = 2, output_count = 2
    integer, parameter :: modes(mode_count) = [-1, 0, 1]
    integer, parameter :: expected_triad_map(3, 3) = reshape([ &
        0, 1, 2, 1, 2, 3, 2, 3, 0], [3, 3])
    integer :: first_mode, second_mode, output_mode, point, left_component
    integer :: right_component, component
    real(dp), parameter :: eps = 1.0e-7_dp
    type(fourier_mode_registry_t) :: registry
    type(fortsparse_status_t) :: status
    complex(dp) :: coupling(output_count, left_count, right_count)
    complex(dp) :: coupling_dot(output_count, left_count, right_count)
    complex(dp) :: coupling_bar(output_count, left_count, right_count)
    complex(dp) :: left(left_count, point_count, mode_count)
    complex(dp) :: right(right_count, point_count, mode_count)
    complex(dp) :: left_dot(left_count, point_count, mode_count)
    complex(dp) :: right_dot(right_count, point_count, mode_count)
    complex(dp) :: left_bar(left_count, point_count, mode_count)
    complex(dp) :: right_bar(right_count, point_count, mode_count)
    complex(dp) :: product(output_count, point_count, mode_count)
    complex(dp) :: product_dot(output_count, point_count, mode_count)
    complex(dp) :: product_plus(output_count, point_count, mode_count)
    complex(dp) :: product_minus(output_count, point_count, mode_count)
    complex(dp) :: product_bar(output_count, point_count, mode_count)
    complex(dp) :: oracle(output_count, point_count, mode_count)
    complex(dp) :: oracle_dot(output_count, point_count, mode_count)
    complex(dp) :: term, term_dot
    complex(dp) :: bad_product(1, point_count, mode_count)
    integer, allocatable :: triad_map(:, :)
    integer :: missing_triads
    real(dp) :: lhs, rhs

    call initialize_fourier_mode_registry( &
        registry, modes, modes, 1, 0.0_dp, 0.0_dp, .false., status=status)
    call check_condition(status%code == 0, &
        "Fourier vector product accepts the retained mode registry")
    call build_fourier_mode_triad_map( &
        registry, triad_map, missing_triads, status)
    if (status%code /= 0) error stop 1
    if (missing_triads /= 2) error stop 1
    if (.not. all(triad_map == expected_triad_map)) error stop 1
    call check_condition(status%code == 0 .and. missing_triads == 2 .and. &
        all(triad_map == expected_triad_map), &
        "Fourier triad map reports retained and omitted sums independently")

    do component = 1, output_count
        do left_component = 1, left_count
            do right_component = 1, right_count
                coupling(component, left_component, right_component) = &
                    cmplx(0.1_dp*real( &
                    component + left_component + right_component, dp), &
                    -0.04_dp*real( &
                    component + 2*left_component + right_component, dp), dp)
                coupling_dot(component, left_component, right_component) = &
                    cmplx(-0.02_dp*real(component + left_component, dp), &
                    0.03_dp*real(right_component, dp), dp)
            end do
        end do
    end do
    do output_mode = 1, mode_count
        do point = 1, point_count
            do left_component = 1, left_count
                left(left_component, point, output_mode) = cmplx( &
                    0.2_dp*real(left_component + point + output_mode, dp), &
                    -0.1_dp*real(left_component + output_mode, dp), dp)
                left_dot(left_component, point, output_mode) = cmplx( &
                    -0.03_dp*real(point + output_mode, dp), &
                    0.02_dp*real(left_component, dp), dp)
            end do
            do right_component = 1, right_count
                right(right_component, point, output_mode) = cmplx( &
                    -0.15_dp*real(right_component + point, dp), &
                    0.07_dp*real(right_component + output_mode, dp), dp)
                right_dot(right_component, point, output_mode) = cmplx( &
                    0.04_dp*real(right_component + output_mode, dp), &
                    -0.01_dp*real(point, dp), dp)
            end do
        end do
    end do

    oracle = cmplx(0.0_dp, 0.0_dp, dp)
    oracle_dot = cmplx(0.0_dp, 0.0_dp, dp)
    do first_mode = 1, mode_count
        do second_mode = 1, mode_count
            output_mode = find_mode(modes, modes(first_mode) + modes(second_mode))
            if (output_mode == 0) cycle
            do point = 1, point_count
                do component = 1, output_count
                    term = cmplx(0.0_dp, 0.0_dp, dp)
                    term_dot = cmplx(0.0_dp, 0.0_dp, dp)
                    do left_component = 1, left_count
                        do right_component = 1, right_count
                            term = term + coupling(component, left_component, &
                                right_component)*left( &
                                left_component, point, first_mode)* &
                                right(right_component, point, second_mode)
                            term_dot = term_dot + &
                                coupling_dot(component, left_component, &
                                right_component)* &
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

    call assemble_fourier_vector_product( &
        registry, coupling, left, right, product, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(product - oracle)) < 1.0e-13_dp, &
        "Fourier vector product matches an independent triad oracle")
    call assemble_fourier_vector_product_jvp( &
        registry, coupling, left, right, coupling_dot, left_dot, right_dot, &
        product_dot, status)
    call assemble_fourier_vector_product( &
        registry, coupling + eps*coupling_dot, left + eps*left_dot, &
        right + eps*right_dot, product_plus, status)
    call assemble_fourier_vector_product( &
        registry, coupling - eps*coupling_dot, left - eps*left_dot, &
        right - eps*right_dot, product_minus, status)
    call check_condition(status%code == 0 .and. &
        maxval(abs(product_dot - oracle_dot)) < 1.0e-13_dp .and. &
        maxval(abs((product_plus - product_minus)/(2.0_dp*eps) - product_dot)) &
        < 1.0e-8_dp, "Fourier vector-product JVP matches product rule")

    do component = 1, output_count
        do point = 1, point_count
            do output_mode = 1, mode_count
                product_bar(component, point, output_mode) = cmplx( &
                    0.06_dp*real(component + point + output_mode, dp), &
                    -0.02_dp*real(component + output_mode, dp), dp)
            end do
        end do
    end do
    call assemble_fourier_vector_product_vjp( &
        registry, coupling, left, right, product_bar, left_bar, right_bar, &
        coupling_bar, status)
    lhs = real(sum(conjg(product_bar)*product_dot), dp)
    rhs = real(sum(conjg(left_bar)*left_dot) + sum(conjg(right_bar)*right_dot) + &
        sum(conjg(coupling_bar)*coupling_dot), dp)
    call check_condition(status%code == 0 .and. abs(lhs - rhs) < 1.0e-12_dp, &
        "Fourier vector-product VJP satisfies the complex dot identity")

    call assemble_fourier_vector_product( &
        registry, coupling, left, right, bad_product, status)
    call check_condition(status%code /= 0, &
        "Fourier vector product rejects an incompatible output shape")
    call check_summary("Fourier vector product")

contains

    pure integer function find_mode(mode_set, requested) result(index)
        integer, intent(in) :: mode_set(:), requested
        integer :: candidate

        index = 0
        do candidate = 1, size(mode_set)
            if (mode_set(candidate) /= requested) cycle
            index = candidate
            return
        end do
    end function find_mode

end program test_fourier_vector_product
