program test_fourier_mode_properties
    use, intrinsic :: iso_fortran_env, only: int32
    use check, only: check_condition, check_property, check_summary, &
        property_random_integer, property_random_unit, property_rng_initialize, &
        property_rng_t
    use fortfem_fourier, only: &
        build_fourier_mode_padded_registry, find_fourier_mode, &
        fourier_mode_registry_t, initialize_fourier_mode_registry, &
        validate_fourier_mode_registry
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    logical :: all_passed
    integer(int32) :: first_failed_seed, shrunk_seed

    call check_property( &
        "random conjugate-packed Fourier sets preserve one-product closure", &
        314159_int32, 20, mode_case, all_passed, first_failed_seed, shrunk_seed)
    call check_condition(all_passed .and. first_failed_seed == 0_int32 .and. &
        shrunk_seed == 0_int32, &
        "random Fourier mode property reports no failure seed")
    call check_summary("random Fourier mode registry properties")
    if (.not. all_passed) error stop 1

contains

    logical function mode_case(case_seed)
        integer(int32), intent(in) :: case_seed
        integer, parameter :: maximum_pairs = 4
        type(property_rng_t) :: rng
        type(fourier_mode_registry_t) :: registry, padded
        type(fortsparse_status_t) :: status
        integer, allocatable :: base_poloidal(:), base_toroidal(:)
        integer, allocatable :: poloidal(:), toroidal(:), radial(:)
        real(dp), allocatable :: normalization(:)
        integer :: pair_count, pair, attempt, first_mode, second_mode
        integer :: candidate_m, candidate_n, mode_count
        integer :: mode_index, field_periods
        real(dp) :: poloidal_phase, toroidal_phase
        logical :: accepted

        mode_case = .false.
        call property_rng_initialize(rng, case_seed)
        pair_count = property_random_integer(rng, 1, maximum_pairs)
        allocate(base_poloidal(pair_count), base_toroidal(pair_count))
        do pair = 1, pair_count
            accepted = .false.
            do attempt = 1, 100
                candidate_m = property_random_integer(rng, -3, 3)
                candidate_n = property_random_integer(rng, -3, 3)
                if (candidate_m == 0 .and. candidate_n == 0) cycle
                accepted = .true.
                do first_mode = 1, pair - 1
                    if (candidate_m == base_poloidal(first_mode) .and. &
                        candidate_n == base_toroidal(first_mode)) accepted = .false.
                    if (candidate_m == -base_poloidal(first_mode) .and. &
                        candidate_n == -base_toroidal(first_mode)) accepted = .false.
                end do
                if (accepted) exit
            end do
            if (.not. accepted) return
            base_poloidal(pair) = candidate_m
            base_toroidal(pair) = candidate_n
        end do

        mode_count = 2*pair_count
        allocate(poloidal(mode_count), toroidal(mode_count), radial(mode_count), &
            normalization(mode_count))
        do pair = 1, pair_count
            poloidal(2*pair - 1) = base_poloidal(pair)
            toroidal(2*pair - 1) = base_toroidal(pair)
            poloidal(2*pair) = -base_poloidal(pair)
            toroidal(2*pair) = -base_toroidal(pair)
            radial(2*pair - 1:2*pair) = abs(base_poloidal(pair))
            normalization(2*pair - 1:2*pair) = 0.75_dp + &
                0.5_dp*property_random_unit(rng)
        end do
        field_periods = property_random_integer(rng, 1, 4)
        poloidal_phase = property_random_unit(rng) - 0.5_dp
        toroidal_phase = property_random_unit(rng) - 0.5_dp
        call initialize_fourier_mode_registry( &
            registry, poloidal, toroidal, field_periods, poloidal_phase, &
            toroidal_phase, .true., radial, normalization, status)
        if (status%code /= 0) return
        if (.not. validate_fourier_mode_registry(registry, status)) return
        if (status%code /= 0) return
        call build_fourier_mode_padded_registry(registry, padded, status)
        if (status%code /= 0) return

        do first_mode = 1, mode_count
            do second_mode = 1, mode_count
                mode_index = find_fourier_mode(padded, &
                    poloidal(first_mode) + poloidal(second_mode), &
                    toroidal(first_mode) + toroidal(second_mode))
                if (mode_index == 0) return
            end do
        end do
        mode_case = .true.
    end function mode_case

end program test_fourier_mode_properties
