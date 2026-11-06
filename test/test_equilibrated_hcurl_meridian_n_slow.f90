program test_equilibrated_hcurl_meridian_n_slow
    !! Full refinement ladders of the per-n meridian H(curl) equilibrated
    !! bound (meridian_hcurl_cases); test_equilibrated_hcurl_meridian_n runs
    !! the coarse levels.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_summary
    use fortfem_equilibrated_hcurl_meridian_n, only: HCURL_GLOBAL, HCURL_PATCHWISE
    use fortfem_kinds, only: dp
    use meridian_hcurl_cases, only: cases_passed, exact_recovery, run_manufactured
    use meridian_test_fields, only: ann_c, ann_r1, ann_r2, rect_chi
    implicit none

    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    ann_c = 3.0_dp
    ann_r1 = 0.4_dp
    ann_r2 = 1.0_dp
    call exact_recovery(3, 3, 2, 0.0_dp)
    call exact_recovery(3, 4, 1, 0.0_dp)
    call exact_recovery(3, 3, 2, 1.7_dp)
    rect_chi = .true.
    call run_manufactured(1, 2, [1, 3], 5, .true., HCURL_GLOBAL, 2.6_dp, .true.)
    call run_manufactured(1, 3, [2], 5, .true., HCURL_GLOBAL, 1.05_dp, .true.)
    call run_manufactured(2, 3, [1, 3], 5, .true., HCURL_GLOBAL, 1.8_dp, .true.)
    call run_manufactured(2, 4, [2], 4, .true., HCURL_GLOBAL, 1.05_dp, .true.)
    call run_manufactured(3, 4, [2], 3, .true., HCURL_GLOBAL, 1.3_dp, .true.)
    call run_manufactured(2, 3, [3], 4, .false., HCURL_GLOBAL, 0.0_dp, .false.)
    call run_manufactured(2, 3, [1, 3, 10], 4, .true., HCURL_PATCHWISE, 0.0_dp, &
        .false.)
    call check_summary("Meridian per-n H(curl) equilibrated bound, full ladders")
    if (.not. cases_passed) error stop 1
end program test_equilibrated_hcurl_meridian_n_slow
