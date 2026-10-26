program test_parametric_bem_thread_determinism
    !! The P2 pair assembly computes pair blocks in parallel and accumulates
    !! them in a fixed order, so every operator is bitwise identical for any
    !! thread count (dense all-pairs and block-circulant paths).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t
    use fortfem_parametric_bem_high_order, only: assemble_parametric_bem_p2, &
        assemble_parametric_bem_p2_circulant, bem_pair_options_t, &
        bem_pair_report_t
    use fortfem_bem_circulant, only: bem_circulant_t
    use parametric_bem_test_geometry, only: test_torus_t
    use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    implicit none

    integer, parameter :: thread_counts(3) = [2, 3, 5]
    type(test_torus_t) :: torus
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    type(bem_circulant_t) :: reference_c, other_c
    type(bem_operators_t) :: reference_d, other_d
    integer, allocatable :: td(:, :), fd(:, :)
    real(dp), allocatable :: vertices(:, :)
    integer :: previous, n, status, nphi
    logical :: passed

    passed = .true.
    previous = omp_get_max_threads()
    options%relative_tolerance = .03_dp
    options%absolute_tolerance = 1.e-7_dp
    options%maximum_degree = 48

    nphi = 4
    call generate_torus_surface_mesh(2._dp, .6_dp, 3, nphi, vertices, &
        torus%faces, torus%parameters)
    call omp_set_num_threads(1)
    call assemble_parametric_bem_p2_circulant(torus, torus%faces, nphi, &
        options, reference_c, td, fd, report, status)
    call record(status == 0, 'serial circulant assembly converges')
    call assemble_parametric_bem_p2(torus, torus%faces, options, reference_d, &
        td, fd, report, status)
    call record(status == 0, 'serial all-pairs assembly converges')
    do n = 1, size(thread_counts)
        call omp_set_num_threads(thread_counts(n))
        call assemble_parametric_bem_p2_circulant(torus, torus%faces, nphi, &
            options, other_c, td, fd, report, status)
        call record(status == 0 .and. &
            all(other_c%v%blocks == reference_c%v%blocks) .and. &
            all(other_c%k%blocks == reference_c%k%blocks) .and. &
            all(other_c%w%blocks == reference_c%w%blocks) .and. &
            all(other_c%mass%blocks == reference_c%mass%blocks), &
            'circulant operators are bitwise independent of threads')
        call assemble_parametric_bem_p2(torus, torus%faces, options, other_d, &
            td, fd, report, status)
        call record(status == 0 .and. all(other_d%v == reference_d%v) .and. &
            all(other_d%k == reference_d%k) .and. &
            all(other_d%w == reference_d%w), &
            'all-pairs operators are bitwise independent of threads')
    end do
    call omp_set_num_threads(previous)
    call check_summary('Thread-deterministic P2 BEM assembly')
    if (.not. passed) error stop 1

contains

    subroutine record(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        call check_condition(condition, label)
        if (.not. condition) passed = .false.
    end subroutine record

end program test_parametric_bem_thread_determinism
