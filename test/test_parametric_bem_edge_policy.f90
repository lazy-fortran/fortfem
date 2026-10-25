module parametric_bem_edge_policy_test_support
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    use fortfem_parametric_bem_high_order, only: bem_shared_edge_chart_t
    use parametric_bem_test_geometry, only: test_sphere_t
    implicit none
    type, extends(bem_surface_t) :: pair_failure_surface_t
        type(test_sphere_t) :: sphere
    contains
        procedure :: evaluate => pair_failure_evaluate
    end type
contains
    pure subroutine pair_failure_evaluate(self, panel, uv, x, a, b, status)
        class(pair_failure_surface_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        call self%sphere%evaluate(panel, uv, x, a, b, status)
        ! Gauss-Duffy precomputation avoids the centroid. Native self/closest
        ! integration uses it, so status 1 originates in the pair integrator.
        if (maxval(abs(uv-1._dp/3)) < 1.e-14_dp) status = 1
    end subroutine

    pure subroutine inactive(surface, i, j, target, source, chart, selected, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i, j, target(3), source(3)
        type(bem_shared_edge_chart_t), intent(out) :: chart
        logical, intent(out) :: selected
        integer, intent(out) :: status
        real(dp) :: x(3), a(3), b(3)
        chart = bem_shared_edge_chart_t()
        selected = .false.
        call surface%evaluate(i, [1._dp/3, 1._dp/3], x, a, b, status)
        if (i == j .or. minval(target) < 1 .or. minval(source) < 1) status = 1
    end subroutine

    pure subroutine select_edge(surface, i, j, target, source, chart, &
                               selected, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i, j, target(3), source(3)
        type(bem_shared_edge_chart_t), intent(out) :: chart
        logical, intent(out) :: selected
        integer, intent(out) :: status
        integer :: shared(3), n, a, endpoint
        call inactive(surface, i, j, target, source, chart, selected, status)
        if (status /= 0) return
        n = 0
        do a = 1, 3
            if (.not. any(source == target(a))) cycle
            n = n + 1
            shared(n) = target(a)
        end do
        if (n /= 2) return
        shared(1:2) = [minval(shared(:2)), maxval(shared(:2))]
        do endpoint = 1, 2
            chart%vertices(endpoint, 1) = findloc(target, shared(endpoint), dim=1)
            chart%vertices(endpoint, 2) = findloc(source, shared(endpoint), dim=1)
        end do
        chart%vertices(3, :) = 6 - sum(chart%vertices(:2, :), dim=1)
        selected = .true.
    end subroutine

    pure subroutine bad_aspect(surface, i, j, target, source, chart, &
                              selected, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i, j, target(3), source(3)
        type(bem_shared_edge_chart_t), intent(out) :: chart
        logical, intent(out) :: selected
        integer, intent(out) :: status
        call select_edge(surface, i, j, target, source, chart, selected, status)
        if (selected) chart%aspect = 0._dp
    end subroutine

    pure subroutine bad_vertices(surface, i, j, target, source, chart, &
                                selected, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i, j, target(3), source(3)
        type(bem_shared_edge_chart_t), intent(out) :: chart
        logical, intent(out) :: selected
        integer, intent(out) :: status
        call select_edge(surface, i, j, target, source, chart, selected, status)
        if (selected) chart%vertices(1, 1) = 0
    end subroutine

    pure subroutine mismatched_endpoints(surface, i, j, target, source, chart, &
                                    selected, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i, j, target(3), source(3)
        type(bem_shared_edge_chart_t), intent(out) :: chart
        logical, intent(out) :: selected
        integer, intent(out) :: status
        call select_edge(surface, i, j, target, source, chart, selected, status)
        if (selected) chart%vertices(1:2, 2) = chart%vertices(2:1:-1, 2)
    end subroutine

    pure subroutine failing(surface, i, j, target, source, chart, selected, status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i, j, target(3), source(3)
        type(bem_shared_edge_chart_t), intent(out) :: chart
        logical, intent(out) :: selected
        integer, intent(out) :: status
        call inactive(surface, i, j, target, source, chart, selected, status)
        status = 2
    end subroutine
end module

program test_parametric_bem_edge_policy
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t
    use fortfem_parametric_bem_high_order, only: assemble_parametric_bem_p2, &
        integrate_parametric_bem_p2_pair, bem_pair_options_t, bem_pair_report_t, &
        bem_shared_edge_chart_t
    use parametric_bem_test_geometry, only: test_sphere_t
    use parametric_bem_edge_policy_test_support, only: inactive, select_edge, &
        bad_aspect, bad_vertices, mismatched_endpoints, failing, &
        pair_failure_surface_t
    implicit none
    type(test_sphere_t) :: sphere
    type(pair_failure_surface_t) :: failure_surface
    type(bem_operators_t) :: baseline, inactive_op, selected_op, failed_op
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: baseline_report, report, raw_report
    type(bem_shared_edge_chart_t) :: chart, swapped
    integer, allocatable :: td(:, :), fd(:, :)
    real(dp) :: v(3, 3), k(3, 6), w(6, 6), c(3, 6), reverse_v(3, 3)
    real(dp), allocatable :: ones(:)
    integer :: status, i, j, shared_pairs, nf
    logical :: passed, selected
    passed = .true.
    call generate_sphere_surface_mesh(1._dp, 0, sphere%vertices, sphere%faces)
    nf = size(sphere%faces, 2)
    ! This is a dispatch/reporting test, not a new accuracy tolerance.
    options%relative_tolerance = .1_dp
    options%absolute_tolerance = 1.e-7_dp
    options%initial_degree = 8
    options%maximum_degree = 24
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, baseline, &
                                  td, fd, baseline_report, status)
    call record(status == 0, 'existing no-policy assembly succeeds')
    if (status /= 0) error stop 1
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, inactive_op, &
                                  td, fd, report, status, edge_selector=inactive)
    call record(status == 0, 'inactive policy preserves all non-self paths')
    if (status /= 0) error stop 1
    call record(maxval(abs(inactive_op%v-baseline%v)) < 1.e-13_dp .and. &
                maxval(abs(inactive_op%k-baseline%k)) < 1.e-13_dp .and. &
                maxval(abs(inactive_op%w-baseline%w)) < 1.e-12_dp, &
                'inactive selector preserves V K W results')
    call record(report%kernel_evaluations == baseline_report%kernel_evaluations &
                .and. report%maximum_degree_used == &
                baseline_report%maximum_degree_used .and. &
                report%error_ratio == baseline_report%error_ratio, &
                'inactive selector preserves native work and estimator reports')
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, selected_op, &
                               td, fd, report, status, edge_selector=select_edge)
    call record(status == 0 .and. report%accepted_pairs == nf**2 .and. &
                report%failed_pairs == 0, 'selected global assembly reports all pairs')
    if (status /= 0) error stop 1
    ! Each selector rejects i==j, so successful assembly proves self bypass.
    call integrate_parametric_bem_p2_pair(sphere, 1, 1, options, &
                                         v, k, w, c, raw_report, status)
    call record(status == 0, 'native self rule remains available')
    if (status /= 0) error stop 1
    call record(maxval(abs(selected_op%v(fd(:, 1), fd(:, 1)) - &
                .5_dp*(v+transpose(v)))) < 1.e-13_dp, &
                'policy assembly keeps the original coincident self block')
    shared_pairs = 0
    do i = 1, nf
        do j = 1, nf
            if (i == j) cycle
            call select_edge(sphere, i, j, sphere%faces(:, i), &
                             sphere%faces(:, j), chart, selected, status)
            if (.not. selected) cycle
            shared_pairs = shared_pairs + 1
        end do
    end do
    outer: do i = 1, nf
        do j = i + 1, nf
            call select_edge(sphere, i, j, sphere%faces(:, i), &
                             sphere%faces(:, j), chart, selected, status)
            if (selected) exit outer
        end do
    end do outer
    call record(selected, 'closed sphere fixture supplies a shared-edge pair')
    if (.not. selected) error stop 1
    call integrate_parametric_bem_p2_pair(sphere, i, j, options, &
                                v, k, w, c, raw_report, status, edge_chart=chart)
    call record(status == 0, 'selected forward chart succeeds in the raw API')
    if (status /= 0) error stop 1
    swapped = chart
    swapped%vertices = chart%vertices(:, 2:1:-1)
    call integrate_parametric_bem_p2_pair(sphere, j, i, options, reverse_v, &
                                 k, w, c, raw_report, status, edge_chart=swapped)
    call record(status == 0, 'swapped panels retain matching endpoint orientation')
    if (status /= 0) error stop 1
    call record(maxval(abs(selected_op%v(fd(:, i), fd(:, j)) - &
                .5_dp*(v+transpose(reverse_v)))) < 1.e-13_dp, &
                'global V dispatch matches both oriented raw chart blocks')
    allocate (ones(size(selected_op%k, 2)))
    ones = 1._dp
    call record(norm2(matmul(selected_op%k+.5_dp*selected_op%mass, ones)) &
                < 1.e-11_dp, 'selected charts preserve target correction convention')
    call record(norm2(matmul(selected_op%w, ones)) < 1.e-10_dp .and. &
                maxval(abs(selected_op%w-transpose(selected_op%w))) < 1.e-13_dp, &
                'selected charts preserve W symmetry and the constant nullspace')
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, failed_op, &
                                td, fd, report, status, edge_selector=bad_aspect)
    call record(status == 1 .and. report%failed_pairs == shared_pairs, &
                'invalid selected aspect fails every affected pair explicitly')
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, failed_op, &
                              td, fd, report, status, edge_selector=bad_vertices)
    call record(status == 1 .and. report%failed_pairs == shared_pairs, &
                'invalid local permutation is rejected before vertex indexing')
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, failed_op, &
                    td, fd, report, status, edge_selector=mismatched_endpoints)
    call record(status == 1 .and. report%failed_pairs == shared_pairs, &
                'mismatched shared endpoints fail rather than silently falling back')
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, failed_op, &
                                   td, fd, report, status, edge_selector=failing)
    call record(status == 2 .and. report%failed_pairs == nf*(nf-1) .and. &
                report%accepted_pairs == nf, 'callback failure preserves pair counts')
    call generate_sphere_surface_mesh(1._dp, 0, &
                       failure_surface%sphere%vertices, failure_surface%sphere%faces)
    call assemble_parametric_bem_p2(failure_surface, failure_surface%sphere%faces, &
                                   options, failed_op, td, fd, report, status)
    call record(status == 3 .and. report%failed_pairs > 0 .and. &
                report%maximum_degree_used > 0, &
                'no-policy native status 1 retains legacy aggregate status 3')
    options%maximum_evaluations = 1
    call assemble_parametric_bem_p2(sphere, sphere%faces, options, failed_op, &
                               td, fd, report, status, edge_selector=select_edge)
    call record(status == 3 .and. report%failed_pairs == nf**2 .and. &
                report%accepted_pairs == 0, 'native quadrature cap remains a failure')
    call check_summary('Optional BEM shared-edge assembly policy')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
