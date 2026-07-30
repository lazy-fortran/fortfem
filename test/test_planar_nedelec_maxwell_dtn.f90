program test_planar_nedelec_maxwell_dtn
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_planar_nedelec_maxwell_dtn_form, &
        build_planar_nedelec_trace_sampling, &
        generate_structured_tetra_box_mesh, solve_tetra_nedelec_curl_mass
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: nx = 5, ny = 5
    real(dp), parameter :: wave_number = 3.7_dp
    real(dp), parameter :: exact_field(3) = [1.25_dp, -0.75_dp, 0.5_dp]
    real(dp), parameter :: bounds(3, 2) = reshape( &
        [0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 1.5_dp, 1.0_dp], [3, 2])
    real(dp), parameter :: origin(3) = [0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: periods(3, 2) = reshape( &
        [2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.5_dp, 0.0_dp], [3, 2])
    complex(dp), allocatable :: form(:, :), sampling(:, :), trace(:)
    complex(dp), allocatable :: weak_expected(:), weak_result(:)
    integer, allocatable :: boundary_dofs(:), tetrahedra(:, :)
    real(dp), allocatable :: coefficients(:), vertices(:, :)
    type(fortsparse_status_t) :: sparse_status
    real(dp) :: weight
    integer :: order, status
    logical :: all_passed

    all_passed = .true.
    call generate_structured_tetra_box_mesh( &
        bounds, [1, 1, 1], vertices, tetrahedra, status)
    if (status /= 0) error stop "structured tetrahedral mesh failed"
    weight = norm2(periods(:, 1))*norm2(periods(:, 2))/real(nx*ny, dp)

    do order = 1, 4
        call solve_tetra_nedelec_curl_mass( &
            vertices, tetrahedra, order, constant_source, 1.0_dp, 1.0_dp, &
            coefficients, sparse_status)
        if (sparse_status%code /= 0) error stop "constant Hcurl solve failed"
        call build_planar_nedelec_trace_sampling( &
            vertices, tetrahedra, order, origin, periods, nx, ny, &
            boundary_dofs, sampling, status)
        if (status /= 0) error stop "Nedelec trace sampling failed"
        trace = matmul( &
            sampling, cmplx(coefficients(boundary_dofs), 0.0_dp, dp))
        call record_condition( &
            maxval(abs(trace(1::2) - exact_field(1))) < 2.0e-10_dp .and. &
            maxval(abs(trace(2::2) - exact_field(2))) < 2.0e-10_dp, &
            "arbitrary-order Nedelec sampling reproduces constant trace")

        call assemble_planar_nedelec_maxwell_dtn_form( &
            vertices, tetrahedra, order, origin, periods, nx, ny, &
            wave_number, boundary_dofs, form, status)
        if (status /= 0) error stop "Nedelec Maxwell DtN assembly failed"
        allocate(weak_expected(size(boundary_dofs)))
        weak_expected = matmul( &
            transpose(sampling), &
            -cmplx(0.0_dp, wave_number*weight, dp)*trace)
        weak_result = matmul( &
            form, cmplx(coefficients(boundary_dofs), 0.0_dp, dp))
        call record_condition( &
            maxval(abs(weak_result - weak_expected)) < 3.0e-8_dp, &
            "pulled-back DtN matches exact normal-incidence capacity map")
        call record_condition( &
            maxval(abs(form - transpose(form))) < 3.0e-10_dp, &
            "Nedelec Maxwell DtN block is complex symmetric")
        deallocate( &
            coefficients, boundary_dofs, sampling, trace, form, &
            weak_expected, weak_result)
    end do
    call check_summary("Planar arbitrary-order Nedelec Maxwell DtN")
    if (.not. all_passed) error stop 1

contains

    pure subroutine constant_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value(3)

        associate(unused => x + y + z)
        end associate
        value = exact_field
    end subroutine constant_source

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_planar_nedelec_maxwell_dtn
