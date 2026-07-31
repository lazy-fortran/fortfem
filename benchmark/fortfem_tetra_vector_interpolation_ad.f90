program fortfem_tetra_vector_interpolation_ad
    use, intrinsic :: iso_fortran_env, only: int64
    use fortfem_api, only: initialize_tetra_nedelec_first_kind, &
        initialize_tetra_rt, interpolate_sampled_physical_tetra_nedelec, &
        interpolate_sampled_physical_tetra_nedelec_jvp, &
        interpolate_sampled_physical_tetra_nedelec_vjp, &
        interpolate_sampled_physical_tetra_rt, &
        interpolate_sampled_physical_tetra_rt_jvp, &
        interpolate_sampled_physical_tetra_rt_vjp, &
        tetra_nedelec_first_kind_t, tetra_nedelec_interpolation_points, &
        tetra_rt_interpolation_points, tetra_rt_t, &
        tetra_vector_sample_gradients_t, tetra_vector_samples_t, &
        zero_tetra_vector_samples_like
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: repetitions = 200
    type(tetra_nedelec_first_kind_t) :: nedelec
    type(tetra_rt_t) :: rt
    type(tetra_vector_samples_t) :: nedelec_points, nedelec_samples
    type(tetra_vector_samples_t) :: nedelec_samples_bar
    type(tetra_vector_samples_t) :: nedelec_parameter_dot
    type(tetra_vector_samples_t) :: rt_points, rt_samples, rt_samples_bar
    type(tetra_vector_samples_t) :: rt_parameter_dot
    type(tetra_vector_sample_gradients_t) :: nedelec_gradients, rt_gradients
    real(dp), allocatable :: dofs(:), dofs_bar(:), dofs_dot(:)
    real(dp) :: checksum, vertices(3, 4), vertices_bar(3, 4)
    real(dp) :: vertices_dot(3, 4)
    integer :: iteration, status
    integer(int64) :: finish, rate, start

    vertices = reshape([ &
        0.1_dp, -0.2_dp, 0.05_dp, 1.3_dp, 0.1_dp, -0.1_dp, &
        -0.15_dp, 1.2_dp, 0.2_dp, 0.2_dp, -0.1_dp, 1.4_dp], [3, 4])
    vertices_dot = 0.01_dp
    checksum = 0.0_dp
    call system_clock(count_rate=rate)

    call initialize_tetra_nedelec_first_kind(3, nedelec, status)
    call require_success("initialize Nedelec element", status)
    call tetra_nedelec_interpolation_points( &
        nedelec, vertices, nedelec_points, status)
    call require_success("create Nedelec interpolation points", status)
    call initialize_samples( &
        nedelec_points, nedelec_samples, nedelec_gradients, &
        nedelec_parameter_dot)
    call interpolate_sampled_physical_tetra_nedelec( &
        nedelec, vertices, nedelec_samples, dofs, status)
    call require_success("warm Nedelec interpolation", status)
    allocate(dofs_dot(size(dofs)), dofs_bar(size(dofs)), source=0.01_dp)
    call system_clock(start)
    do iteration = 1, repetitions
        call interpolate_sampled_physical_tetra_nedelec( &
            nedelec, vertices, nedelec_samples, dofs, status)
        checksum = checksum + dofs(1)
    end do
    call require_success("benchmark Nedelec interpolation", status)
    call system_clock(finish)
    call report("nedelec_primal", start, finish, rate)
    call system_clock(start)
    do iteration = 1, repetitions
        call interpolate_sampled_physical_tetra_nedelec_jvp( &
            nedelec, vertices, nedelec_samples, nedelec_gradients, &
            vertices_dot, nedelec_parameter_dot, dofs_dot, status)
        checksum = checksum + dofs_dot(1)
    end do
    call require_success("benchmark Nedelec interpolation JVP", status)
    call system_clock(finish)
    call report("nedelec_jvp", start, finish, rate)
    call system_clock(start)
    do iteration = 1, repetitions
        call interpolate_sampled_physical_tetra_nedelec_vjp( &
            nedelec, vertices, nedelec_samples, nedelec_gradients, dofs_bar, &
            vertices_bar, nedelec_samples_bar, status)
        checksum = checksum + vertices_bar(1, 1)
    end do
    call require_success("benchmark Nedelec interpolation VJP", status)
    call system_clock(finish)
    call report("nedelec_vjp", start, finish, rate)

    deallocate(dofs, dofs_dot, dofs_bar)
    call initialize_tetra_rt(2, rt, status)
    call require_success("initialize RT element", status)
    call tetra_rt_interpolation_points(rt, vertices, rt_points, status)
    call require_success("create RT interpolation points", status)
    call initialize_samples( &
        rt_points, rt_samples, rt_gradients, rt_parameter_dot)
    call interpolate_sampled_physical_tetra_rt( &
        rt, vertices, rt_samples, dofs, status)
    call require_success("warm RT interpolation", status)
    allocate(dofs_dot(size(dofs)), dofs_bar(size(dofs)), source=0.01_dp)
    call system_clock(start)
    do iteration = 1, repetitions
        call interpolate_sampled_physical_tetra_rt( &
            rt, vertices, rt_samples, dofs, status)
        checksum = checksum + dofs(1)
    end do
    call require_success("benchmark RT interpolation", status)
    call system_clock(finish)
    call report("rt_primal", start, finish, rate)
    call system_clock(start)
    do iteration = 1, repetitions
        call interpolate_sampled_physical_tetra_rt_jvp( &
            rt, vertices, rt_samples, rt_gradients, vertices_dot, &
            rt_parameter_dot, dofs_dot, status)
        checksum = checksum + dofs_dot(1)
    end do
    call require_success("benchmark RT interpolation JVP", status)
    call system_clock(finish)
    call report("rt_jvp", start, finish, rate)
    call system_clock(start)
    do iteration = 1, repetitions
        call interpolate_sampled_physical_tetra_rt_vjp( &
            rt, vertices, rt_samples, rt_gradients, dofs_bar, vertices_bar, &
            rt_samples_bar, status)
        checksum = checksum + vertices_bar(1, 1)
    end do
    call require_success("benchmark RT interpolation VJP", status)
    call system_clock(finish)
    call report("rt_vjp", start, finish, rate)
    write(*, "(a,es12.4)") "checksum,", checksum

contains

    subroutine require_success(operation, operation_status)
        character(*), intent(in) :: operation
        integer, intent(in) :: operation_status

        if (operation_status /= 0) then
            write(*, "(a,a,i0)") trim(operation), " failed with status ", &
                operation_status
            error stop 1
        end if
    end subroutine require_success

    subroutine initialize_samples(points, samples, gradients, parameter_dot)
        type(tetra_vector_samples_t), intent(in) :: points
        type(tetra_vector_samples_t), intent(out) :: samples, parameter_dot
        type(tetra_vector_sample_gradients_t), intent(out) :: gradients

        call zero_tetra_vector_samples_like(points, samples)
        call zero_tetra_vector_samples_like(points, parameter_dot)
        allocate(gradients%edge(3, 3, size(points%edge, 2), &
            size(points%edge, 3)), source=0.0_dp)
        allocate(gradients%face(3, 3, size(points%face, 2), &
            size(points%face, 3)), source=0.0_dp)
        allocate(gradients%volume(3, 3, size(points%volume, 2)), source=0.0_dp)
        samples%edge = 0.25_dp
        samples%face = 0.25_dp
        samples%volume = 0.25_dp
        parameter_dot%edge = 0.01_dp
        parameter_dot%face = 0.01_dp
        parameter_dot%volume = 0.01_dp
    end subroutine initialize_samples

    subroutine report(label, first, last, clock_rate)
        character(*), intent(in) :: label
        integer(int64), intent(in) :: first, last, clock_rate
        real(dp) :: microseconds

        microseconds = 1.0e6_dp*real(last - first, dp)/ &
            real(clock_rate*repetitions, dp)
        write(*, "(a,a,f12.3)") label, ",microseconds_per_call,", microseconds
    end subroutine report

end program fortfem_tetra_vector_interpolation_ad
