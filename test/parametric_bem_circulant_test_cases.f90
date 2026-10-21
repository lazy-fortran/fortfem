module parametric_bem_circulant_test_cases
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use check, only: check_condition, check_summary
    use fortfem_torus_surface_mesh, only: generate_torus_surface_mesh
    use fortfem_parametric_laplace_bem_3d, only: bem_operators_t
    use fortfem_parametric_bem_high_order, only: assemble_parametric_bem_p2, &
             assemble_parametric_bem_p2_circulant, bem_pair_options_t, bem_pair_report_t
    use fortfem_bem_circulant, only: bem_circulant_t, circulant_matrix_t, &
                                     assignment(=)
    use fortnum_linalg, only: dense_solve
    use parametric_bem_test_geometry, only: test_torus_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads, omp_set_num_threads
    implicit none
    private
    public :: run_circulant_case
contains
    subroutine run_circulant_case(nphi, scaling)
        integer, intent(in) :: nphi
        logical, intent(in) :: scaling
        type(test_torus_t) :: torus
        type(bem_operators_t) :: dense, expanded
        type(bem_circulant_t) :: circular
        type(bem_pair_options_t) :: options
        type(bem_pair_report_t) :: report, full_report
        type(circulant_matrix_t) :: jump
        integer, allocatable :: td(:, :), fd(:, :)
        real(dp), allocatable :: vertices(:, :), u(:), rhs(:), flux(:), exact(:), y(:)
        real(dp) :: start, full_time, circular_time, difference, x(3), a(3), b(3)
        real(dp), parameter :: nodes(2, 6) = reshape( &
                               [0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 1._dp, &
                                .5_dp, 0._dp, .5_dp, .5_dp, 0._dp, .5_dp], [2, 6])
        integer :: status, i, j, threads, previous_threads
        real(dp) :: serial_time
        logical :: passed
        passed = .true.
        call generate_torus_surface_mesh(2._dp, .6_dp, 3, nphi, vertices, &
                                         torus%faces, torus%parameters)
        options%relative_tolerance = .03_dp
        options%absolute_tolerance = 1.e-7_dp
        options%maximum_degree = 48
        start = omp_get_wtime()
        call assemble_parametric_bem_p2(torus, torus%faces, options, dense, &
                                        td, fd, full_report, status)
        full_time = omp_get_wtime() - start
        call record(status == 0, 'all-pairs torus assembly converges')
        if (status /= 0) error stop 1
        start = omp_get_wtime()
        call assemble_parametric_bem_p2_circulant(torus, torus%faces, nphi, options, &
                                                  circular, td, fd, report, status)
        circular_time = omp_get_wtime() - start
        call record(status == 0, 'sector assembly converges')
        if (status /= 0) error stop 1
        call circular%expand(expanded, status)
        call record(status == 0, 'circulant operator expands in original DOF order')
        difference = max(relative(dense%v, expanded%v), relative(dense%k, expanded%k), &
                     relative(dense%w, expanded%w), relative(dense%mass, expanded%mass))
        write (*, '(a,es16.8)') 'maximum relative matrix difference ', difference
        call record(difference < 1.e-10_dp, 'all four expanded forms equal all-pairs')
        call record(maxval(abs(dense%area - expanded%area)) < 1.e-11_dp, &
                    'panel areas equal all-pairs')
        call record(report%accepted_pairs*nphi == full_report%accepted_pairs, &
                    'exactly one sector of target pairs is evaluated')
        write (*, '(a,2i5,3f14.6)') 'threads sectors dense circulant ratio ', &
          omp_get_max_threads(), nphi, full_time, circular_time, full_time/circular_time
        write (*, '(a,2i18)') 'kernel evaluations dense circulant ', &
            full_report%kernel_evaluations, report%kernel_evaluations
        allocate (u(size(dense%w, 1)), rhs(size(dense%v, 1)), flux(size(dense%v, 1)), &
                  exact(size(dense%v, 1)), y(size(dense%v, 1)))
        do i = 1, size(torus%faces, 2)
            do j = 1, 6
                call torus%evaluate(i, nodes(:, j), x, a, b, status)
                u(td(j, i)) = x(3)
            end do
        end do
        ! u=z is an exact interior harmonic, with flux n_z. Its nonconstant
        ! trace checks seam DOFs as well as the zero sector harmonic.
        jump = circular%k
        jump%blocks = circular%k%blocks + .5_dp*circular%mass%blocks
        expanded%k = dense%k + .5_dp*dense%mass
        rhs = matmul(expanded%k, u)
        call dense_solve(dense%v, rhs, exact, status)
        call record(status == 0, 'dense interior harmonic solve succeeds')
        call jump%apply(u, y, status)
        call record(status == 0, 'rectangular block apply succeeds')
        call record(norm2(y - rhs) < 1.e-10_dp, 'block apply equals dense apply')
        call circular%v%solve(rhs, flux, status)
        call record(status == 0, 'per-harmonic interior solve succeeds')
        call record(norm2(flux - exact)/max(1._dp, norm2(exact)) < 1.e-10_dp, &
                    'per-harmonic interior flux equals dense solve')
        ! The trace of a point source inside the solid torus is exterior harmonic.
        do i = 1, size(torus%faces, 2)
            do j = 1, 6
                call torus%evaluate(i, nodes(:, j), x, a, b, status)
                u(td(j, i)) = 1._dp/norm2(x - [2._dp, 0._dp, 0._dp])
            end do
        end do
        expanded%k = dense%k - .5_dp*dense%mass
        rhs = matmul(expanded%k, u)
        call dense_solve(dense%v, rhs, exact, status)
        call record(status == 0, 'dense exterior point-source solve succeeds')
        call circular%v%solve(rhs, flux, status)
        call record(status == 0, 'per-harmonic exterior solve succeeds')
        call record(norm2(flux - exact)/max(1._dp, norm2(exact)) < 1.e-10_dp, &
                    'per-harmonic exterior flux equals dense solve for nonzero modes')
        if (scaling) then
            previous_threads = omp_get_max_threads()
            threads = 1
            do i = 1, 4
                call omp_set_num_threads(threads)
                start = omp_get_wtime()
                call assemble_parametric_bem_p2_circulant(torus, torus%faces, nphi, &
                                              options, circular, td, fd, report, status)
                circular_time = omp_get_wtime() - start
                call record(status == 0, 'thread sweep assembly converges')
                if (status /= 0) error stop 1
                if (threads == 1) serial_time = circular_time
                call circular%expand(expanded, status)
                difference = max(relative(dense%v, expanded%v), &
                           relative(dense%k, expanded%k), relative(dense%w, expanded%w))
                call record(difference < 1.e-10_dp, 'thread count preserves all forms')
            write (*, '(a,i3,3es16.8)') 'scaling threads seconds speedup difference ', &
                    threads, circular_time, serial_time/circular_time, difference
                threads = 2*threads
            end do
            call omp_set_num_threads(previous_threads)
        end if
        call assemble_parametric_bem_p2_circulant(torus, torus%faces, 0, options, &
                                                  circular, td, fd, report, status)
        call record(status == 1, 'invalid sector count is rejected')
        torus%deformation(1, 1) = 1.1_dp
        call assemble_parametric_bem_p2_circulant(torus, torus%faces, nphi, options, &
                                                  circular, td, fd, report, status)
        call record(status == 1, 'nonrotational geometry is rejected')
        call check_summary('Block-circulant P2 BEM')
        if (.not. passed) error stop 1
    contains
        real(dp) function relative(a, b)
            real(dp), intent(in) :: a(:, :), b(:, :)
            relative = maxval(abs(a - b))/max(tiny(1._dp), maxval(abs(a)))
        end function
        subroutine record(ok, label)
            logical, intent(in) :: ok
            character(*), intent(in) :: label
            passed = passed .and. ok
            call check_condition(ok, label)
        end subroutine
    end subroutine
end module
