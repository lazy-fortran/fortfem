program tetra_h1_poisson
    use fortfem_api, only: evaluate_tetra_lagrange_solution, &
        evaluate_tetra_lagrange_solution_at_point, &
        solve_tetra_lagrange_poisson, tetra_duffy_quadrature
    use fortfem_generated_tetra_h1_oracle, only: generated_tetra_h1_oracle
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortplot, only: add_3d_plot, add_scatter, figure, legend, plot, &
        savefig, set_yscale, title, xlabel, ylabel
    use fortsparse, only: fortsparse_status_t
    implicit none

    character(len=*), parameter :: output_directory = &
        "output/example/tetra_h1_poisson"
    integer, parameter :: tetrahedra(4, 4) = reshape( &
        [1, 2, 3, 5, 1, 4, 2, 5, 1, 3, 4, 5, 2, 4, 3, 5], [4, 4])
    real(dp), parameter :: vertices(3, 5) = reshape( &
        [0.0_dp, 0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, 0.0_dp, &
        0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, &
        0.25_dp, 0.25_dp, 0.25_dp], [3, 5])

    real(dp) :: degrees(4), h1_errors(4), l2_errors(4)
    real(dp), allocatable :: plot_solution(:)
    type(fortsparse_status_t) :: plot_status
    integer :: command_status, degree, unit

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"
    call initialize_gallery_sequence()
    do degree = 1, 4
        degrees(degree) = real(degree, dp)
        call solve_and_measure(degree, l2_errors(degree), h1_errors(degree))
    end do
    if (any(l2_errors(2:) >= l2_errors(:3))) then
        error stop "tetrahedral H1 L2 errors did not decrease"
    end if
    if (any(h1_errors(2:) >= h1_errors(:3))) then
        error stop "tetrahedral H1 gradient errors did not decrease"
    end if
    if (l2_errors(4) >= 1.0e-12_dp) then
        error stop "quartic tetrahedral H1 value was not reproduced"
    end if
    if (h1_errors(4) >= 1.0e-11_dp) then
        error stop "quartic tetrahedral H1 gradient was not reproduced"
    end if

    call solve_tetra_lagrange_poisson( &
        vertices, tetrahedra, 4, bubble_source, zero_boundary, &
        plot_solution, plot_status)
    if (plot_status%code /= 0) error stop "tetrahedral plot solve failed"
    call render_solution(plot_solution)
    call record_gallery_stage("physical_solution")

    open (newunit=unit, file=output_directory//"/convergence.csv", &
        status="replace", action="write")
    write (unit, "(a)") "degree,l2_error,gradient_error"
    do degree = 1, 4
        write (unit, "(i0,2(',',es24.16))") &
            degree, l2_errors(degree), h1_errors(degree)
    end do
    close (unit)

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot( &
        degrees, l2_errors, label="L2 value error", linestyle="-", marker="o")
    call plot( &
        degrees, h1_errors, label="gradient error", linestyle="--", marker="s")
    call set_yscale("log")
    call xlabel("tetrahedral Lagrange degree")
    call ylabel("error norm")
    call title("FortSym oracle vs FortSparse tetrahedral H1 solve")
    call legend()
    call savefig(output_directory//"/tetra_h1_poisson_convergence.png")
    call record_gallery_stage("diagnostics")

contains

    subroutine initialize_gallery_sequence()
        integer :: sequence_unit

        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write")
        close (sequence_unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(len=*), intent(in) :: stage
        integer :: sequence_unit

        open (newunit=sequence_unit, &
            file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write")
        write (sequence_unit, "(a)") trim(stage)
        close (sequence_unit)
    end subroutine record_gallery_stage

    subroutine render_solution(solution)
        real(dp), intent(in) :: solution(:)

        integer, parameter :: sample_side = 24
        real(dp), allocatable :: x_plot(:), y_plot(:), z_plot(:), values(:)
        real(dp) :: point(3), value, gradient(3)
        integer :: cell, count, i, j, k, local_status

        allocate( &
            x_plot(4*sample_side**3), y_plot(4*sample_side**3), &
            z_plot(4*sample_side**3), values(4*sample_side**3))
        count = 0
        do cell = 1, size(tetrahedra, 2)
            do k = 0, sample_side
                do j = 0, sample_side
                    do i = 0, sample_side
                        if (i + j + k > sample_side) cycle
                        point = vertices(:, tetrahedra(1, cell)) + &
                            real(i, dp)/real(sample_side, dp)* &
                            (vertices(:, tetrahedra(2, cell)) - &
                            vertices(:, tetrahedra(1, cell))) + &
                            real(j, dp)/real(sample_side, dp)* &
                            (vertices(:, tetrahedra(3, cell)) - &
                            vertices(:, tetrahedra(1, cell))) + &
                            real(k, dp)/real(sample_side, dp)* &
                            (vertices(:, tetrahedra(4, cell)) - &
                            vertices(:, tetrahedra(1, cell)))
                        call evaluate_tetra_lagrange_solution_at_point( &
                            vertices, tetrahedra, 4, solution, cell, point, &
                            value, gradient, local_status)
                        if (local_status /= 0) cycle
                        count = count + 1
                        x_plot(count) = point(1)
                        y_plot(count) = point(2)
                        z_plot(count) = point(3)
                        values(count) = value
                    end do
                end do
            end do
        end do

        call figure(figsize=[7.5_dp, 6.0_dp])
        call add_scatter( &
            x_plot(:count), y_plot(:count), z_plot(:count), &
            c=values(:count), cmap="viridis", marker=".", &
            markersize=2.5_dp, label="computed quartic solution samples")
        call render_tetrahedron_edges()
        call title("Tetrahedral H1 Poisson computed solution")
        call savefig(output_directory//"/tetra_h1_poisson_solution_3d.png")
    end subroutine render_solution

    subroutine render_tetrahedron_edges()
        integer, parameter :: edge_count = 6
        integer, parameter :: edges(2, edge_count) = reshape([ &
            1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, edge_count])
        integer :: cell, edge

        do cell = 1, size(tetrahedra, 2)
            do edge = 1, edge_count
                call add_3d_plot( &
                    [vertices(1, tetrahedra(edges(1, edge), cell)), &
                    vertices(1, tetrahedra(edges(2, edge), cell))], &
                    [vertices(2, tetrahedra(edges(1, edge), cell)), &
                    vertices(2, tetrahedra(edges(2, edge), cell))], &
                    [vertices(3, tetrahedra(edges(1, edge), cell)), &
                    vertices(3, tetrahedra(edges(2, edge), cell))], &
                    color="dimgray", linewidth=1.1_dp)
            end do
        end do
    end subroutine render_tetrahedron_edges

    subroutine solve_and_measure(degree, l2_error, h1_error)
        integer, intent(in) :: degree
        real(dp), intent(out) :: l2_error, h1_error

        type(fortsparse_status_t) :: status
        real(dp), allocatable :: solution(:)

        call solve_tetra_lagrange_poisson( &
            vertices, tetrahedra, degree, bubble_source, zero_boundary, &
            solution, status)
        if (status%code /= 0) error stop "tetrahedral H1 solve failed"
        call measure_error(degree, solution, l2_error, h1_error)
    end subroutine solve_and_measure

    subroutine measure_error(degree, solution, l2_error, h1_error)
        integer, intent(in) :: degree
        real(dp), intent(in) :: solution(:)
        real(dp), intent(out) :: l2_error, h1_error

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: exact(5), gradient(3), jacobian(3, 3), point(3)
        real(dp) :: value, weight
        integer :: cell, node, quadrature_status, sample, status

        call tetra_duffy_quadrature(10, x, y, z, weights, quadrature_status)
        if (quadrature_status /= 0) error stop "error quadrature failed"
        l2_error = 0.0_dp
        h1_error = 0.0_dp
        do cell = 1, size(tetrahedra, 2)
            do node = 1, 3
                jacobian(:, node) = &
                    vertices(:, tetrahedra(node + 1, cell)) - &
                    vertices(:, tetrahedra(1, cell))
            end do
            do sample = 1, size(weights)
                point = vertices(:, tetrahedra(1, cell)) + &
                    matmul(jacobian, [x(sample), y(sample), z(sample)])
                call generated_tetra_h1_oracle( &
                    point(1), point(2), point(3), exact)
                call evaluate_tetra_lagrange_solution( &
                    vertices, tetrahedra, degree, solution, cell, &
                    [x(sample), y(sample), z(sample)], value, gradient, status)
                if (status /= 0) error stop "H1 field reconstruction failed"
                weight = det3(jacobian)*weights(sample)
                l2_error = l2_error + weight*(value - exact(1))**2
                h1_error = h1_error + weight*sum((gradient - exact(2:4))**2)
            end do
        end do
        l2_error = sqrt(l2_error)
        h1_error = sqrt(h1_error)
    end subroutine measure_error

    pure subroutine bubble_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        real(dp) :: oracle(5)

        call generated_tetra_h1_oracle(x, y, z, oracle)
        value = oracle(5)
    end subroutine bubble_source

    pure subroutine zero_boundary(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        value = 0.0_dp*(x + y + z)
    end subroutine zero_boundary

end program tetra_h1_poisson
