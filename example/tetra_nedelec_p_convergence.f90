program tetra_nedelec_p_convergence
    use fortfem_api, only: &
        evaluate_tetra_nedelec_first_kind, &
        initialize_tetra_nedelec_first_kind, &
        interpolate_reference_tetra_nedelec, tetra_duffy_quadrature, &
        tetra_nedelec_dof_count, tetra_nedelec_first_kind_t
    use fortfem_kinds, only: dp
    use fortplot, only: add_3d_plot, add_scatter, colorbar, figure, plot, &
        quiver, savefig, set_yscale, title, xlabel, ylabel
    implicit none

    character(*), parameter :: output_directory = &
        "output/example/tetra_nedelec_p_convergence"
    type(tetra_nedelec_first_kind_t) :: basis
    real(dp) :: degrees(4), errors(4)
    integer :: command_status, order, status

    call execute_command_line( &
        "mkdir -p "//output_directory, exitstat=command_status)
    if (command_status /= 0) error stop "cannot create example output directory"
    call initialize_gallery_sequence()

    print "(a)", "order  L2 error for manufactured curl-free field"
    do order = 1, 4
        degrees(order) = real(order, dp)
        call initialize_tetra_nedelec_first_kind(order, basis, status)
        if (status /= 0) error stop "basis initialization failed"
        call interpolation_error(basis, errors(order), status)
        if (status /= 0) error stop "interpolation failed"
        print "(i5,2x,es14.6)", order, errors(order)
    end do
    if (any(errors(2:4) >= errors(1:3))) then
        error stop "interpolation error did not decrease"
    end if
    if (errors(4) >= 2.0e-11_dp) then
        error stop "order four did not reproduce the cubic field"
    end if

    call render_field()
    call record_gallery_stage("physical_solution")

    call figure(figsize=[9.0_dp, 5.5_dp])
    call plot(degrees, errors, label="Nedelec interpolation error", &
        marker="o")
    call set_yscale("log")
    call xlabel("Nedelec polynomial order")
    call ylabel("L2 vector error")
    call title("Tetrahedral Nedelec p-convergence")
    call savefig(output_directory//"/p_convergence_1d.png")
    call record_gallery_stage("diagnostics")

contains

    subroutine initialize_gallery_sequence()
        integer :: unit, local_status

        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="replace", action="write", iostat=local_status)
        if (local_status /= 0) error stop "cannot initialize gallery sequence"
        close (unit)
    end subroutine initialize_gallery_sequence

    subroutine record_gallery_stage(stage)
        character(*), intent(in) :: stage
        integer :: unit, local_status

        open (newunit=unit, file=output_directory//"/gallery_sequence.txt", &
            status="old", position="append", action="write", &
            iostat=local_status)
        if (local_status /= 0) error stop "cannot record gallery sequence"
        write (unit, "(a)", iostat=local_status) stage
        close (unit)
        if (local_status /= 0) error stop "cannot write gallery sequence"
    end subroutine record_gallery_stage

    subroutine render_field()
        integer, parameter :: sample_side = 12
        real(dp), allocatable :: dofs(:), x_plot(:), y_plot(:), z_plot(:)
        real(dp), allocatable :: field_magnitude(:), basis_values(:, :)
        real(dp), allocatable :: basis_curls(:, :)
        real(dp) :: point(3), field_value(3), arrow_end(3), arrow_side(3)
        real(dp) :: field_scale, arrow_length
        integer :: count, i, j, k, local_status

        call interpolate_reference_tetra_nedelec( &
            basis, cubic_gradient, dofs, local_status)
        if (local_status /= 0) error stop "plot interpolation failed"
        allocate( &
            x_plot((sample_side + 1)**3), &
            y_plot((sample_side + 1)**3), &
            z_plot((sample_side + 1)**3), &
            field_magnitude((sample_side + 1)**3), &
            basis_values(3, size(dofs)), basis_curls(3, size(dofs)))
        count = 0
        do k = 0, sample_side
            do j = 0, sample_side
                do i = 0, sample_side
                    if (i + j + k > sample_side) cycle
                    point = [ &
                        real(i, dp)/real(sample_side, dp), &
                        real(j, dp)/real(sample_side, dp), &
                        real(k, dp)/real(sample_side, dp)]
                    call evaluate_tetra_nedelec_first_kind( &
                        basis, point, basis_values, basis_curls, local_status)
                    if (local_status /= 0) cycle
                    field_value = matmul(basis_values, dofs)
                    count = count + 1
                    x_plot(count) = point(1)
                    y_plot(count) = point(2)
                    z_plot(count) = point(3)
                    field_magnitude(count) = norm2(field_value)
                end do
            end do
        end do
        call figure(figsize=[7.5_dp, 6.0_dp])
        call add_scatter( &
            x_plot(:count), y_plot(:count), z_plot(:count), &
            c=field_magnitude(:count), cmap="viridis", marker=".", &
            markersize=4.0_dp, label="computed Nedelec field magnitude")
        call render_tetrahedron_edges()
        field_scale = maxval(field_magnitude(:count))
        field_scale = max(field_scale, epsilon(1.0_dp))
        do k = 0, sample_side, 3
            do j = 0, sample_side, 3
                do i = 0, sample_side, 3
                    if (i + j + k > sample_side) cycle
                    point = [ &
                        real(i, dp)/real(sample_side, dp), &
                        real(j, dp)/real(sample_side, dp), &
                        real(k, dp)/real(sample_side, dp)]
                    call evaluate_tetra_nedelec_first_kind( &
                        basis, point, basis_values, basis_curls, local_status)
                    if (local_status /= 0) cycle
                    field_value = matmul(basis_values, dofs)
                    arrow_length = 0.16_dp*norm2(field_value)/field_scale
                    if (arrow_length <= 1.0e-12_dp) cycle
                    arrow_end = point + arrow_length*field_value/ &
                        max(norm2(field_value), epsilon(1.0_dp))
                    arrow_side = arrow_end - point
                    call add_3d_plot( &
                        [point(1), arrow_end(1)], [point(2), arrow_end(2)], &
                        [point(3), arrow_end(3)], color="black", linewidth=1.2_dp)
                    ! A short V-shaped head makes direction readable in the
                    ! projected 3-D renderer without adding a second API.
                    call add_3d_plot( &
                        [arrow_end(1), arrow_end(1) - 0.25_dp*arrow_side(1) + &
                        0.08_dp, arrow_end(1)], &
                        [arrow_end(2), arrow_end(2) - 0.25_dp*arrow_side(2), &
                        arrow_end(2) + 0.08_dp], &
                        [arrow_end(3), arrow_end(3) - 0.25_dp*arrow_side(3), &
                        arrow_end(3)], color="black", linewidth=1.0_dp)
                end do
            end do
        end do
        call title("Tetrahedral Nedelec solution: magnitude and direction")
        call savefig(output_directory//"/nedelec_field_3d.png")
        call render_field_slice(basis, dofs)
    end subroutine render_field

    subroutine render_tetrahedron_edges()
        real(dp), parameter :: vertices(3, 4) = reshape([ &
            0.0_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 1.0_dp], [3, 4])
        integer, parameter :: edges(2, 6) = reshape([ &
            1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4], [2, 6])
        integer :: edge

        do edge = 1, size(edges, 2)
            call add_3d_plot( &
                [vertices(1, edges(1, edge)), &
                vertices(1, edges(2, edge))], &
                [vertices(2, edges(1, edge)), &
                vertices(2, edges(2, edge))], &
                [vertices(3, edges(1, edge)), &
                vertices(3, edges(2, edge))], &
                color="dimgray", linewidth=2.0_dp)
        end do
    end subroutine render_tetrahedron_edges

    subroutine render_field_slice(basis, dofs)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: dofs(:)

        integer, parameter :: slice_side = 24, vector_stride = 3
        real(dp), allocatable :: x(:), y(:), magnitude(:)
        real(dp), allocatable :: arrow_x(:), arrow_y(:), arrow_u(:), arrow_v(:)
        real(dp), allocatable :: basis_values(:, :), basis_curls(:, :)
        real(dp) :: point(3), field_value(3), field_norm
        integer :: count, arrow_count, i, j, local_status

        allocate( &
            x(slice_side**2), y(slice_side**2), magnitude(slice_side**2), &
            arrow_x((slice_side/vector_stride + 1)**2), &
            arrow_y((slice_side/vector_stride + 1)**2), &
            arrow_u((slice_side/vector_stride + 1)**2), &
            arrow_v((slice_side/vector_stride + 1)**2), &
            basis_values(3, size(dofs)), basis_curls(3, size(dofs)))
        count = 0
        arrow_count = 0
        do j = 1, slice_side
            do i = 1, slice_side
                point = [ &
                    (real(i, dp) - 0.5_dp)/real(slice_side, dp), &
                    (real(j, dp) - 0.5_dp)/real(slice_side, dp), 0.25_dp]
                if (point(1) + point(2) + point(3) > 1.0_dp) cycle
                call evaluate_tetra_nedelec_first_kind( &
                    basis, point, basis_values, basis_curls, local_status)
                if (local_status /= 0) cycle
                field_value = matmul(basis_values, dofs)
                count = count + 1
                x(count) = point(1)
                y(count) = point(2)
                magnitude(count) = norm2(field_value)
                if (mod(i - 1, vector_stride) == 0 .and. &
                    mod(j - 1, vector_stride) == 0) then
                    arrow_count = arrow_count + 1
                    arrow_x(arrow_count) = point(1)
                    arrow_y(arrow_count) = point(2)
                    arrow_u(arrow_count) = field_value(1)
                    arrow_v(arrow_count) = field_value(2)
                end if
            end do
        end do
        do i = 1, arrow_count
            field_norm = sqrt(arrow_u(i)**2 + arrow_v(i)**2)
            if (field_norm > epsilon(1.0_dp)) then
                arrow_u(i) = 0.07_dp*arrow_u(i)/field_norm
                arrow_v(i) = 0.07_dp*arrow_v(i)/field_norm
            end if
        end do
        call figure(figsize=[8.0_dp, 6.5_dp])
        call add_scatter(x(:count), y(:count), c=magnitude(:count), &
            cmap="viridis", marker=".", markersize=9.0_dp, &
            label="computed field magnitude at z=0.25")
        call plot([0.0_dp, 0.75_dp, 0.0_dp, 0.0_dp], &
            [0.0_dp, 0.0_dp, 0.75_dp, 0.0_dp], &
            color=[0.25_dp, 0.25_dp, 0.25_dp], &
            linewidth=2.0_dp, label="tetrahedron slice")
        call quiver( &
            arrow_x(:arrow_count), arrow_y(:arrow_count), &
            arrow_u(:arrow_count), arrow_v(:arrow_count), scale=1.0_dp, &
            scale_units="xy", angles="xy", color="black", width=0.003_dp, &
            headwidth=3.0_dp)
        call colorbar(label="|E_h|")
        call xlabel("x")
        call ylabel("y")
        call title("Tetrahedral Nedelec solution on the z=0.25 slice")
        call savefig(output_directory//"/nedelec_field_slice_2d.png")
    end subroutine render_field_slice

    subroutine interpolation_error(basis, error, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(out) :: error
        integer, intent(out) :: status

        real(dp), allocatable :: curls(:, :), dofs(:), values(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: difference(3), point(3)
        integer :: node

        call interpolate_reference_tetra_nedelec( &
            basis, cubic_gradient, dofs, status)
        if (status /= 0) return
        allocate( &
            values(3, tetra_nedelec_dof_count(basis)), &
            curls(3, tetra_nedelec_dof_count(basis)))
        call tetra_duffy_quadrature(12, x, y, z, weights, status)
        if (status /= 0) return
        error = 0.0_dp
        do node = 1, size(weights)
            point = [x(node), y(node), z(node)]
            call evaluate_tetra_nedelec_first_kind( &
                basis, point, values, curls, status)
            if (status /= 0) return
            difference = matmul(values, dofs) - cubic_gradient_value(point)
            error = error + weights(node) * dot_product(difference, difference)
        end do
        error = sqrt(error)
    end subroutine interpolation_error

    pure subroutine cubic_gradient(point, value)
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: value(3)

        value = cubic_gradient_value(point)
    end subroutine cubic_gradient

    pure function cubic_gradient_value(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp) :: value(3)

        ! A cubic scalar potential gives a curl-free field with all three
        ! components visible in the vector plot:
        ! phi = x^4 + (2/3)y^3 + z^3.
        value = [4.0_dp * point(1)**3, &
            2.0_dp * point(2)**2, 3.0_dp * point(3)**2]
    end function cubic_gradient_value

end program tetra_nedelec_p_convergence
