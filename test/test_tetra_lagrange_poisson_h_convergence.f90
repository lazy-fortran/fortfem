program test_tetra_lagrange_poisson_h_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_tetra_lagrange_solution_prepared, &
        initialize_tetra_lagrange_solution_evaluator, &
        solve_tetra_lagrange_poisson, tetra_duffy_quadrature, &
        tetra_lagrange_solution_evaluator_t
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortsparse, only: fortsparse_status_t
    implicit none

    integer, parameter :: refinements(3, 4) = reshape([ &
        2, 3, 4, &
        1, 2, 3, &
        1, 2, 3, &
        1, 2, 3], [3, 4])
    real(dp) :: errors(2, 3, 4), rates(2)
    logical :: all_passed
    integer :: degree, level

    all_passed = .true.
    do degree = 1, 4
        do level = 1, size(refinements, 1)
            call solve_and_measure( &
                refinements(level, degree), degree, errors(:, level, degree))
        end do
        rates = log(errors(:, 2, degree)/errors(:, 3, degree))/ &
            log(real(refinements(3, degree), dp)/ &
            real(refinements(2, degree), dp))
        write (*, '(a,i0,a,2(es12.4,1x),a,2(f8.3,1x))') &
            "H1 degree ", degree, " fine errors ", errors(:, 3, degree), &
            " rates ", rates
        call record_condition( &
            all(errors(:, 2:, degree) < errors(:, :2, degree)), &
            "tetrahedral H1 errors decrease under mesh refinement")
        call record_condition( &
            rates(1) > 0.75_dp*real(degree + 1, dp), &
            "tetrahedral H1 L2 error reaches the expected h order")
        call record_condition( &
            rates(2) > 0.75_dp*real(degree, dp), &
            "tetrahedral H1 gradient error reaches the expected h order")
    end do

    call check_summary("tetrahedral H1 Poisson h convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine solve_and_measure(refinement, degree, error)
        integer, intent(in) :: refinement, degree
        real(dp), intent(out) :: error(2)

        integer, allocatable :: tetrahedra(:, :)
        real(dp), allocatable :: solution(:), vertices(:, :)
        type(fortsparse_status_t) :: status

        call build_cube_mesh(refinement, vertices, tetrahedra)
        call solve_tetra_lagrange_poisson( &
            vertices, tetrahedra, degree, manufactured_source, &
            manufactured_boundary, solution, status)
        call record_condition( &
            status%code == 0, "refined tetrahedral H1 solve succeeds")
        if (status%code /= 0) then
            error = huge(1.0_dp)
            return
        end if
        call measure_error(vertices, tetrahedra, degree, solution, error)
    end subroutine solve_and_measure

    subroutine measure_error(vertices, tetrahedra, degree, solution, error)
        real(dp), intent(in) :: vertices(:, :), solution(:)
        integer, intent(in) :: tetrahedra(:, :), degree
        real(dp), intent(out) :: error(2)

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        type(tetra_lagrange_solution_evaluator_t) :: evaluator
        real(dp) :: exact_gradient(3), gradient(3), jacobian(3, 3)
        real(dp) :: point(3), value, weight
        integer :: cell, node, quadrature_status, sample, status

        call tetra_duffy_quadrature( &
            2*degree + 6, x, y, z, weights, quadrature_status)
        if (quadrature_status /= 0) error stop "H1 error quadrature failed"
        call initialize_tetra_lagrange_solution_evaluator( &
            vertices, tetrahedra, degree, evaluator, status)
        if (status /= 0) error stop "H1 evaluator initialization failed"
        error = 0.0_dp
        do cell = 1, size(tetrahedra, 2)
            do node = 1, 3
                jacobian(:, node) = &
                    vertices(:, tetrahedra(node + 1, cell)) - &
                    vertices(:, tetrahedra(1, cell))
            end do
            do sample = 1, size(weights)
                point = vertices(:, tetrahedra(1, cell)) + &
                    matmul(jacobian, [x(sample), y(sample), z(sample)])
                call evaluate_tetra_lagrange_solution_prepared( &
                    evaluator, solution, cell, &
                    [x(sample), y(sample), z(sample)], value, gradient, status)
                if (status /= 0) error stop "H1 reconstruction failed"
                exact_gradient = manufactured_gradient(point)
                weight = det3(jacobian)*weights(sample)
                error(1) = error(1) + &
                    weight*(value - manufactured_value(point))**2
                error(2) = error(2) + &
                    weight*sum((gradient - exact_gradient)**2)
            end do
        end do
        error = sqrt(error)
    end subroutine measure_error

    subroutine build_cube_mesh(refinement, vertices, tetrahedra)
        integer, intent(in) :: refinement
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: tetrahedra(:, :)

        integer :: cube(8), i, j, k, cell, vertex
        real(dp) :: jacobian(3, 3)

        allocate (vertices(3, (refinement + 1)**3))
        do k = 0, refinement
            do j = 0, refinement
                do i = 0, refinement
                    vertex = vertex_index(i, j, k, refinement)
                    vertices(:, vertex) = real([i, j, k], dp)/ &
                        real(refinement, dp)
                end do
            end do
        end do
        allocate (tetrahedra(4, 6*refinement**3))
        cell = 0
        do k = 0, refinement - 1
            do j = 0, refinement - 1
                do i = 0, refinement - 1
                    cube = [ &
                        vertex_index(i, j, k, refinement), &
                        vertex_index(i + 1, j, k, refinement), &
                        vertex_index(i, j + 1, k, refinement), &
                        vertex_index(i + 1, j + 1, k, refinement), &
                        vertex_index(i, j, k + 1, refinement), &
                        vertex_index(i + 1, j, k + 1, refinement), &
                        vertex_index(i, j + 1, k + 1, refinement), &
                        vertex_index(i + 1, j + 1, k + 1, refinement)]
                    call add_cube_tetrahedra(cube, tetrahedra, cell)
                end do
            end do
        end do
        do cell = 1, size(tetrahedra, 2)
            jacobian(:, 1) = vertices(:, tetrahedra(2, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 2) = vertices(:, tetrahedra(3, cell)) - &
                vertices(:, tetrahedra(1, cell))
            jacobian(:, 3) = vertices(:, tetrahedra(4, cell)) - &
                vertices(:, tetrahedra(1, cell))
            if (det3(jacobian) < 0.0_dp) then
                vertex = tetrahedra(3, cell)
                tetrahedra(3, cell) = tetrahedra(4, cell)
                tetrahedra(4, cell) = vertex
            end if
        end do
    end subroutine build_cube_mesh

    pure integer function vertex_index(i, j, k, refinement) result(index)
        integer, intent(in) :: i, j, k, refinement

        index = 1 + i + (refinement + 1)*(j + (refinement + 1)*k)
    end function vertex_index

    subroutine add_cube_tetrahedra(cube, tetrahedra, cell_count)
        integer, intent(in) :: cube(8)
        integer, intent(inout) :: tetrahedra(:, :), cell_count

        integer :: local_tetrahedra(4, 6), local_cell

        local_tetrahedra(:, 1) = cube([1, 2, 4, 8])
        local_tetrahedra(:, 2) = cube([1, 2, 6, 8])
        local_tetrahedra(:, 3) = cube([1, 3, 4, 8])
        local_tetrahedra(:, 4) = cube([1, 3, 7, 8])
        local_tetrahedra(:, 5) = cube([1, 5, 6, 8])
        local_tetrahedra(:, 6) = cube([1, 5, 7, 8])
        do local_cell = 1, 6
            cell_count = cell_count + 1
            tetrahedra(:, cell_count) = local_tetrahedra(:, local_cell)
        end do
    end subroutine add_cube_tetrahedra

    pure real(dp) function manufactured_value(point) result(value)
        real(dp), intent(in) :: point(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = product(sin(pi*point))
    end function manufactured_value

    pure function manufactured_gradient(point) result(gradient)
        real(dp), intent(in) :: point(3)
        real(dp) :: gradient(3)
        real(dp), parameter :: pi = acos(-1.0_dp)

        gradient = pi*[ &
            cos(pi*point(1))*sin(pi*point(2))*sin(pi*point(3)), &
            sin(pi*point(1))*cos(pi*point(2))*sin(pi*point(3)), &
            sin(pi*point(1))*sin(pi*point(2))*cos(pi*point(3))]
    end function manufactured_gradient

    pure subroutine manufactured_source(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value
        real(dp), parameter :: pi = acos(-1.0_dp)

        value = 3.0_dp*pi*pi*manufactured_value([x, y, z])
    end subroutine manufactured_source

    pure subroutine manufactured_boundary(x, y, z, value)
        real(dp), intent(in) :: x, y, z
        real(dp), intent(out) :: value

        value = manufactured_value([x, y, z])
    end subroutine manufactured_boundary

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_lagrange_poisson_h_convergence
