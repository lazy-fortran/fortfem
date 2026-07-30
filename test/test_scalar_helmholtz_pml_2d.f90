program test_scalar_helmholtz_pml_2d
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_scalar_helmholtz_pml_p1_2d
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 2.0_dp
    real(dp), parameter :: physical_end = 1.0_dp
    real(dp), parameter :: outer_end = 2.0_dp
    real(dp), parameter :: sigma_max = 3.0_dp
    integer, parameter :: counts(2) = [9, 17]
    complex(dp), allocatable :: boundary_values(:), exact(:), load(:)
    complex(dp), allocatable :: solution(:), stretch(:, :)
    integer, allocatable :: boundary_nodes(:), triangles(:, :)
    real(dp), allocatable :: vertices(:, :)
    real(dp) :: errors(2)
    integer :: level, status
    logical :: all_passed

    all_passed = .true.
    do level = 1, 2
        call build_mesh(counts(level), vertices, triangles, boundary_nodes)
        call build_stretch(vertices, triangles, stretch)
        call exact_solution(vertices(1, :), exact)
        allocate(load(size(vertices, 2)))
        allocate(boundary_values(size(boundary_nodes)))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        boundary_values = exact(boundary_nodes)
        call solve_scalar_helmholtz_pml_p1_2d( &
            vertices, triangles, stretch, wave_number, load, boundary_nodes, &
            boundary_values, solution, status)
        if (status /= 0) error stop "2D scalar Helmholtz PML solve failed"
        errors(level) = sqrt(sum(abs(solution - exact)**2))/ &
            sqrt(sum(abs(exact)**2))
        deallocate( &
            vertices, triangles, boundary_nodes, stretch, exact, load, &
            boundary_values, solution)
    end do

    call record_condition(errors(2) < 0.35_dp*errors(1), &
        "2D scalar PML converges under uniform refinement")
    call record_condition(errors(2) < 3.0e-3_dp, &
        "2D scalar PML matches a complex-stretched plane wave")
    call check_summary("Scalar Helmholtz PML in 2D")
    if (.not. all_passed) error stop 1

contains

    subroutine build_mesh(grid_count, vertices, triangles, boundary_nodes)
        integer, intent(in) :: grid_count
        real(dp), allocatable, intent(out) :: vertices(:, :)
        integer, allocatable, intent(out) :: triangles(:, :)
        integer, allocatable, intent(out) :: boundary_nodes(:)

        logical, allocatable :: boundary(:)
        integer :: cell, column, node, row, triangle

        allocate(vertices(2, grid_count*grid_count))
        allocate(triangles(3, 2*(grid_count - 1)**2))
        allocate(boundary(grid_count*grid_count))
        boundary = .false.
        node = 0
        do row = 0, grid_count - 1
            do column = 0, grid_count - 1
                node = node + 1
                vertices(:, node) = [ &
                    outer_end*real(column, dp)/real(grid_count - 1, dp), &
                    real(row, dp)/real(grid_count - 1, dp)]
                boundary(node) = row == 0 .or. row == grid_count - 1 .or. &
                    column == 0 .or. column == grid_count - 1
            end do
        end do
        triangle = 0
        do row = 1, grid_count - 1
            do column = 1, grid_count - 1
                cell = column + (row - 1)*grid_count
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    cell, cell + 1, cell + grid_count + 1]
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    cell, cell + grid_count + 1, cell + grid_count]
            end do
        end do
        allocate(boundary_nodes(count(boundary)))
        boundary_nodes = pack([(node, node=1, size(boundary))], boundary)
    end subroutine build_mesh

    subroutine build_stretch(vertices, triangles, stretch)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        complex(dp), allocatable, intent(out) :: stretch(:, :)

        real(dp) :: midpoint, sigma
        integer :: element

        allocate(stretch(2, size(triangles, 2)))
        stretch = cmplx(1.0_dp, 0.0_dp, dp)
        do element = 1, size(triangles, 2)
            midpoint = sum(vertices(1, triangles(:, element)))/3.0_dp
            if (midpoint <= physical_end) cycle
            sigma = sigma_max*((midpoint - physical_end)/ &
                (outer_end - physical_end))**2
            stretch(1, element) = cmplx(1.0_dp, sigma/wave_number, dp)
        end do
    end subroutine build_stretch

    subroutine exact_solution(x, values)
        real(dp), intent(in) :: x(:)
        complex(dp), allocatable, intent(out) :: values(:)

        complex(dp) :: stretched_x
        integer :: node

        allocate(values(size(x)))
        do node = 1, size(x)
            stretched_x = cmplx(x(node), 0.0_dp, dp)
            if (x(node) > physical_end) then
                stretched_x = cmplx(x(node), sigma_max* &
                    (x(node) - physical_end)**3/( &
                    3.0_dp*wave_number*(outer_end - physical_end)**2), dp)
            end if
            values(node) = exp(cmplx(0.0_dp, wave_number, dp)*stretched_x)
        end do
    end subroutine exact_solution

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_scalar_helmholtz_pml_2d
