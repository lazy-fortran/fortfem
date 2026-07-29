program test_scalar_helmholtz_planar_dtn_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_scalar_helmholtz_planar_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: refinement_count = 3
    integer, parameter :: side_points(refinement_count) = [5, 9, 17]
    real(dp) :: error(refinement_count)
    integer :: refinement
    logical :: all_passed

    all_passed = .true.
    do refinement = 1, refinement_count
        call plane_wave_error(side_points(refinement), error(refinement))
    end do
    call record_condition(error(2) < 0.35_dp*error(1) .and. &
        error(3) < 0.35_dp*error(2), &
        "Scalar Helmholtz FEM-DtN error converges quadratically")
    call record_condition(error(3) < 2.0e-3_dp, &
        "Refined scalar Helmholtz FEM-DtN solution reaches target accuracy")

    call check_summary("Scalar Helmholtz planar DtN convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine plane_wave_error(n, rms_error)
        integer, intent(in) :: n
        real(dp), intent(out) :: rms_error

        real(dp), parameter :: wavenumber = 2.0_dp
        complex(dp), allocatable :: dirichlet_values(:), load(:), solution(:)
        real(dp), allocatable :: vertices(:, :)
        integer, allocatable :: dirichlet_nodes(:), dtn_nodes(:)
        integer, allocatable :: triangles(:, :)
        complex(dp) :: exact
        real(dp) :: spacing
        integer :: dirichlet_count, node, status

        allocate(vertices(2, n*n), triangles(3, 2*(n - 1)**2))
        call build_square_mesh(n, vertices, triangles)
        allocate(dtn_nodes(n))
        do node = 1, n
            dtn_nodes(node) = node*n
        end do
        dirichlet_count = 3*n - 2
        allocate(dirichlet_nodes(dirichlet_count))
        call build_dirichlet_nodes(n, dirichlet_nodes)
        allocate(dirichlet_values(dirichlet_count))
        do node = 1, dirichlet_count
            dirichlet_values(node) = exp(cmplx(0.0_dp, wavenumber* &
                vertices(1, dirichlet_nodes(node)), dp))
        end do
        allocate(load(n*n), solution(n*n))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        spacing = 1.0_dp/real(n - 1, dp)
        call solve_scalar_helmholtz_planar_dtn_p1( &
            vertices, triangles, dtn_nodes, wavenumber, &
            spacing*real(n, dp), load, dirichlet_nodes, dirichlet_values, &
            solution, status)
        call record_condition(status == 0, &
            "Scalar Helmholtz FEM-DtN solve succeeds")

        rms_error = 0.0_dp
        do node = 1, n*n
            exact = exp(cmplx(0.0_dp, &
                wavenumber*vertices(1, node), dp))
            rms_error = rms_error + abs(solution(node) - exact)**2
        end do
        rms_error = sqrt(rms_error/real(n*n, dp))
    end subroutine plane_wave_error

    subroutine build_square_mesh(n, vertices, triangles)
        integer, intent(in) :: n
        real(dp), intent(out) :: vertices(:, :)
        integer, intent(out) :: triangles(:, :)

        real(dp) :: spacing
        integer :: column, lower_left, row, triangle, vertex

        spacing = 1.0_dp/real(n - 1, dp)
        vertex = 0
        do row = 0, n - 1
            do column = 0, n - 1
                vertex = vertex + 1
                vertices(:, vertex) = spacing*real([column, row], dp)
            end do
        end do
        triangle = 0
        do row = 1, n - 1
            do column = 1, n - 1
                lower_left = column + (row - 1)*n
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + 1, lower_left + n + 1]
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + n + 1, lower_left + n]
            end do
        end do
    end subroutine build_square_mesh

    subroutine build_dirichlet_nodes(n, nodes)
        integer, intent(in) :: n
        integer, intent(out) :: nodes(:)

        integer :: column, entry, row

        entry = 0
        do row = 1, n
            entry = entry + 1
            nodes(entry) = (row - 1)*n + 1
        end do
        do column = 2, n
            entry = entry + 1
            nodes(entry) = column
        end do
        do column = 2, n
            entry = entry + 1
            nodes(entry) = (n - 1)*n + column
        end do
    end subroutine build_dirichlet_nodes

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_scalar_helmholtz_planar_dtn_convergence
