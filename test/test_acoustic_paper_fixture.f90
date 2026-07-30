program test_acoustic_paper_fixture
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 6, ny = 24
    integer, parameter :: node_count = (nx + 1)*(ny + 1)
    integer, parameter :: interface_count = ny + 1
    integer, parameter :: absorbing_count = 2*nx + ny
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: frequency = 5.0e6_dp
    real(dp), parameter :: angular_frequency = 2.0_dp*pi*frequency
    real(dp), parameter :: length_x = 17.0e-3_dp
    real(dp), parameter :: length_y = 240.0e-3_dp
    real(dp), parameter :: fluid_density = 1000.0_dp
    real(dp), parameter :: solid_density = 8000.0_dp
    real(dp), parameter :: poisson_ratio = 0.29_dp
    complex(dp), parameter :: sound_speed = &
        1480.0_dp*cmplx(1.0_dp, 0.01_dp, dp)
    complex(dp), parameter :: young_modulus = &
        210.0e9_dp*cmplx(1.0_dp, 0.01_dp, dp)

    complex(dp), allocatable :: dirichlet_values(:), incident_pressure(:)
    complex(dp), allocatable :: load(:), solution(:)
    real(dp) :: absorbing_normals(2, absorbing_count)
    real(dp) :: interface_normals(2, interface_count)
    real(dp) :: vertices(2, node_count)
    integer :: absorbing_edges(2, absorbing_count)
    integer, allocatable :: dirichlet_dofs(:)
    integer :: interface_nodes(interface_count)
    integer :: triangles(3, 2*nx*ny)
    ! paper_acoustics 6300ab0 fem_dtn_solver, nx=6, ny=24,
    ! boundary_samples=24, nfmax=11:
    real(dp), parameter :: paper_center_amplitude = &
        1.590446791623552e-16_dp
    real(dp) :: center_amplitude, profile_variation
    integer :: center, node, status
    logical :: all_passed

    all_passed = .true.
    call build_mesh()
    allocate(load(2*node_count), solution(2*node_count))
    allocate(incident_pressure(interface_count))
    allocate(dirichlet_dofs(0), dirichlet_values(0))
    load = cmplx(0.0_dp, 0.0_dp, dp)
    incident_pressure = cmplx(1.0_dp, 0.0_dp, dp)
    call solve_elasticity_planar_acoustic_dtn_p1( &
        vertices, triangles, interface_nodes, interface_normals, &
        angular_frequency, sound_speed, fluid_density, &
        length_y*real(interface_count, dp)/real(ny, dp), ny/2, &
        young_modulus, poisson_ratio, solid_density, load, &
        incident_pressure, dirichlet_dofs, dirichlet_values, solution, &
        status, absorbing_edges=absorbing_edges, &
        absorbing_normals=absorbing_normals)
    call record_condition(status == 0, &
        "Acoustic paper steel-water fixture solves")

    center = (interface_count + 1)/2
    center_amplitude = 0.0_dp
    do node = center - 1, center + 2
        center_amplitude = center_amplitude + &
            abs(solution(2*interface_nodes(node) - 1))/4.0_dp
    end do
    call record_condition(abs(center_amplitude - paper_center_amplitude)/ &
        paper_center_amplitude < 2.0e-3_dp, &
        "Acoustic paper center displacement matches its Fortran reference")

    profile_variation = 0.0_dp
    do node = 3, interface_count - 2
        profile_variation = max(profile_variation, abs( &
            abs(solution(2*interface_nodes(node) - 1)) - center_amplitude))
    end do
    call record_condition(profile_variation/center_amplitude < 0.15_dp, &
        "Acoustic paper interface has the expected plane-wave plateau")

    call check_summary("Acoustic paper FEM-DtN fixture")
    if (.not. all_passed) error stop 1

contains

    subroutine build_mesh()
        real(dp) :: dx, dy
        integer :: column, edge, lower_left, row, triangle, vertex

        dx = length_x/real(nx, dp)
        dy = length_y/real(ny, dp)
        vertex = 0
        do row = 0, ny
            do column = 0, nx
                vertex = vertex + 1
                vertices(:, vertex) = [ &
                    dx*real(column, dp), dy*real(row, dp)]
            end do
        end do
        triangle = 0
        do row = 1, ny
            do column = 1, nx
                lower_left = column + (row - 1)*(nx + 1)
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + 1, lower_left + nx + 2]
                triangle = triangle + 1
                triangles(:, triangle) = [ &
                    lower_left, lower_left + nx + 2, lower_left + nx + 1]
            end do
        end do
        do row = 1, interface_count
            interface_nodes(row) = 1 + (row - 1)*(nx + 1)
        end do
        interface_normals = 0.0_dp
        interface_normals(1, :) = 1.0_dp

        edge = 0
        do column = 1, nx
            edge = edge + 1
            absorbing_edges(:, edge) = [column, column + 1]
            absorbing_normals(:, edge) = [0.0_dp, -1.0_dp]
        end do
        do row = 1, ny
            edge = edge + 1
            absorbing_edges(:, edge) = [ &
                row*(nx + 1), (row + 1)*(nx + 1)]
            absorbing_normals(:, edge) = [1.0_dp, 0.0_dp]
        end do
        do column = nx, 1, -1
            edge = edge + 1
            absorbing_edges(:, edge) = [ &
                ny*(nx + 1) + column + 1, ny*(nx + 1) + column]
            absorbing_normals(:, edge) = [0.0_dp, 1.0_dp]
        end do
    end subroutine build_mesh

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_acoustic_paper_fixture
