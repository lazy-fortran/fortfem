program test_acoustic_paper_fixture
    !! Steel plate (17 mm x 240 mm, 6 x 24 cells) wetted by water on x = 0 and
    !! absorbing on its other faces, 5 MHz, geometry and materials of the
    !! paper_acoustics steel-water fixture.
    !!
    !! FortFEM uses exp(-i omega t), so a lossy medium has Im(c) < 0 and
    !! Im(E) < 0. The paper code uses exp(+i omega t) and writes the same 1 %
    !! loss as 1 + 0.01i, which under FortFEM's convention would be gain; the
    !! parameters here are its complex conjugates.
    !!
    !! Oracle: away from the lateral faces the field is a laterally uniform
    !! compressional wave, so the wet-face displacement is that of the 1D P1
    !! layer with the same cells. Its nodal values are discrete plane waves
    !! u_j = A lambda**j + B lambda**(-j), lambda = exp(i theta),
    !! cos(theta) = (1 - q/3)/(1 + q/6), q = rho_s omega**2 dx**2/M,
    !! M = lambda_L + 2 mu (plane strain), and A, B follow from the two end rows
    !!   (d/2 - i omega Z_f) u_0 + o u_1 = p      (fluid: DtN mode 0)
    !!   (d/2 - i omega Z_p) u_N + o u_{N-1} = 0  (absorbing face)
    !! with d = 2 (M/dx - rho_s omega**2 dx/3), o = -(M/dx + rho_s omega**2 dx/6),
    !! Z_f = rho_f c and Z_p = Re sqrt(rho_s M). The -i omega Z terms are the
    !! outgoing (absorbing) impedances under exp(-i omega t). With 6 cells the
    !! plate is in the P1 stop band (q > 12), so this is a check of the
    !! assembly and the convention, not of the continuum solution.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_planar_acoustic_dtn_p1
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
        1480.0_dp*cmplx(1.0_dp, -0.01_dp, dp)
    complex(dp), parameter :: young_modulus = &
        210.0e9_dp*cmplx(1.0_dp, -0.01_dp, dp)

    complex(dp), allocatable :: dirichlet_values(:), incident_pressure(:)
    complex(dp), allocatable :: load(:), solution(:)
    real(dp) :: absorbing_normals(2, absorbing_count)
    real(dp) :: interface_normals(2, interface_count)
    real(dp) :: vertices(2, node_count)
    integer :: absorbing_edges(2, absorbing_count)
    integer, allocatable :: dirichlet_dofs(:)
    integer :: interface_nodes(interface_count)
    integer :: triangles(3, 2*nx*ny)
    complex(dp) :: layer_displacement
    real(dp) :: center_amplitude, profile_variation, layer_error
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
    center_amplitude = abs(solution(2*interface_nodes(center) - 1))
    layer_displacement = layer_wet_face_displacement()
    layer_error = abs(solution(2*interface_nodes(center) - 1) - &
        layer_displacement)/abs(layer_displacement)
    write (*, "(a,es10.3)") "centre displacement error against the layer ", &
        layer_error
    call record_condition(layer_error < 1.0e-4_dp, &
        "Acoustic paper centre displacement matches the plane-wave layer")

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
        ! Outward from the solid, as the incident load expects; the paper
        ! fixture's +1 only flips the sign of u, since the DtN block is n n^T.
        interface_normals = 0.0_dp
        interface_normals(1, :) = -1.0_dp

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

    complex(dp) function layer_wet_face_displacement() result(u0)
        complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
        complex(dp) :: modulus, diagonal, off_diagonal, lambda, theta
        complex(dp) :: fluid_term, solid_term, system(2, 2), a, b
        real(dp) :: dx, inertia, poisson_factor

        dx = length_x/real(nx, dp)
        inertia = solid_density*angular_frequency**2
        poisson_factor = (1.0_dp - poisson_ratio)/ &
            ((1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
        modulus = young_modulus*poisson_factor
        diagonal = 2.0_dp*(modulus/dx - inertia*dx/3.0_dp)
        off_diagonal = -(modulus/dx + inertia*dx/6.0_dp)
        theta = acos((1.0_dp - inertia*dx**2/(3.0_dp*modulus))/ &
            (1.0_dp + inertia*dx**2/(6.0_dp*modulus)))
        lambda = exp(i*theta)
        fluid_term = diagonal/2.0_dp - &
            i*angular_frequency*fluid_density*sound_speed
        solid_term = diagonal/2.0_dp - i*angular_frequency* &
            real(sqrt(solid_density*modulus), dp)
        system(1, :) = [fluid_term + off_diagonal*lambda, &
            fluid_term + off_diagonal/lambda]
        system(2, :) = [(solid_term*lambda + off_diagonal)*lambda**(nx - 1), &
            (solid_term/lambda + off_diagonal)*lambda**(1 - nx)]
        a = system(2, 2)/(system(1, 1)*system(2, 2) - system(1, 2)*system(2, 1))
        b = -system(2, 1)/(system(1, 1)*system(2, 2) - system(1, 2)*system(2, 1))
        u0 = a + b
    end function layer_wet_face_displacement

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_acoustic_paper_fixture
