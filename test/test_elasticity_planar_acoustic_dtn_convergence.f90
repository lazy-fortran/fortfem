program test_elasticity_planar_acoustic_dtn_convergence
    use check, only: check_condition, check_summary
    use fortfem_api, only: solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: refinement_count = 3
    integer, parameter :: side_points(refinement_count) = [5, 9, 17]
    real(dp) :: error(refinement_count)
    integer :: refinement
    logical :: all_passed

    all_passed = .true.
    do refinement = 1, refinement_count
        call longitudinal_wave_error(side_points(refinement), &
            error(refinement))
    end do
    call record_condition(error(2) < 0.35_dp*error(1) .and. &
        error(3) < 0.35_dp*error(2), &
        "Elastic FEM-acoustic-DtN error converges quadratically")
    call record_condition(error(3) < 2.0e-3_dp, &
        "Refined elastic FEM-acoustic-DtN solution reaches target accuracy")

    call check_summary("Elasticity with planar acoustic DtN convergence")
    if (.not. all_passed) error stop 1

contains

    subroutine longitudinal_wave_error(n, rms_error)
        integer, intent(in) :: n
        real(dp), intent(out) :: rms_error

        real(dp), parameter :: angular_frequency = 2.3_dp
        real(dp), parameter :: young_modulus = 10.0_dp
        real(dp), parameter :: poisson_ratio = 0.25_dp
        real(dp), parameter :: solid_density = 1.0_dp
        real(dp), parameter :: fluid_density = 1.0_dp
        real(dp), parameter :: sound_speed = 1.0_dp
        complex(dp), allocatable :: dirichlet_values(:)
        complex(dp), allocatable :: incident_pressure(:), load(:), solution(:)
        real(dp), allocatable :: interface_normals(:, :), vertices(:, :)
        integer, allocatable :: dirichlet_dofs(:), interface_nodes(:)
        integer, allocatable :: triangles(:, :)
        complex(dp) :: incident
        real(dp) :: compressional_modulus, lambda, mu, solid_wavenumber
        real(dp) :: spacing
        integer :: dof, node, status

        mu = young_modulus/(2.0_dp*(1.0_dp + poisson_ratio))
        lambda = young_modulus*poisson_ratio/( &
            (1.0_dp + poisson_ratio)*(1.0_dp - 2.0_dp*poisson_ratio))
        compressional_modulus = lambda + 2.0_dp*mu
        solid_wavenumber = angular_frequency*sqrt( &
            solid_density/compressional_modulus)
        allocate(vertices(2, n*n), triangles(3, 2*(n - 1)**2))
        call build_square_mesh(n, vertices, triangles)
        allocate(interface_nodes(n), interface_normals(2, n))
        do node = 1, n
            interface_nodes(node) = node*n
        end do
        interface_normals = 0.0_dp
        interface_normals(1, :) = 1.0_dp
        allocate(incident_pressure(n))
        incident = cmplx( &
            -compressional_modulus*solid_wavenumber* &
            cos(solid_wavenumber), &
            fluid_density*angular_frequency*sound_speed* &
            sin(solid_wavenumber), dp)
        incident_pressure = incident

        call build_dirichlet_dofs(n, dirichlet_dofs)
        allocate(dirichlet_values(size(dirichlet_dofs)))
        do dof = 1, size(dirichlet_dofs)
            node = (dirichlet_dofs(dof) + 1)/2
            if (mod(dirichlet_dofs(dof), 2) == 1) then
                dirichlet_values(dof) = cmplx( &
                    sin(solid_wavenumber*vertices(1, node)), 0.0_dp, dp)
            else
                dirichlet_values(dof) = cmplx(0.0_dp, 0.0_dp, dp)
            end if
        end do
        allocate(load(2*n*n), solution(2*n*n))
        load = cmplx(0.0_dp, 0.0_dp, dp)
        spacing = 1.0_dp/real(n - 1, dp)
        call solve_elasticity_planar_acoustic_dtn_p1( &
            vertices, triangles, interface_nodes, interface_normals, &
            angular_frequency, sound_speed, fluid_density, &
            spacing*real(n, dp), 2, young_modulus, poisson_ratio, &
            solid_density, load, incident_pressure, dirichlet_dofs, &
            dirichlet_values, solution, status)
        call record_condition(status == 0, &
            "Elastic FEM-acoustic-DtN solve succeeds")

        rms_error = 0.0_dp
        do node = 1, n*n
            rms_error = rms_error + abs(solution(2*node - 1) - &
                sin(solid_wavenumber*vertices(1, node)))**2 + &
                abs(solution(2*node))**2
        end do
        rms_error = sqrt(rms_error/real(2*n*n, dp))
    end subroutine longitudinal_wave_error

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

    subroutine build_dirichlet_dofs(n, dofs)
        integer, intent(in) :: n
        integer, allocatable, intent(out) :: dofs(:)

        logical, allocatable :: constrained(:)
        integer :: column, count_dofs, node, row

        allocate(constrained(2*n*n))
        constrained = .false.
        do row = 1, n
            node = (row - 1)*n + 1
            constrained(2*node - 1:2*node) = .true.
        end do
        do column = 1, n
            node = column
            constrained(2*node - 1:2*node) = .true.
            node = (n - 1)*n + column
            constrained(2*node - 1:2*node) = .true.
        end do
        count_dofs = count(constrained)
        allocate(dofs(count_dofs))
        count_dofs = 0
        do node = 1, 2*n*n
            if (.not. constrained(node)) cycle
            count_dofs = count_dofs + 1
            dofs(count_dofs) = node
        end do
    end subroutine build_dirichlet_dofs

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition
end program test_elasticity_planar_acoustic_dtn_convergence
