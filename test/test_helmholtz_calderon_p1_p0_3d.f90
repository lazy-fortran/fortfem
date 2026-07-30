program test_helmholtz_calderon_p1_p0_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_helmholtz_calderon_p1_p0_3d, generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: wave_number = 0.7_dp
    complex(dp), allocatable :: adjoint(:, :), double_layer(:, :)
    complex(dp), allocatable :: hypersingular(:, :), single_layer(:, :)
    complex(dp), allocatable :: dirichlet(:), neumann(:), residual(:)
    real(dp), allocatable :: mass(:, :), vertices(:, :)
    integer, allocatable :: triangles(:, :)
    real(dp) :: errors(0:2), symmetry_error
    integer :: level, status
    logical :: all_passed

    all_passed = .true.
    do level = 0, 2
        call generate_sphere_surface_mesh(1.0_dp, level, vertices, triangles)
        call assemble_helmholtz_calderon_p1_p0_3d( &
            vertices, triangles, wave_number, 8, single_layer, double_layer, &
            adjoint, hypersingular, status)
        if (status /= 0) error stop "3D Helmholtz Calderon assembly failed"
        call assemble_trace_mass(vertices, triangles, mass)
        allocate(dirichlet(size(vertices, 2)))
        allocate(neumann(size(triangles, 2)))
        dirichlet = exp(cmplx(0.0_dp, wave_number, dp))
        neumann = cmplx(-1.0_dp, wave_number, dp)*dirichlet(1)
        residual = matmul(0.5_dp*mass - double_layer, dirichlet) + &
            matmul(single_layer, neumann)
        errors(level) = norm2(abs(residual))/norm2(abs(matmul(mass, dirichlet)))
        symmetry_error = max( &
            maxval(abs(single_layer - transpose(single_layer))), &
            maxval(abs(adjoint - transpose(double_layer))), &
            maxval(abs(hypersingular - transpose(hypersingular))))
        call record_condition(symmetry_error < 2.0e-12_dp, &
            "3D Helmholtz Calderon blocks have complex-bilinear symmetry")
        deallocate(dirichlet, neumann, vertices, triangles)
        deallocate(single_layer, double_layer, adjoint, hypersingular, mass)
    end do

    call record_condition( &
        errors(1) < 0.8_dp*errors(0) .and. &
        errors(2) < 0.8_dp*errors(1), &
        "outgoing spherical Calderon residual decreases under refinement")
    call record_condition(errors(2) < 8.0e-2_dp, &
        "refined Calderon projector matches the analytical outgoing monopole")
    call check_summary("Three-dimensional Helmholtz P1/P0 Calderon operators")
    if (.not. all_passed) error stop 1

contains

    subroutine assemble_trace_mass(mesh_vertices, cells, trace_mass)
        real(dp), intent(in) :: mesh_vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), allocatable, intent(out) :: trace_mass(:, :)

        real(dp) :: area
        integer :: node, triangle

        allocate(trace_mass(size(cells, 2), size(mesh_vertices, 2)))
        trace_mass = 0.0_dp
        do triangle = 1, size(cells, 2)
            area = 0.5_dp*norm2(cross_product( &
                mesh_vertices(:, cells(2, triangle)) - &
                mesh_vertices(:, cells(1, triangle)), &
                mesh_vertices(:, cells(3, triangle)) - &
                mesh_vertices(:, cells(1, triangle))))
            do node = 1, 3
                trace_mass(triangle, cells(node, triangle)) = area/3.0_dp
            end do
        end do
    end subroutine assemble_trace_mass

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_calderon_p1_p0_3d
