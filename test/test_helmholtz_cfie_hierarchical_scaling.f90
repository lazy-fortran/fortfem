program test_helmholtz_cfie_hierarchical_scaling
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: &
        apply_helmholtz_cfie_p0_hierarchical_3d, &
        assemble_helmholtz_double_layer_p0_3d, &
        assemble_helmholtz_single_layer_p0_3d
    use fortfem_core, only: &
        generate_sphere_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: dense_action(:), density(:), fast_action(:)
    complex(dp), allocatable :: double_layer(:, :), matrix(:, :)
    complex(dp), allocatable :: single_layer(:, :)
    real(dp), allocatable :: areas(:), vertices(:, :)
    integer, allocatable :: triangles(:, :)
    integer :: interactions(3), level, panel, panel_counts(3), status
    real(dp) :: relative_error, scaling_exponent
    logical :: all_passed

    all_passed = .true.
    do level = 2, 4
        call generate_sphere_surface_mesh(1.0_dp, level, vertices, triangles)
        panel_counts(level - 1) = size(triangles, 2)
        allocate(density(size(triangles, 2)))
        density = cmplx(1.0_dp, 0.25_dp, dp)
        call apply_helmholtz_cfie_p0_hierarchical_3d( &
            vertices, triangles, density, acos(-1.0_dp), acos(-1.0_dp), &
            0.45_dp, 6, fast_action, status, interactions(level - 1))
        call record_condition(status == 0, &
            "Hierarchical sphere CFIE action succeeds")
        if (level == 2) then
            call assemble_helmholtz_single_layer_p0_3d( &
                vertices, triangles, acos(-1.0_dp), 8, single_layer, status)
            call assemble_helmholtz_double_layer_p0_3d( &
                vertices, triangles, acos(-1.0_dp), 8, double_layer, status)
            areas = triangle_areas(vertices, triangles)
            matrix = double_layer - &
                cmplx(0.0_dp, acos(-1.0_dp), dp)*single_layer
            do concurrent (panel = 1:size(areas))
                matrix(panel, panel) = matrix(panel, panel) + 0.5_dp*areas(panel)
            end do
            dense_action = matmul(matrix, density)
            relative_error = norm2(abs(fast_action - dense_action))/ &
                norm2(abs(dense_action))
            print '(a,es12.4)', "CFIE action relative error ", relative_error
            call record_condition(relative_error < 1.0e-1_dp, &
                "Hierarchical CFIE action agrees with dense Galerkin assembly")
        end if
        deallocate(density)
    end do

    scaling_exponent = log(real(interactions(3), dp)/interactions(2))/ &
        log(real(panel_counts(3), dp)/panel_counts(2))
    print '(a,3(i0,1x),a,3(i0,1x),a,f6.3)', &
        "CFIE panels ", panel_counts, " interactions ", interactions, &
        "scaling exponent ", scaling_exponent
    call record_condition(scaling_exponent < 1.8_dp, &
        "Hierarchical CFIE interactions grow subquadratically")

    call check_summary("Hierarchical Helmholtz CFIE scaling")
    if (.not. all_passed) error stop 1

contains

    pure function triangle_areas(vertices, triangles) result(result)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: triangles(:, :)
        real(dp) :: result(size(triangles, 2))

        real(dp) :: first(3), second(3)
        integer :: triangle

        do triangle = 1, size(triangles, 2)
            first = vertices(:, triangles(2, triangle)) - &
                vertices(:, triangles(1, triangle))
            second = vertices(:, triangles(3, triangle)) - &
                vertices(:, triangles(1, triangle))
            result(triangle) = 0.5_dp*norm2([ &
                first(2)*second(3) - first(3)*second(2), &
                first(3)*second(1) - first(1)*second(3), &
                first(1)*second(2) - first(2)*second(1)])
        end do
    end function triangle_areas

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_helmholtz_cfie_hierarchical_scaling
