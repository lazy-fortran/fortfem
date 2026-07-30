program test_maxwell_efie_bc_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: &
        assemble_maxwell_bc_potential_operators_3d, &
        assemble_maxwell_bc_scalar_potential_3d, &
        assemble_maxwell_efie_bc_3d, &
        build_maxwell_bc_panel_divergence, &
        build_maxwell_bc_to_refined_rwg, build_maxwell_bc_transformation
    use fortfem_kinds, only: dp
    implicit none

    complex(dp), allocatable :: full_matrix(:, :), matrix(:, :)
    complex(dp), allocatable :: scalar_potential(:, :), scaled_matrix(:, :)
    complex(dp), allocatable :: vector_potential(:, :)
    integer, allocatable :: refined_triangles(:, :)
    real(dp), allocatable :: divergence(:, :), refined_vertices(:, :)
    real(dp), allocatable :: refined_transformation(:, :), transformation(:, :)
    real(dp) :: integral, vertices(3, 4)
    integer :: basis, panel, status, triangles(3, 4)
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
    vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
    triangles(:, 1) = [1, 3, 2]
    triangles(:, 2) = [1, 2, 4]
    triangles(:, 3) = [1, 4, 3]
    triangles(:, 4) = [2, 3, 4]

    call build_maxwell_bc_panel_divergence( &
        vertices, triangles, divergence, status)
    call build_maxwell_bc_transformation( &
        vertices, triangles, refined_vertices, refined_triangles, &
        transformation, status)
    call record_condition(status == 0 .and. &
        size(divergence, 1) == size(refined_triangles, 2) .and. &
        size(divergence, 2) == size(transformation, 2), &
        "BC divergence is defined on every refined panel and primal edge")
    do basis = 1, size(divergence, 2)
        integral = 0.0_dp
        do panel = 1, size(refined_triangles, 2)
            integral = integral + divergence(panel, basis)*triangle_area( &
                refined_vertices(:, refined_triangles(:, panel)))
        end do
        call record_condition(abs(integral) < 3.0e-14_dp, &
            "each BC basis has globally conservative surface divergence")
    end do

    call assemble_maxwell_bc_scalar_potential_3d( &
        vertices, triangles, 1.2_dp, 2, 1.0e-6_dp, 2, matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(matrix - transpose(matrix))) < 3.0e-14_dp, &
        "BC scalar-potential block is complex symmetric")
    call assemble_maxwell_bc_scalar_potential_3d( &
        2.0_dp*vertices, triangles, 0.6_dp, 2, 1.0e-6_dp, 2, &
        scaled_matrix, status)
    call record_condition(maxval(abs(scaled_matrix - 0.5_dp*matrix)) < &
        2.0e-11_dp, &
        "dual-normalized BC scalar block obeys inverse-length scaling")
    call build_maxwell_bc_to_refined_rwg( &
        vertices, triangles, refined_vertices, refined_triangles, &
        refined_transformation, status)
    call record_condition(status == 0 .and. &
        size(refined_transformation, 1) == 36 .and. &
        size(refined_transformation, 2) == 6, &
        "localized BC traces compress into the conforming refined RWG space")
    call assemble_maxwell_bc_potential_operators_3d( &
        vertices, triangles, 1.2_dp, 2, 1.0e-6_dp, 2, vector_potential, &
        scalar_potential, status)
    call record_condition(maxval(abs(scalar_potential - matrix)) < 3.0e-13_dp, &
        "refined-RWG and panel-divergence BC scalar blocks agree")
    call record_condition(maxval(abs( &
        vector_potential - transpose(vector_potential))) < 3.0e-13_dp, &
        "BC vector-potential block is complex symmetric")
    call assemble_maxwell_efie_bc_3d( &
        vertices, triangles, 1.2_dp, 1.7_dp, 2, 1.0e-5_dp, 1, &
        full_matrix, status)
    call record_condition(status == 0 .and. &
        maxval(abs(full_matrix - transpose(full_matrix))) < 3.0e-13_dp .and. &
        all(abs(full_matrix) < huge(1.0_dp)), &
        "BC EFIE is a finite complex-symmetric operator")

    call check_summary("Three-dimensional Maxwell BC electric block")
    if (.not. all_passed) error stop 1

contains

    pure function triangle_area(points) result(area)
        real(dp), intent(in) :: points(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            points(:, 2) - points(:, 1), points(:, 3) - points(:, 1)))
    end function triangle_area

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

end program test_maxwell_efie_bc_3d
