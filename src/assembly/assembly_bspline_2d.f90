module fortfem_assembly_bspline_2d
    !! Sparse scalar isogeometric assembly on one rational tensor patch.
    use fortfem_bspline_feec, only: &
        build_bspline_derivative_matrix, evaluate_bspline_basis, &
        evaluate_nurbs_surface_geometry
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: &
        csc_from_triplet, csc_matmul, csc_t, csc_transpose, &
        FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: assemble_bspline_h1_operator_csc
    public :: assemble_bspline_hcurl_operator_csc
    public :: assemble_bspline_hdiv_operator_csc
    public :: assemble_bspline_h1_hcurl_gradient_csc
    public :: assemble_bspline_hcurl_h1_adjoint_gradient_csc
    public :: build_bspline_feec_2d_operators_csc
    public :: scalar_weight_2d
    public :: tensor_weight_2d

    abstract interface
        pure subroutine scalar_weight_2d(point, value)
            import :: dp
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value
        end subroutine scalar_weight_2d

        pure subroutine tensor_weight_2d(point, value)
            import :: dp
            real(dp), intent(in) :: point(2)
            real(dp), intent(out) :: value(2, 2)
        end subroutine tensor_weight_2d
    end interface

contains

    subroutine assemble_bspline_hcurl_h1_adjoint_gradient_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: weak_gradient

        call assemble_bspline_h1_hcurl_gradient_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, weak_gradient, status)
        if (status%code /= 0) return
        call csc_transpose(weak_gradient, matrix, status)
    end subroutine assemble_bspline_hcurl_h1_adjoint_gradient_csc

    subroutine assemble_bspline_h1_hcurl_gradient_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: curl_incidence, gradient_incidence, hcurl_mass

        call build_bspline_feec_2d_operators_csc( &
            knots_x, knots_y, degree_x, degree_y, gradient_incidence, &
            curl_incidence, status)
        if (status%code /= 0) return
        call assemble_bspline_hcurl_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, hcurl_mass, status, curl_coefficient=0.0_dp, &
            mass_coefficient=1.0_dp)
        if (status%code /= 0) return
        call csc_matmul(hcurl_mass, gradient_incidence, matrix, status)
    end subroutine assemble_bspline_h1_hcurl_gradient_csc

    subroutine assemble_bspline_hdiv_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, divergence_coefficient, &
            mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: dx(:), dx_reduced(:), dy(:), dy_reduced(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), triplet_values(:)
        real(dp), allocatable :: vx(:), vx_reduced(:), vy(:), vy_reduced(:)
        real(dp) :: determinant, divergence_column, divergence_row
        real(dp) :: divergence_weight, geometry_jacobian(2, 2)
        real(dp) :: geometry_point(2), inverse(2, 2), mass_weight
        real(dp) :: physical_column(2), physical_row(2), physical_weight
        real(dp) :: reference_column(2), reference_row(2)
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, point_x, point_y, span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric Hdiv assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        divergence_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*degree_y + degree_x*(degree_y + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            qw_x(quadrature_order), qw_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, qw_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, qw_x)
                call build_hdiv_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), vy, dy, &
                        local_status)
                    if (local_status /= 0) return
                    call evaluate_bspline_basis( &
                        knots_y(2:size(knots_y) - 1), degree_y - 1, &
                        nodes_y(point_y), vy_reduced, dy_reduced, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), vx, dx, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_bspline_basis( &
                            knots_x(2:size(knots_x) - 1), degree_x - 1, &
                            nodes_x(point_x), vx_reduced, dx_reduced, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant*qw_x(point_x)*qw_y(point_y)
                        do local_column = 1, local_count
                            call hdiv_local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, vx, dx, vx_reduced, vy, dy, &
                                vy_reduced, reference_column, divergence_column)
                            physical_column = matmul( &
                                geometry_jacobian, reference_column)/determinant
                            divergence_column = divergence_column/determinant
                            do local_row = 1, local_count
                                call hdiv_local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, vx, dx, vx_reduced, vy, dy, &
                                    vy_reduced, reference_row, divergence_row)
                                physical_row = matmul( &
                                    geometry_jacobian, reference_row)/determinant
                                divergence_row = divergence_row/determinant
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(mass_weight*dot_product( &
                                    physical_row, physical_column) + &
                                    divergence_weight*divergence_row* &
                                    divergence_column)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            nx*(ny - 1) + (nx - 1)*ny, &
            nx*(ny - 1) + (nx - 1)*ny, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_hdiv_operator_csc

    pure subroutine build_hdiv_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx, ny
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof, x_component_count

        local_dof = 0
        x_component_count = nx*(ny - 1)
        do basis_y = span_y - degree_y, span_y - 1
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*nx
            end do
        end do
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x - 1
                local_dof = local_dof + 1
                local_dofs(local_dof) = &
                    x_component_count + basis_x + (basis_y - 1)*(nx - 1)
            end do
        end do
    end subroutine build_hdiv_local_dofs

    pure subroutine hdiv_local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, vx, dx, &
            vx_reduced, vy, dy, vy_reduced, value, divergence)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: vx(:), dx(:), vx_reduced(:)
        real(dp), intent(in) :: vy(:), dy(:), vy_reduced(:)
        real(dp), intent(out) :: value(2), divergence

        integer :: basis_x, basis_y, offset, x_local_count

        value = 0.0_dp
        x_local_count = (degree_x + 1)*degree_y
        if (local_dof <= x_local_count) then
            offset = local_dof - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x + 1)
            basis_y = span_y - degree_y + offset/(degree_x + 1)
            value(1) = vx(basis_x)*vy_reduced(basis_y)
            divergence = dx(basis_x)*vy_reduced(basis_y)
        else
            offset = local_dof - x_local_count - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x)
            basis_y = span_y - degree_y + offset/degree_x
            value(2) = vx_reduced(basis_x)*vy(basis_y)
            divergence = vx_reduced(basis_x)*dy(basis_y)
        end if
    end subroutine hdiv_local_basis_data

    subroutine assemble_bspline_hcurl_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, curl_coefficient, &
            mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: curl_coefficient, mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: dx(:), dx_reduced(:), dy(:), dy_reduced(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), triplet_values(:)
        real(dp), allocatable :: vx(:), vx_reduced(:), vy(:), vy_reduced(:)
        real(dp) :: curl_column, curl_row, curl_weight, determinant
        real(dp) :: geometry_jacobian(2, 2), geometry_point(2), inverse(2, 2)
        real(dp) :: mass_weight, physical_column(2), physical_row(2)
        real(dp) :: reference_column(2), reference_row(2), physical_weight
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, point_x, point_y, span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric Hcurl assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        curl_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(curl_coefficient)) curl_weight = curl_coefficient
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = degree_x*(degree_y + 1) + (degree_x + 1)*degree_y
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            qw_x(quadrature_order), qw_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, qw_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, qw_x)
                call build_hcurl_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), vy, dy, &
                        local_status)
                    if (local_status /= 0) return
                    call evaluate_bspline_basis( &
                        knots_y(2:size(knots_y) - 1), degree_y - 1, &
                        nodes_y(point_y), vy_reduced, dy_reduced, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), vx, dx, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_bspline_basis( &
                            knots_x(2:size(knots_x) - 1), degree_x - 1, &
                            nodes_x(point_x), vx_reduced, dx_reduced, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant*qw_x(point_x)*qw_y(point_y)
                        do local_column = 1, local_count
                            call hcurl_local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, vx, dx, vx_reduced, vy, dy, &
                                vy_reduced, reference_column, curl_column)
                            physical_column = &
                                matmul(transpose(inverse), reference_column)
                            curl_column = curl_column/determinant
                            do local_row = 1, local_count
                                call hcurl_local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, vx, dx, vx_reduced, vy, dy, &
                                    vy_reduced, reference_row, curl_row)
                                physical_row = &
                                    matmul(transpose(inverse), reference_row)
                                curl_row = curl_row/determinant
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(mass_weight*dot_product( &
                                    physical_row, physical_column) + &
                                    curl_weight*curl_row*curl_column)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            (nx - 1)*ny + nx*(ny - 1), &
            (nx - 1)*ny + nx*(ny - 1), rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_hcurl_operator_csc

    pure subroutine build_hcurl_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, ny, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx, ny
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof, x_component_count

        local_dof = 0
        x_component_count = (nx - 1)*ny
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x - 1
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*(nx - 1)
            end do
        end do
        do basis_y = span_y - degree_y, span_y - 1
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = &
                    x_component_count + basis_x + (basis_y - 1)*nx
            end do
        end do
    end subroutine build_hcurl_local_dofs

    pure subroutine hcurl_local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, vx, dx, &
            vx_reduced, vy, dy, vy_reduced, value, curl)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: vx(:), dx(:), vx_reduced(:)
        real(dp), intent(in) :: vy(:), dy(:), vy_reduced(:)
        real(dp), intent(out) :: value(2), curl

        integer :: basis_x, basis_y, offset, x_local_count

        value = 0.0_dp
        x_local_count = degree_x*(degree_y + 1)
        if (local_dof <= x_local_count) then
            offset = local_dof - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x)
            basis_y = span_y - degree_y + offset/degree_x
            value(1) = vx_reduced(basis_x)*vy(basis_y)
            curl = -vx_reduced(basis_x)*dy(basis_y)
        else
            offset = local_dof - x_local_count - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x + 1)
            basis_y = span_y - degree_y + offset/(degree_x + 1)
            value(2) = vx(basis_x)*vy_reduced(basis_y)
            curl = dx(basis_x)*vy_reduced(basis_y)
        end if
    end subroutine hcurl_local_basis_data

    subroutine build_bspline_feec_2d_operators_csc( &
            knots_x, knots_y, degree_x, degree_y, gradient, curl, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        type(csc_t), intent(out) :: gradient, curl
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), rows(:)
        real(dp), allocatable :: derivative_x(:, :), derivative_y(:, :)
        real(dp), allocatable :: values(:)
        integer :: column, entry, ix, iy, nx, ny, row, x_component_count
        integer :: local_status

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric FEEC incidence assembly failed")
        call build_bspline_derivative_matrix( &
            knots_x, degree_x, derivative_x, local_status)
        if (local_status /= 0) return
        call build_bspline_derivative_matrix( &
            knots_y, degree_y, derivative_y, local_status)
        if (local_status /= 0) return
        nx = size(derivative_x, 2)
        ny = size(derivative_y, 2)
        x_component_count = (nx - 1)*ny
        allocate(rows(2*(x_component_count + nx*(ny - 1))))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iy = 1, ny
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = ix, ix + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = column + (iy - 1)*nx
                    values(entry) = derivative_x(ix, column)
                end do
            end do
        end do
        do iy = 1, ny - 1
            do ix = 1, nx
                row = x_component_count + ix + (iy - 1)*nx
                do column = iy, iy + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = ix + (column - 1)*nx
                    values(entry) = derivative_y(iy, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            x_component_count + nx*(ny - 1), nx*ny, rows(:entry), &
            columns(:entry), values(:entry), gradient, status)
        if (status%code /= 0) return

        deallocate(rows, columns, values)
        allocate(rows(4*(nx - 1)*(ny - 1)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iy = 1, ny - 1
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = iy, iy + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = ix + (column - 1)*(nx - 1)
                    values(entry) = -derivative_y(iy, column)
                end do
                do column = ix, ix + 1
                    entry = entry + 1
                    rows(entry) = row
                    columns(entry) = x_component_count + &
                        column + (iy - 1)*nx
                    values(entry) = derivative_x(ix, column)
                end do
            end do
        end do
        call csc_from_triplet( &
            (nx - 1)*(ny - 1), x_component_count + nx*(ny - 1), &
            rows, columns, values, curl, status)
    end subroutine build_bspline_feec_2d_operators_csc

    subroutine assemble_bspline_h1_operator_csc( &
            knots_x, knots_y, degree_x, degree_y, control_points, weights, &
            quadrature_order, matrix, status, stiffness_coefficient, &
            mass_coefficient, stiffness_weight_function, mass_weight_function, &
            stiffness_tensor_function)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :), weights(:, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient
        procedure(scalar_weight_2d), optional :: stiffness_weight_function
        procedure(scalar_weight_2d), optional :: mass_weight_function
        procedure(tensor_weight_2d), optional :: stiffness_tensor_function

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: derivative_x(:), derivative_y(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: quadrature_weights_x(:)
        real(dp), allocatable :: quadrature_weights_y(:), triplet_values(:)
        real(dp), allocatable :: value_x(:), value_y(:)
        real(dp) :: determinant, geometry_jacobian(2, 2), geometry_point(2)
        real(dp) :: basis_column, basis_row, gradient_column(2), gradient_row(2)
        real(dp) :: inverse(2, 2)
        real(dp) :: mass_weight, physical_weight, stiffness_weight
        real(dp) :: mass_weight_at_point, stiffness_weight_at_point
        real(dp) :: stiffness_tensor_at_point(2, 2)
        integer :: entry, local_column, local_count, local_row
        integer :: local_status, max_entries, nx, ny, point_x, point_y
        integer :: span_x, span_y

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse isogeometric H1 assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1) return
        if (any(shape(weights) /= [nx, ny])) return
        if (size(control_points, 1) /= 2 .or. &
            any(shape(control_points(1, :, :)) /= [nx, ny])) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*(degree_y + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), &
            local_dofs(local_count), local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            quadrature_weights_x(quadrature_order), &
            quadrature_weights_y(quadrature_order))
        entry = 0
        do span_y = degree_y + 1, ny
            if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                nodes_y, quadrature_weights_y)
            do span_x = degree_x + 1, nx
                if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_x(span_x), knots_x(span_x + 1), &
                    nodes_x, quadrature_weights_x)
                call build_local_dofs( &
                    span_x, span_y, degree_x, degree_y, nx, local_dofs)
                local_matrix = 0.0_dp
                do point_y = 1, quadrature_order
                    call evaluate_bspline_basis( &
                        knots_y, degree_y, nodes_y(point_y), value_y, &
                        derivative_y, local_status)
                    if (local_status /= 0) return
                    do point_x = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_x, degree_x, nodes_x(point_x), value_x, &
                            derivative_x, local_status)
                        if (local_status /= 0) return
                        call evaluate_nurbs_surface_geometry( &
                            knots_x, knots_y, degree_x, degree_y, &
                            control_points, weights, nodes_x(point_x), &
                            nodes_y(point_y), geometry_point, &
                            geometry_jacobian, local_status)
                        if (local_status /= 0) return
                        call inverse_2d( &
                            geometry_jacobian, inverse, determinant, local_status)
                        if (local_status /= 0 .or. determinant <= 0.0_dp) return
                        physical_weight = determinant* &
                            quadrature_weights_x(point_x)* &
                            quadrature_weights_y(point_y)
                        stiffness_weight_at_point = stiffness_weight
                        mass_weight_at_point = mass_weight
                        stiffness_tensor_at_point = reshape( &
                            [1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
                        if (present(stiffness_weight_function)) then
                            call stiffness_weight_function( &
                                geometry_point, stiffness_weight_at_point)
                        end if
                        if (present(mass_weight_function)) then
                            call mass_weight_function( &
                                geometry_point, mass_weight_at_point)
                        end if
                        if (present(stiffness_tensor_function)) then
                            call stiffness_tensor_function( &
                                geometry_point, stiffness_tensor_at_point)
                        end if
                        do local_column = 1, local_count
                            call local_basis_data( &
                                local_column, span_x, span_y, degree_x, &
                                degree_y, value_x, derivative_x, value_y, &
                                derivative_y, basis_column, gradient_column)
                            gradient_column = &
                                matmul(transpose(inverse), gradient_column)
                            do local_row = 1, local_count
                                call local_basis_data( &
                                    local_row, span_x, span_y, degree_x, &
                                    degree_y, value_x, derivative_x, value_y, &
                                    derivative_y, basis_row, gradient_row)
                                gradient_row = &
                                    matmul(transpose(inverse), gradient_row)
                                local_matrix(local_row, local_column) = &
                                    local_matrix(local_row, local_column) + &
                                    physical_weight*(stiffness_weight_at_point* &
                                    dot_product(gradient_row, matmul( &
                                    stiffness_tensor_at_point, &
                                    gradient_column)) + &
                                    mass_weight_at_point*basis_row*basis_column)
                            end do
                        end do
                    end do
                end do
                do local_column = 1, local_count
                    do local_row = 1, local_count
                        entry = entry + 1
                        rows(entry) = local_dofs(local_row)
                        columns(entry) = local_dofs(local_column)
                        triplet_values(entry) = &
                            local_matrix(local_row, local_column)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            nx*ny, nx*ny, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_h1_operator_csc

    subroutine local_basis_data( &
            local_dof, span_x, span_y, degree_x, degree_y, value_x, &
            derivative_x, value_y, derivative_y, basis_value, gradient)
        integer, intent(in) :: local_dof, span_x, span_y, degree_x, degree_y
        real(dp), intent(in) :: value_x(:), derivative_x(:)
        real(dp), intent(in) :: value_y(:), derivative_y(:)
        real(dp), intent(out) :: basis_value
        real(dp), intent(out) :: gradient(2)

        integer :: basis_x, basis_y, offset_x, offset_y

        offset_x = modulo(local_dof - 1, degree_x + 1)
        offset_y = (local_dof - 1)/(degree_x + 1)
        basis_x = span_x - degree_x + offset_x
        basis_y = span_y - degree_y + offset_y
        basis_value = value_x(basis_x)*value_y(basis_y)
        gradient = [ &
            derivative_x(basis_x)*value_y(basis_y), &
            value_x(basis_x)*derivative_y(basis_y)]
    end subroutine local_basis_data

    pure subroutine build_local_dofs( &
            span_x, span_y, degree_x, degree_y, nx, local_dofs)
        integer, intent(in) :: span_x, span_y, degree_x, degree_y, nx
        integer, intent(out) :: local_dofs(:)

        integer :: basis_x, basis_y, local_dof

        local_dof = 0
        do basis_y = span_y - degree_y, span_y
            do basis_x = span_x - degree_x, span_x
                local_dof = local_dof + 1
                local_dofs(local_dof) = basis_x + (basis_y - 1)*nx
            end do
        end do
    end subroutine build_local_dofs

    pure integer function positive_span_count( &
            knots, degree, basis_count) result(count)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree, basis_count
        integer :: span

        count = 0
        do span = degree + 1, basis_count
            if (knots(span + 1) > knots(span)) count = count + 1
        end do
    end function positive_span_count

    pure subroutine inverse_2d(matrix, inverse, determinant, status)
        real(dp), intent(in) :: matrix(2, 2)
        real(dp), intent(out) :: inverse(2, 2), determinant
        integer, intent(out) :: status

        determinant = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        status = 1
        inverse = 0.0_dp
        if (abs(determinant) <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(matrix))**2)) return
        inverse = reshape([ &
            matrix(2, 2), -matrix(2, 1), &
            -matrix(1, 2), matrix(1, 1)], [2, 2])/determinant
        status = 0
    end subroutine inverse_2d

end module fortfem_assembly_bspline_2d
