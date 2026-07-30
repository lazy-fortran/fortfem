module fortfem_assembly_bspline_3d
    !! Direct sparse incidence assembly for the 3D spline de Rham complex.
    use fortfem_bspline_feec, only: &
        build_bspline_derivative_matrix, evaluate_bspline_basis, &
        evaluate_nurbs_volume_geometry
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortsparse, only: &
        csc_from_triplet, csc_matmul, csc_t, csc_transpose, &
        FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: build_bspline_feec_3d_operators_csc
    public :: assemble_bspline_h1_operator_3d_csc
    public :: assemble_bspline_hdiv_operator_3d_csc
    public :: assemble_bspline_l2_mass_3d_csc
    public :: assemble_bspline_hdiv_l2_divergence_3d_csc
    public :: assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc

contains

    subroutine assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: weak_divergence

        call assemble_bspline_hdiv_l2_divergence_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, weak_divergence, status)
        if (status%code /= 0) return
        call csc_transpose(weak_divergence, matrix, status)
    end subroutine assemble_bspline_l2_hdiv_adjoint_divergence_3d_csc

    subroutine assemble_bspline_hdiv_l2_divergence_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        type(csc_t) :: curl, divergence, gradient, l2_mass

        call build_bspline_feec_3d_operators_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            gradient, curl, divergence, status)
        if (status%code /= 0) return
        call assemble_bspline_l2_mass_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, l2_mass, status)
        if (status%code /= 0) return
        call csc_matmul(l2_mass, divergence, matrix, status)
    end subroutine assemble_bspline_hdiv_l2_divergence_3d_csc

    subroutine assemble_bspline_l2_mass_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, matrix, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: derivatives(:), local_matrix(:, :)
        real(dp), allocatable :: nodes_x(:), nodes_y(:), nodes_z(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), qw_z(:)
        real(dp), allocatable :: triplet_values(:), vx(:), vy(:), vz(:)
        real(dp) :: determinant, geometry_jacobian(3, 3), geometry_point(3)
        real(dp) :: inverse(3, 3), physical_weight
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, nz, point_x, point_y, point_z
        integer :: span_x, span_y, span_z

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse 3D isogeometric L2 mass assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. degree_z < 1) return
        if (quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        nz = size(knots_z) - degree_z - 1
        if (any(shape(weights) /= [nx, ny, nz])) return
        if (size(control_points, 1) /= 3 .or. &
            any(shape(control_points(1, :, :, :)) /= [nx, ny, nz])) return
        local_count = degree_x*degree_y*degree_z
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)* &
            positive_span_count(knots_z, degree_z, nz)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            nodes_z(quadrature_order), qw_x(quadrature_order), &
            qw_y(quadrature_order), qw_z(quadrature_order))
        entry = 0
        do span_z = degree_z + 1, nz
            if (knots_z(span_z + 1) <= knots_z(span_z)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_z(span_z), knots_z(span_z + 1), &
                nodes_z, qw_z)
            do span_y = degree_y + 1, ny
                if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                    nodes_y, qw_y)
                do span_x = degree_x + 1, nx
                    if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                    call gauss_legendre_ab( &
                        quadrature_order, knots_x(span_x), &
                        knots_x(span_x + 1), nodes_x, qw_x)
                    call build_l2_local_dofs_3d( &
                        span_x, span_y, span_z, degree_x, degree_y, degree_z, &
                        nx, ny, local_dofs)
                    local_matrix = 0.0_dp
                    do point_z = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_z(2:size(knots_z) - 1), degree_z - 1, &
                            nodes_z(point_z), vz, derivatives, local_status)
                        if (local_status /= 0) return
                        do point_y = 1, quadrature_order
                            call evaluate_bspline_basis( &
                                knots_y(2:size(knots_y) - 1), degree_y - 1, &
                                nodes_y(point_y), vy, derivatives, local_status)
                            if (local_status /= 0) return
                            do point_x = 1, quadrature_order
                                call evaluate_bspline_basis( &
                                    knots_x(2:size(knots_x) - 1), degree_x - 1, &
                                    nodes_x(point_x), vx, derivatives, &
                                    local_status)
                                if (local_status /= 0) return
                                call evaluate_nurbs_volume_geometry( &
                                    knots_x, knots_y, knots_z, degree_x, &
                                    degree_y, degree_z, control_points, weights, &
                                    nodes_x(point_x), nodes_y(point_y), &
                                    nodes_z(point_z), geometry_point, &
                                    geometry_jacobian, local_status)
                                if (local_status /= 0) return
                                call inverse_3d( &
                                    geometry_jacobian, inverse, determinant, &
                                    local_status)
                                if (local_status /= 0 .or. &
                                    determinant <= 0.0_dp) return
                                physical_weight = qw_x(point_x)*qw_y(point_y)* &
                                    qw_z(point_z)/determinant
                                do local_column = 1, local_count
                                    do local_row = 1, local_count
                                        local_matrix(local_row, local_column) = &
                                            local_matrix( &
                                            local_row, local_column) + &
                                            physical_weight*l2_local_value_3d( &
                                            local_row, span_x, span_y, span_z, &
                                            degree_x, degree_y, degree_z, &
                                            vx, vy, vz)* &
                                            l2_local_value_3d( &
                                            local_column, span_x, span_y, &
                                            span_z, degree_x, degree_y, &
                                            degree_z, vx, vy, vz)
                                    end do
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
        end do
        call csc_from_triplet( &
            (nx - 1)*(ny - 1)*(nz - 1), &
            (nx - 1)*(ny - 1)*(nz - 1), rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_l2_mass_3d_csc

    pure function l2_local_value_3d( &
            local_dof, span_x, span_y, span_z, degree_x, degree_y, degree_z, &
            vx, vy, vz) result(value)
        integer, intent(in) :: local_dof, span_x, span_y, span_z
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), intent(in) :: vx(:), vy(:), vz(:)
        real(dp) :: value

        integer :: basis_x, basis_y, basis_z, offset

        offset = local_dof - 1
        basis_x = span_x - degree_x + modulo(offset, degree_x)
        offset = offset/degree_x
        basis_y = span_y - degree_y + modulo(offset, degree_y)
        basis_z = span_z - degree_z + offset/degree_y
        value = vx(basis_x)*vy(basis_y)*vz(basis_z)
    end function l2_local_value_3d

    pure subroutine build_l2_local_dofs_3d( &
            span_x, span_y, span_z, degree_x, degree_y, degree_z, nx, ny, dofs)
        integer, intent(in) :: span_x, span_y, span_z
        integer, intent(in) :: degree_x, degree_y, degree_z, nx, ny
        integer, intent(out) :: dofs(:)

        integer :: basis_x, basis_y, basis_z, local_dof

        local_dof = 0
        do basis_z = span_z - degree_z, span_z - 1
            do basis_y = span_y - degree_y, span_y - 1
                do basis_x = span_x - degree_x, span_x - 1
                    local_dof = local_dof + 1
                    dofs(local_dof) = index_3d( &
                        basis_x, basis_y, basis_z, nx - 1, ny - 1)
                end do
            end do
        end do
    end subroutine build_l2_local_dofs_3d

    subroutine assemble_bspline_hdiv_operator_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, matrix, status, &
            divergence_coefficient, mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: divergence_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: dx(:), dxr(:), dy(:), dyr(:), dz(:), dzr(:)
        real(dp), allocatable :: local_matrix(:, :), nodes_x(:), nodes_y(:)
        real(dp), allocatable :: nodes_z(:), qw_x(:), qw_y(:), qw_z(:)
        real(dp), allocatable :: triplet_values(:), vx(:), vxr(:), vy(:)
        real(dp), allocatable :: vyr(:), vz(:), vzr(:)
        real(dp) :: determinant, divergence_column, divergence_row
        real(dp) :: divergence_weight, geometry_jacobian(3, 3)
        real(dp) :: geometry_point(3), inverse(3, 3), mass_weight
        real(dp) :: physical_column(3), physical_row(3), physical_weight
        real(dp) :: reference_column(3), reference_row(3)
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, nz, point_x, point_y, point_z
        integer :: span_x, span_y, span_z

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse 3D isogeometric Hdiv assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. degree_z < 1) return
        if (quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        nz = size(knots_z) - degree_z - 1
        if (any(shape(weights) /= [nx, ny, nz])) return
        if (size(control_points, 1) /= 3 .or. &
            any(shape(control_points(1, :, :, :)) /= [nx, ny, nz])) return
        divergence_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(divergence_coefficient)) then
            divergence_weight = divergence_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = &
            (degree_x + 1)*degree_y*degree_z + &
            degree_x*(degree_y + 1)*degree_z + &
            degree_x*degree_y*(degree_z + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)* &
            positive_span_count(knots_z, degree_z, nz)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            nodes_z(quadrature_order), qw_x(quadrature_order), &
            qw_y(quadrature_order), qw_z(quadrature_order))
        entry = 0
        do span_z = degree_z + 1, nz
            if (knots_z(span_z + 1) <= knots_z(span_z)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_z(span_z), knots_z(span_z + 1), &
                nodes_z, qw_z)
            do span_y = degree_y + 1, ny
                if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                    nodes_y, qw_y)
                do span_x = degree_x + 1, nx
                    if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                    call gauss_legendre_ab( &
                        quadrature_order, knots_x(span_x), &
                        knots_x(span_x + 1), nodes_x, qw_x)
                    call build_hdiv_local_dofs_3d( &
                        span_x, span_y, span_z, degree_x, degree_y, degree_z, &
                        nx, ny, nz, local_dofs)
                    local_matrix = 0.0_dp
                    do point_z = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_z, degree_z, nodes_z(point_z), vz, dz, &
                            local_status)
                        if (local_status /= 0) return
                        call evaluate_bspline_basis( &
                            knots_z(2:size(knots_z) - 1), degree_z - 1, &
                            nodes_z(point_z), vzr, dzr, local_status)
                        if (local_status /= 0) return
                        do point_y = 1, quadrature_order
                            call evaluate_bspline_basis( &
                                knots_y, degree_y, nodes_y(point_y), vy, dy, &
                                local_status)
                            if (local_status /= 0) return
                            call evaluate_bspline_basis( &
                                knots_y(2:size(knots_y) - 1), degree_y - 1, &
                                nodes_y(point_y), vyr, dyr, local_status)
                            if (local_status /= 0) return
                            do point_x = 1, quadrature_order
                                call evaluate_bspline_basis( &
                                    knots_x, degree_x, nodes_x(point_x), vx, &
                                    dx, local_status)
                                if (local_status /= 0) return
                                call evaluate_bspline_basis( &
                                    knots_x(2:size(knots_x) - 1), &
                                    degree_x - 1, nodes_x(point_x), vxr, dxr, &
                                    local_status)
                                if (local_status /= 0) return
                                call evaluate_nurbs_volume_geometry( &
                                    knots_x, knots_y, knots_z, degree_x, &
                                    degree_y, degree_z, control_points, weights, &
                                    nodes_x(point_x), nodes_y(point_y), &
                                    nodes_z(point_z), geometry_point, &
                                    geometry_jacobian, local_status)
                                if (local_status /= 0) return
                                call inverse_3d( &
                                    geometry_jacobian, inverse, determinant, &
                                    local_status)
                                if (local_status /= 0 .or. &
                                    determinant <= 0.0_dp) return
                                physical_weight = determinant*qw_x(point_x)* &
                                    qw_y(point_y)*qw_z(point_z)
                                do local_column = 1, local_count
                                    call hdiv_local_basis_data_3d( &
                                        local_column, span_x, span_y, span_z, &
                                        degree_x, degree_y, degree_z, vx, dx, &
                                        vxr, vy, dy, vyr, vz, dz, vzr, &
                                        reference_column, divergence_column)
                                    physical_column = matmul( &
                                        geometry_jacobian, reference_column)/ &
                                        determinant
                                    divergence_column = &
                                        divergence_column/determinant
                                    do local_row = 1, local_count
                                        call hdiv_local_basis_data_3d( &
                                            local_row, span_x, span_y, span_z, &
                                            degree_x, degree_y, degree_z, vx, &
                                            dx, vxr, vy, dy, vyr, vz, dz, &
                                            vzr, reference_row, divergence_row)
                                        physical_row = matmul( &
                                            geometry_jacobian, reference_row)/ &
                                            determinant
                                        divergence_row = &
                                            divergence_row/determinant
                                        local_matrix(local_row, local_column) = &
                                            local_matrix( &
                                            local_row, local_column) + &
                                            physical_weight*(mass_weight* &
                                            dot_product( &
                                            physical_row, physical_column) + &
                                            divergence_weight*divergence_row* &
                                            divergence_column)
                                    end do
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
        end do
        call csc_from_triplet( &
            nx*(ny - 1)*(nz - 1) + (nx - 1)*ny*(nz - 1) + &
            (nx - 1)*(ny - 1)*nz, &
            nx*(ny - 1)*(nz - 1) + (nx - 1)*ny*(nz - 1) + &
            (nx - 1)*(ny - 1)*nz, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_hdiv_operator_3d_csc

    pure subroutine build_hdiv_local_dofs_3d( &
            span_x, span_y, span_z, degree_x, degree_y, degree_z, nx, ny, nz, &
            dofs)
        integer, intent(in) :: span_x, span_y, span_z
        integer, intent(in) :: degree_x, degree_y, degree_z, nx, ny, nz
        integer, intent(out) :: dofs(:)

        integer :: basis_x, basis_y, basis_z, bx_count, by_count, local_dof

        local_dof = 0
        bx_count = nx*(ny - 1)*(nz - 1)
        by_count = (nx - 1)*ny*(nz - 1)
        do basis_z = span_z - degree_z, span_z - 1
            do basis_y = span_y - degree_y, span_y - 1
                do basis_x = span_x - degree_x, span_x
                    local_dof = local_dof + 1
                    dofs(local_dof) = &
                        index_3d(basis_x, basis_y, basis_z, nx, ny - 1)
                end do
            end do
        end do
        do basis_z = span_z - degree_z, span_z - 1
            do basis_y = span_y - degree_y, span_y
                do basis_x = span_x - degree_x, span_x - 1
                    local_dof = local_dof + 1
                    dofs(local_dof) = bx_count + &
                        index_3d(basis_x, basis_y, basis_z, nx - 1, ny)
                end do
            end do
        end do
        do basis_z = span_z - degree_z, span_z
            do basis_y = span_y - degree_y, span_y - 1
                do basis_x = span_x - degree_x, span_x - 1
                    local_dof = local_dof + 1
                    dofs(local_dof) = bx_count + by_count + &
                        index_3d(basis_x, basis_y, basis_z, nx - 1, ny - 1)
                end do
            end do
        end do
        associate(unused => nz)
        end associate
    end subroutine build_hdiv_local_dofs_3d

    pure subroutine hdiv_local_basis_data_3d( &
            local_dof, span_x, span_y, span_z, degree_x, degree_y, degree_z, &
            vx, dx, vxr, vy, dy, vyr, vz, dz, vzr, value, divergence)
        integer, intent(in) :: local_dof, span_x, span_y, span_z
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), intent(in) :: vx(:), dx(:), vxr(:)
        real(dp), intent(in) :: vy(:), dy(:), vyr(:)
        real(dp), intent(in) :: vz(:), dz(:), vzr(:)
        real(dp), intent(out) :: value(3), divergence

        integer :: basis_x, basis_y, basis_z, offset, x_count, y_count

        value = 0.0_dp
        x_count = (degree_x + 1)*degree_y*degree_z
        y_count = degree_x*(degree_y + 1)*degree_z
        if (local_dof <= x_count) then
            offset = local_dof - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x + 1)
            offset = offset/(degree_x + 1)
            basis_y = span_y - degree_y + modulo(offset, degree_y)
            basis_z = span_z - degree_z + offset/degree_y
            value(1) = vx(basis_x)*vyr(basis_y)*vzr(basis_z)
            divergence = dx(basis_x)*vyr(basis_y)*vzr(basis_z)
        else if (local_dof <= x_count + y_count) then
            offset = local_dof - x_count - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x)
            offset = offset/degree_x
            basis_y = span_y - degree_y + modulo(offset, degree_y + 1)
            basis_z = span_z - degree_z + offset/(degree_y + 1)
            value(2) = vxr(basis_x)*vy(basis_y)*vzr(basis_z)
            divergence = vxr(basis_x)*dy(basis_y)*vzr(basis_z)
        else
            offset = local_dof - x_count - y_count - 1
            basis_x = span_x - degree_x + modulo(offset, degree_x)
            offset = offset/degree_x
            basis_y = span_y - degree_y + modulo(offset, degree_y)
            basis_z = span_z - degree_z + offset/degree_y
            value(3) = vxr(basis_x)*vyr(basis_y)*vz(basis_z)
            divergence = vxr(basis_x)*vyr(basis_y)*dz(basis_z)
        end if
    end subroutine hdiv_local_basis_data_3d

    subroutine assemble_bspline_h1_operator_3d_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            control_points, weights, quadrature_order, matrix, status, &
            stiffness_coefficient, mass_coefficient)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z, quadrature_order
        real(dp), intent(in) :: control_points(:, :, :, :), weights(:, :, :)
        type(csc_t), intent(out) :: matrix
        type(fortsparse_status_t), intent(out) :: status
        real(dp), intent(in), optional :: stiffness_coefficient
        real(dp), intent(in), optional :: mass_coefficient

        integer, allocatable :: columns(:), local_dofs(:), rows(:)
        real(dp), allocatable :: dx(:), dy(:), dz(:), local_matrix(:, :)
        real(dp), allocatable :: nodes_x(:), nodes_y(:), nodes_z(:)
        real(dp), allocatable :: qw_x(:), qw_y(:), qw_z(:)
        real(dp), allocatable :: triplet_values(:), vx(:), vy(:), vz(:)
        real(dp) :: basis_column, basis_row, determinant
        real(dp) :: geometry_jacobian(3, 3), geometry_point(3)
        real(dp) :: gradient_column(3), gradient_row(3), inverse(3, 3)
        real(dp) :: mass_weight, physical_weight, stiffness_weight
        integer :: entry, local_column, local_count, local_row, local_status
        integer :: max_entries, nx, ny, nz, point_x, point_y, point_z
        integer :: span_x, span_y, span_z

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse 3D isogeometric H1 assembly failed")
        if (degree_x < 1 .or. degree_y < 1 .or. degree_z < 1) return
        if (quadrature_order < 1) return
        nx = size(knots_x) - degree_x - 1
        ny = size(knots_y) - degree_y - 1
        nz = size(knots_z) - degree_z - 1
        if (nx < degree_x + 1 .or. ny < degree_y + 1 .or. &
            nz < degree_z + 1) return
        if (any(shape(weights) /= [nx, ny, nz])) return
        if (size(control_points, 1) /= 3 .or. &
            any(shape(control_points(1, :, :, :)) /= [nx, ny, nz])) return
        stiffness_weight = 1.0_dp
        mass_weight = 0.0_dp
        if (present(stiffness_coefficient)) then
            stiffness_weight = stiffness_coefficient
        end if
        if (present(mass_coefficient)) mass_weight = mass_coefficient
        local_count = (degree_x + 1)*(degree_y + 1)*(degree_z + 1)
        max_entries = positive_span_count(knots_x, degree_x, nx)* &
            positive_span_count(knots_y, degree_y, ny)* &
            positive_span_count(knots_z, degree_z, nz)*local_count**2
        allocate( &
            rows(max_entries), columns(max_entries), &
            triplet_values(max_entries), local_dofs(local_count), &
            local_matrix(local_count, local_count), &
            nodes_x(quadrature_order), nodes_y(quadrature_order), &
            nodes_z(quadrature_order), qw_x(quadrature_order), &
            qw_y(quadrature_order), qw_z(quadrature_order))
        entry = 0
        do span_z = degree_z + 1, nz
            if (knots_z(span_z + 1) <= knots_z(span_z)) cycle
            call gauss_legendre_ab( &
                quadrature_order, knots_z(span_z), knots_z(span_z + 1), &
                nodes_z, qw_z)
            do span_y = degree_y + 1, ny
                if (knots_y(span_y + 1) <= knots_y(span_y)) cycle
                call gauss_legendre_ab( &
                    quadrature_order, knots_y(span_y), knots_y(span_y + 1), &
                    nodes_y, qw_y)
                do span_x = degree_x + 1, nx
                    if (knots_x(span_x + 1) <= knots_x(span_x)) cycle
                    call gauss_legendre_ab( &
                        quadrature_order, knots_x(span_x), &
                        knots_x(span_x + 1), nodes_x, qw_x)
                    call build_h1_local_dofs_3d( &
                        span_x, span_y, span_z, degree_x, degree_y, degree_z, &
                        nx, ny, local_dofs)
                    local_matrix = 0.0_dp
                    do point_z = 1, quadrature_order
                        call evaluate_bspline_basis( &
                            knots_z, degree_z, nodes_z(point_z), vz, dz, &
                            local_status)
                        if (local_status /= 0) return
                        do point_y = 1, quadrature_order
                            call evaluate_bspline_basis( &
                                knots_y, degree_y, nodes_y(point_y), vy, dy, &
                                local_status)
                            if (local_status /= 0) return
                            do point_x = 1, quadrature_order
                                call evaluate_bspline_basis( &
                                    knots_x, degree_x, nodes_x(point_x), vx, dx, &
                                    local_status)
                                if (local_status /= 0) return
                                call evaluate_nurbs_volume_geometry( &
                                    knots_x, knots_y, knots_z, degree_x, &
                                    degree_y, degree_z, control_points, weights, &
                                    nodes_x(point_x), nodes_y(point_y), &
                                    nodes_z(point_z), geometry_point, &
                                    geometry_jacobian, local_status)
                                if (local_status /= 0) return
                                call inverse_3d( &
                                    geometry_jacobian, inverse, determinant, &
                                    local_status)
                                if (local_status /= 0 .or. &
                                    determinant <= 0.0_dp) return
                                physical_weight = determinant*qw_x(point_x)* &
                                    qw_y(point_y)*qw_z(point_z)
                                do local_column = 1, local_count
                                    call h1_local_basis_data_3d( &
                                        local_column, span_x, span_y, span_z, &
                                        degree_x, degree_y, degree_z, vx, dx, &
                                        vy, dy, vz, dz, basis_column, &
                                        gradient_column)
                                    gradient_column = matmul( &
                                        transpose(inverse), gradient_column)
                                    do local_row = 1, local_count
                                        call h1_local_basis_data_3d( &
                                            local_row, span_x, span_y, span_z, &
                                            degree_x, degree_y, degree_z, vx, &
                                            dx, vy, dy, vz, dz, basis_row, &
                                            gradient_row)
                                        gradient_row = matmul( &
                                            transpose(inverse), gradient_row)
                                        local_matrix(local_row, local_column) = &
                                            local_matrix( &
                                            local_row, local_column) + &
                                            physical_weight*( &
                                            stiffness_weight*dot_product( &
                                            gradient_row, gradient_column) + &
                                            mass_weight*basis_row*basis_column)
                                    end do
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
        end do
        call csc_from_triplet( &
            nx*ny*nz, nx*ny*nz, rows(:entry), columns(:entry), &
            triplet_values(:entry), matrix, status)
    end subroutine assemble_bspline_h1_operator_3d_csc

    pure subroutine h1_local_basis_data_3d( &
            local_dof, span_x, span_y, span_z, degree_x, degree_y, degree_z, &
            vx, dx, vy, dy, vz, dz, value, gradient)
        integer, intent(in) :: local_dof, span_x, span_y, span_z
        integer, intent(in) :: degree_x, degree_y, degree_z
        real(dp), intent(in) :: vx(:), dx(:), vy(:), dy(:), vz(:), dz(:)
        real(dp), intent(out) :: value, gradient(3)

        integer :: basis_x, basis_y, basis_z, offset

        offset = local_dof - 1
        basis_x = span_x - degree_x + modulo(offset, degree_x + 1)
        offset = offset/(degree_x + 1)
        basis_y = span_y - degree_y + modulo(offset, degree_y + 1)
        basis_z = span_z - degree_z + offset/(degree_y + 1)
        value = vx(basis_x)*vy(basis_y)*vz(basis_z)
        gradient = [ &
            dx(basis_x)*vy(basis_y)*vz(basis_z), &
            vx(basis_x)*dy(basis_y)*vz(basis_z), &
            vx(basis_x)*vy(basis_y)*dz(basis_z)]
    end subroutine h1_local_basis_data_3d

    pure subroutine build_h1_local_dofs_3d( &
            span_x, span_y, span_z, degree_x, degree_y, degree_z, nx, ny, dofs)
        integer, intent(in) :: span_x, span_y, span_z
        integer, intent(in) :: degree_x, degree_y, degree_z, nx, ny
        integer, intent(out) :: dofs(:)

        integer :: basis_x, basis_y, basis_z, local_dof

        local_dof = 0
        do basis_z = span_z - degree_z, span_z
            do basis_y = span_y - degree_y, span_y
                do basis_x = span_x - degree_x, span_x
                    local_dof = local_dof + 1
                    dofs(local_dof) = index_3d( &
                        basis_x, basis_y, basis_z, nx, ny)
                end do
            end do
        end do
    end subroutine build_h1_local_dofs_3d

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

    pure subroutine inverse_3d(matrix, inverse, determinant, status)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp), intent(out) :: inverse(3, 3), determinant
        integer, intent(out) :: status

        determinant = &
            matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2)) - &
            matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1))
        inverse = 0.0_dp
        status = 1
        if (abs(determinant) <= 128.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(matrix))**3)) return
        inverse(1, 1) = matrix(2, 2)*matrix(3, 3) - &
            matrix(2, 3)*matrix(3, 2)
        inverse(1, 2) = matrix(1, 3)*matrix(3, 2) - &
            matrix(1, 2)*matrix(3, 3)
        inverse(1, 3) = matrix(1, 2)*matrix(2, 3) - &
            matrix(1, 3)*matrix(2, 2)
        inverse(2, 1) = matrix(2, 3)*matrix(3, 1) - &
            matrix(2, 1)*matrix(3, 3)
        inverse(2, 2) = matrix(1, 1)*matrix(3, 3) - &
            matrix(1, 3)*matrix(3, 1)
        inverse(2, 3) = matrix(1, 3)*matrix(2, 1) - &
            matrix(1, 1)*matrix(2, 3)
        inverse(3, 1) = matrix(2, 1)*matrix(3, 2) - &
            matrix(2, 2)*matrix(3, 1)
        inverse(3, 2) = matrix(1, 2)*matrix(3, 1) - &
            matrix(1, 1)*matrix(3, 2)
        inverse(3, 3) = matrix(1, 1)*matrix(2, 2) - &
            matrix(1, 2)*matrix(2, 1)
        inverse = inverse/determinant
        status = 0
    end subroutine inverse_3d

    subroutine build_bspline_feec_3d_operators_csc( &
            knots_x, knots_y, knots_z, degree_x, degree_y, degree_z, &
            gradient, curl, divergence, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:), knots_z(:)
        integer, intent(in) :: degree_x, degree_y, degree_z
        type(csc_t), intent(out) :: gradient, curl, divergence
        type(fortsparse_status_t), intent(out) :: status

        real(dp), allocatable :: dx(:, :), dy(:, :), dz(:, :), values(:)
        integer, allocatable :: columns(:), rows(:)
        integer :: bx_count, by_count, column, entry, ex_count, ey_count
        integer :: ix, iy, iz, local_status, nx, ny, nz, row

        call status_set( &
            status, FORTSPARSE_INVALID_MATRIX, &
            "Sparse 3D isogeometric FEEC incidence assembly failed")
        call build_bspline_derivative_matrix( &
            knots_x, degree_x, dx, local_status)
        if (local_status /= 0) return
        call build_bspline_derivative_matrix( &
            knots_y, degree_y, dy, local_status)
        if (local_status /= 0) return
        call build_bspline_derivative_matrix( &
            knots_z, degree_z, dz, local_status)
        if (local_status /= 0) return
        nx = size(dx, 2)
        ny = size(dy, 2)
        nz = size(dz, 2)
        ex_count = (nx - 1)*ny*nz
        ey_count = nx*(ny - 1)*nz
        bx_count = nx*(ny - 1)*(nz - 1)
        by_count = (nx - 1)*ny*(nz - 1)

        allocate(rows(2*(ex_count + ey_count + nx*ny*(nz - 1))))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iz = 1, nz
            do iy = 1, ny
                do ix = 1, nx - 1
                    row = index_3d(ix, iy, iz, nx - 1, ny)
                    do column = ix, ix + 1
                        call append_entry( &
                            row, index_3d(column, iy, iz, nx, ny), &
                            dx(ix, column), rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        do iz = 1, nz
            do iy = 1, ny - 1
                do ix = 1, nx
                    row = ex_count + index_3d(ix, iy, iz, nx, ny - 1)
                    do column = iy, iy + 1
                        call append_entry( &
                            row, index_3d(ix, column, iz, nx, ny), &
                            dy(iy, column), rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        do iz = 1, nz - 1
            do iy = 1, ny
                do ix = 1, nx
                    row = ex_count + ey_count + &
                        index_3d(ix, iy, iz, nx, ny)
                    do column = iz, iz + 1
                        call append_entry( &
                            row, index_3d(ix, iy, column, nx, ny), &
                            dz(iz, column), rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            ex_count + ey_count + nx*ny*(nz - 1), nx*ny*nz, rows, columns, &
            values, gradient, status)
        if (status%code /= 0) return

        deallocate(rows, columns, values)
        allocate(rows(4*(bx_count + by_count + (nx - 1)*(ny - 1)*nz)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx
                    row = index_3d(ix, iy, iz, nx, ny - 1)
                    do column = iy, iy + 1
                        call append_entry( &
                            row, ex_count + ey_count + &
                            index_3d(ix, column, iz, nx, ny), dy(iy, column), &
                            rows, columns, values, entry)
                    end do
                    do column = iz, iz + 1
                        call append_entry( &
                            row, ex_count + &
                            index_3d(ix, iy, column, nx, ny - 1), &
                            -dz(iz, column), rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        do iz = 1, nz - 1
            do iy = 1, ny
                do ix = 1, nx - 1
                    row = bx_count + index_3d(ix, iy, iz, nx - 1, ny)
                    do column = iz, iz + 1
                        call append_entry( &
                            row, index_3d(ix, iy, column, nx - 1, ny), &
                            dz(iz, column), rows, columns, values, entry)
                    end do
                    do column = ix, ix + 1
                        call append_entry( &
                            row, ex_count + ey_count + &
                            index_3d(column, iy, iz, nx, ny), -dx(ix, column), &
                            rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        do iz = 1, nz
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    row = bx_count + by_count + &
                        index_3d(ix, iy, iz, nx - 1, ny - 1)
                    do column = ix, ix + 1
                        call append_entry( &
                            row, ex_count + &
                            index_3d(column, iy, iz, nx, ny - 1), &
                            dx(ix, column), rows, columns, values, entry)
                    end do
                    do column = iy, iy + 1
                        call append_entry( &
                            row, index_3d(ix, column, iz, nx - 1, ny), &
                            -dy(iy, column), rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            bx_count + by_count + (nx - 1)*(ny - 1)*nz, &
            ex_count + ey_count + nx*ny*(nz - 1), rows, columns, values, &
            curl, status)
        if (status%code /= 0) return

        deallocate(rows, columns, values)
        allocate(rows(6*(nx - 1)*(ny - 1)*(nz - 1)))
        allocate(columns(size(rows)), values(size(rows)))
        entry = 0
        do iz = 1, nz - 1
            do iy = 1, ny - 1
                do ix = 1, nx - 1
                    row = index_3d(ix, iy, iz, nx - 1, ny - 1)
                    do column = ix, ix + 1
                        call append_entry( &
                            row, index_3d(column, iy, iz, nx, ny - 1), &
                            dx(ix, column), rows, columns, values, entry)
                    end do
                    do column = iy, iy + 1
                        call append_entry( &
                            row, bx_count + &
                            index_3d(ix, column, iz, nx - 1, ny), &
                            dy(iy, column), rows, columns, values, entry)
                    end do
                    do column = iz, iz + 1
                        call append_entry( &
                            row, bx_count + by_count + &
                            index_3d(ix, iy, column, nx - 1, ny - 1), &
                            dz(iz, column), rows, columns, values, entry)
                    end do
                end do
            end do
        end do
        call csc_from_triplet( &
            (nx - 1)*(ny - 1)*(nz - 1), &
            bx_count + by_count + (nx - 1)*(ny - 1)*nz, rows, columns, &
            values, divergence, status)
    end subroutine build_bspline_feec_3d_operators_csc

    pure subroutine append_entry( &
            row, column, value, rows, columns, values, entry)
        integer, intent(in) :: row, column
        real(dp), intent(in) :: value
        integer, intent(inout) :: rows(:), columns(:), entry
        real(dp), intent(inout) :: values(:)

        entry = entry + 1
        rows(entry) = row
        columns(entry) = column
        values(entry) = value
    end subroutine append_entry

    pure integer function index_3d( &
            ix, iy, iz, count_x, count_y) result(index)
        integer, intent(in) :: ix, iy, iz, count_x, count_y

        index = ix + (iy - 1)*count_x + (iz - 1)*count_x*count_y
    end function index_3d

end module fortfem_assembly_bspline_3d
