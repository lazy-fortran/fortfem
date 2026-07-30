module fortfem_assembly_bspline_3d
    !! Direct sparse incidence assembly for the 3D spline de Rham complex.
    use fortfem_bspline_feec, only: build_bspline_derivative_matrix
    use fortfem_kinds, only: dp
    use fortsparse, only: &
        csc_from_triplet, csc_t, FORTSPARSE_INVALID_MATRIX, &
        fortsparse_status_t, status_set
    implicit none
    private

    public :: build_bspline_feec_3d_operators_csc

contains

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
