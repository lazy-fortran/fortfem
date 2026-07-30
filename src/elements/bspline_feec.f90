module fortfem_bspline_feec
    !! Tensor-product B-spline de Rham complexes on the parametric domain.
    !!
    !! The degree-reduced component spaces and commuting coefficient
    !! derivatives follow Buffa, Rivas, Sangalli, and Vázquez,
    !! SIAM J. Numer. Anal. 49 (2011), doi:10.1137/100786708.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: build_bspline_derivative_matrix
    public :: build_bspline_feec_2d_operators
    public :: evaluate_bspline_basis

contains

    subroutine evaluate_bspline_basis( &
            knots, degree, coordinate, values, derivatives, status)
        real(dp), intent(in) :: knots(:), coordinate
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: values(:), derivatives(:)
        integer, intent(out) :: status

        real(dp), allocatable :: lower(:), work(:)
        real(dp) :: denominator, evaluation_point, left, right
        integer :: basis, basis_count, level

        status = 1
        if (.not. valid_knot_vector(knots, degree)) return
        basis_count = size(knots) - degree - 1
        left = knots(degree + 1)
        right = knots(basis_count + 1)
        if (coordinate < left .or. coordinate > right) return
        evaluation_point = coordinate
        if (coordinate == right) evaluation_point = nearest(right, -1.0_dp)
        allocate(work(size(knots) - 1))
        work = 0.0_dp
        do basis = 1, size(work)
            if (knots(basis) <= evaluation_point .and. &
                evaluation_point < knots(basis + 1)) work(basis) = 1.0_dp
        end do
        do level = 1, degree
            if (level == degree) then
                allocate(lower(basis_count + 1))
                lower = work(:basis_count + 1)
            end if
            do basis = 1, size(knots) - level - 1
                left = 0.0_dp
                denominator = knots(basis + level) - knots(basis)
                if (denominator > 0.0_dp) left = &
                    (evaluation_point - knots(basis))*work(basis)/denominator
                right = 0.0_dp
                denominator = knots(basis + level + 1) - knots(basis + 1)
                if (denominator > 0.0_dp) right = &
                    (knots(basis + level + 1) - evaluation_point)* &
                    work(basis + 1)/denominator
                work(basis) = left + right
            end do
        end do
        allocate(values(basis_count), derivatives(basis_count))
        values = work(:basis_count)
        if (degree == 0) then
            derivatives = 0.0_dp
        else
            do basis = 1, basis_count
                derivatives(basis) = 0.0_dp
                denominator = knots(basis + degree) - knots(basis)
                if (denominator > 0.0_dp) derivatives(basis) = &
                    real(degree, dp)*lower(basis)/denominator
                denominator = &
                    knots(basis + degree + 1) - knots(basis + 1)
                if (denominator > 0.0_dp) derivatives(basis) = &
                    derivatives(basis) - &
                    real(degree, dp)*lower(basis + 1)/denominator
            end do
        end if
        status = 0
    end subroutine evaluate_bspline_basis

    subroutine build_bspline_derivative_matrix( &
            knots, degree, matrix, status)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp) :: scale
        integer :: basis, basis_count

        status = 1
        if (.not. valid_knot_vector(knots, degree) .or. degree < 1) return
        basis_count = size(knots) - degree - 1
        allocate(matrix(basis_count - 1, basis_count))
        matrix = 0.0_dp
        do basis = 1, basis_count - 1
            if (knots(basis + degree + 1) <= knots(basis + 1)) return
            scale = real(degree, dp)/ &
                (knots(basis + degree + 1) - knots(basis + 1))
            matrix(basis, basis) = -scale
            matrix(basis, basis + 1) = scale
        end do
        status = 0
    end subroutine build_bspline_derivative_matrix

    subroutine build_bspline_feec_2d_operators( &
            knots_x, knots_y, degree_x, degree_y, gradient, curl, status)
        real(dp), intent(in) :: knots_x(:), knots_y(:)
        integer, intent(in) :: degree_x, degree_y
        real(dp), allocatable, intent(out) :: gradient(:, :), curl(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: derivative_x(:, :), derivative_y(:, :)
        integer :: column, ix, iy, nx, ny, row, x_component_count

        status = 1
        call build_bspline_derivative_matrix( &
            knots_x, degree_x, derivative_x, status)
        if (status /= 0) return
        call build_bspline_derivative_matrix( &
            knots_y, degree_y, derivative_y, status)
        if (status /= 0) return
        nx = size(derivative_x, 2)
        ny = size(derivative_y, 2)
        x_component_count = (nx - 1)*ny
        allocate( &
            gradient(x_component_count + nx*(ny - 1), nx*ny), &
            curl((nx - 1)*(ny - 1), x_component_count + nx*(ny - 1)))
        gradient = 0.0_dp
        curl = 0.0_dp
        do iy = 1, ny
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = 1, nx
                    gradient(row, column + (iy - 1)*nx) = &
                        derivative_x(ix, column)
                end do
            end do
        end do
        do iy = 1, ny - 1
            do ix = 1, nx
                row = x_component_count + ix + (iy - 1)*nx
                do column = 1, ny
                    gradient(row, ix + (column - 1)*nx) = &
                        derivative_y(iy, column)
                end do
            end do
        end do
        do iy = 1, ny - 1
            do ix = 1, nx - 1
                row = ix + (iy - 1)*(nx - 1)
                do column = 1, ny
                    curl(row, ix + (column - 1)*(nx - 1)) = &
                        -derivative_y(iy, column)
                end do
                do column = 1, nx
                    curl(row, x_component_count + &
                        column + (iy - 1)*nx) = derivative_x(ix, column)
                end do
            end do
        end do
        status = 0
    end subroutine build_bspline_feec_2d_operators

    pure logical function valid_knot_vector(knots, degree) result(valid)
        real(dp), intent(in) :: knots(:)
        integer, intent(in) :: degree

        valid = degree >= 0 .and. size(knots) >= 2*degree + 2
        if (.not. valid) return
        valid = all(knots(2:) >= knots(:size(knots) - 1))
        if (.not. valid) return
        valid = knots(degree + 1) < knots(size(knots) - degree)
    end function valid_knot_vector

end module fortfem_bspline_feec
