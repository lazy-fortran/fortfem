module fortfem_tetra_rt_arbitrary_order
    use fortfem_kinds, only: dp
    use fortfem_generated_tetra_rt_candidates_degree_0, only: &
        evaluate_rt_candidates_degree_0
    use fortfem_generated_tetra_rt_candidates_degree_1, only: &
        evaluate_rt_candidates_degree_1
    use fortfem_generated_tetra_rt_candidates_degree_2, only: &
        evaluate_rt_candidates_degree_2
    use fortfem_generated_tetra_rt_candidates_degree_3, only: &
        evaluate_rt_candidates_degree_3
    use fortfem_generated_tetra_rt_candidates_degree_4, only: &
        evaluate_rt_candidates_degree_4
    use fortfem_generated_tetra_rt_coefficients, only: &
        load_tetra_rt_coefficients
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_special_jacobi, only: tetrahedron_koornwinder, &
        tetrahedron_koornwinder_gradient, triangle_dubiner
    implicit none
    private

    type :: tetra_rt_t
        integer :: degree = -1
        integer :: dof_count = 0
        real(dp), allocatable :: coefficients(:, :)
    end type tetra_rt_t

    interface assignment(=)
        module procedure assign_tetra_rt
    end interface

    public :: assignment(=)
    public :: evaluate_tetra_rt
    public :: initialize_tetra_rt
    public :: tetra_rt_dof_count
    public :: tetra_rt_t

contains

    subroutine initialize_tetra_rt(degree, basis, status)
        integer, intent(in) :: degree
        type(tetra_rt_t), intent(out) :: basis
        integer, intent(out) :: status

        basis%degree = -1
        basis%dof_count = 0
        if (allocated(basis%coefficients)) deallocate(basis%coefficients)
        status = 1
        if (degree < 0) return
        basis%dof_count = 3*(degree + 1)*(degree + 2)*(degree + 3)/6 + &
            (degree + 1)*(degree + 2)/2
        if (degree <= 4) then
            call load_tetra_rt_coefficients( &
                degree, basis%coefficients, status)
        else
            call build_runtime_coefficients( &
                degree, basis%coefficients, status)
        end if
        if (status /= 0) return
        if (size(basis%coefficients, 1) /= basis%dof_count .or. &
            size(basis%coefficients, 2) /= basis%dof_count) then
            status = 2
            return
        end if
        basis%degree = degree
        status = 0
    end subroutine initialize_tetra_rt

    subroutine evaluate_tetra_rt( &
            basis, point, values, divergences, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), divergences(:)
        integer, intent(out) :: status

        real(dp), allocatable :: candidate_divergences(:)
        real(dp), allocatable :: candidate_values(:, :)
        real(dp) :: tolerance

        values = 0.0_dp
        divergences = 0.0_dp
        status = 1
        if (basis%degree < 0 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 3 .or. &
            size(values, 2) /= basis%dof_count) return
        if (size(divergences) /= basis%dof_count) return
        tolerance = 64.0_dp*epsilon(1.0_dp)
        if (any(point < -tolerance) .or. &
            sum(point) > 1.0_dp + tolerance) return
        allocate(candidate_values(3, basis%dof_count))
        allocate(candidate_divergences(basis%dof_count))
        call evaluate_candidates( &
            basis, point, candidate_values, candidate_divergences)
        values = matmul(candidate_values, basis%coefficients)
        divergences = matmul( &
            candidate_divergences, basis%coefficients)
        status = 0
    end subroutine evaluate_tetra_rt

    pure integer function tetra_rt_dof_count(basis) result(dof_count)
        type(tetra_rt_t), intent(in) :: basis

        dof_count = basis%dof_count
    end function tetra_rt_dof_count

    subroutine evaluate_candidates(basis, point, values, divergences)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), divergences(:)

        select case (basis%degree)
        case (0)
            call evaluate_rt_candidates_degree_0( &
                point(1), point(2), point(3), values, divergences)
        case (1)
            call evaluate_rt_candidates_degree_1( &
                point(1), point(2), point(3), values, divergences)
        case (2)
            call evaluate_rt_candidates_degree_2( &
                point(1), point(2), point(3), values, divergences)
        case (3)
            call evaluate_rt_candidates_degree_3( &
                point(1), point(2), point(3), values, divergences)
        case (4)
            call evaluate_rt_candidates_degree_4( &
                point(1), point(2), point(3), values, divergences)
        case default
            call evaluate_runtime_candidates( &
                basis%degree, point, values, divergences)
        end select
    end subroutine evaluate_candidates

    subroutine build_runtime_coefficients(degree, coefficients, status)
        integer, intent(in) :: degree
        real(dp), allocatable, intent(out) :: coefficients(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: candidates(:, :), divergences(:)
        real(dp), allocatable :: column_scale(:), row_scale(:)
        real(dp), allocatable :: matrix(:, :), right_hand_side(:, :)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: area_normal(3), point(3)
        integer :: component, diagonal, face, info, moment, node
        integer :: total_degree, x_degree, y_degree, z_degree

        status = 1
        allocate(matrix( &
            (degree + 1)*(degree + 2)*(degree + 4)/2, &
            (degree + 1)*(degree + 2)*(degree + 4)/2))
        allocate(right_hand_side(size(matrix, 1), size(matrix, 2)))
        allocate(candidates(3, size(matrix, 2)))
        allocate(divergences(size(matrix, 2)))
        matrix = 0.0_dp
        moment = 0
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do total_degree = 0, degree
                do x_degree = 0, total_degree
                    y_degree = total_degree - x_degree
                    moment = moment + 1
                    do node = 1, size(weights)
                        call reference_face( &
                            face, x(node), y(node), point, area_normal)
                        call evaluate_runtime_candidates( &
                            degree, point, candidates, divergences)
                        matrix(moment, :) = matrix(moment, :) + &
                            weights(node)*moment_value_triangle( &
                            degree, x_degree, y_degree, &
                            x(node), y(node))* &
                            matmul(area_normal, candidates)
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total_degree = 0, degree - 1
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            call evaluate_runtime_candidates( &
                                degree, point, candidates, divergences)
                            matrix(moment, :) = matrix(moment, :) + &
                                weights(node)*moment_value_tetrahedron( &
                                degree, x_degree, y_degree, z_degree, &
                                point)* &
                                candidates(component, :)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(matrix, 1)) return
        allocate(row_scale(size(matrix, 1)), column_scale(size(matrix, 2)))
        do moment = 1, size(matrix, 1)
            row_scale(moment) = 1.0_dp/maxval(abs(matrix(moment, :)))
            matrix(moment, :) = row_scale(moment)*matrix(moment, :)
        end do
        do component = 1, size(matrix, 2)
            column_scale(component) = &
                1.0_dp/maxval(abs(matrix(:, component)))
            matrix(:, component) = &
                column_scale(component)*matrix(:, component)
        end do
        right_hand_side = 0.0_dp
        do diagonal = 1, size(matrix, 1)
            right_hand_side(diagonal, diagonal) = row_scale(diagonal)
        end do
        allocate(coefficients(size(matrix, 1), size(matrix, 2)))
        call dense_solve(matrix, right_hand_side, coefficients, info)
        if (info /= 0) return
        do component = 1, size(coefficients, 1)
            coefficients(component, :) = &
                column_scale(component)*coefficients(component, :)
        end do
        status = 0
    end subroutine build_runtime_coefficients

    pure subroutine evaluate_runtime_candidates( &
            degree, point, values, divergences)
        integer, intent(in) :: degree
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), divergences(:)

        integer :: candidate, component, powers(3)
        integer :: total_degree, x_degree, y_degree, z_degree
        real(dp) :: gradient(3), value

        values = 0.0_dp
        divergences = 0.0_dp
        candidate = 0
        do component = 1, 3
            do total_degree = 0, degree
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        candidate = candidate + 1
                        if (degree <= 5) then
                            powers = [x_degree, y_degree, z_degree]
                            value = monomial(point, powers)
                            gradient(component) = monomial_derivative( &
                                point, powers, component)
                        else
                            value = tetrahedron_koornwinder( &
                                x_degree, y_degree, z_degree, &
                                point(1), point(2), point(3))
                            call tetrahedron_koornwinder_gradient( &
                                x_degree, y_degree, z_degree, &
                                point(1), point(2), point(3), gradient)
                        end if
                        values(component, candidate) = value
                        divergences(candidate) = gradient(component)
                    end do
                end do
            end do
        end do
        total_degree = degree
        do x_degree = 0, total_degree
            do y_degree = 0, total_degree - x_degree
                z_degree = total_degree - x_degree - y_degree
                candidate = candidate + 1
                if (degree <= 5) then
                    powers = [x_degree, y_degree, z_degree]
                    value = monomial(point, powers)
                    do component = 1, 3
                        gradient(component) = monomial_derivative( &
                            point, powers, component)
                    end do
                else
                    value = tetrahedron_koornwinder( &
                        x_degree, y_degree, z_degree, &
                        point(1), point(2), point(3))
                    call tetrahedron_koornwinder_gradient( &
                        x_degree, y_degree, z_degree, &
                        point(1), point(2), point(3), gradient)
                end if
                values(:, candidate) = point*value
                divergences(candidate) = &
                    3.0_dp*value + dot_product(point, gradient)
            end do
        end do
    end subroutine evaluate_runtime_candidates

    pure real(dp) function moment_value_triangle( &
            degree, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: degree, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (degree <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function moment_value_triangle

    pure real(dp) function moment_value_tetrahedron( &
            degree, first_degree, second_degree, third_degree, point) &
            result(value)
        integer, intent(in) :: degree
        integer, intent(in) :: first_degree, second_degree, third_degree
        real(dp), intent(in) :: point(3)

        if (degree <= 5) then
            value = point(1)**first_degree*point(2)**second_degree* &
                point(3)**third_degree
        else
            value = tetrahedron_koornwinder( &
                first_degree, second_degree, third_degree, &
                point(1), point(2), point(3))
        end if
    end function moment_value_tetrahedron

    pure real(dp) function monomial(point, powers) result(value)
        real(dp), intent(in) :: point(3)
        integer, intent(in) :: powers(3)

        value = point(1)**powers(1)*point(2)**powers(2)* &
            point(3)**powers(3)
    end function monomial

    pure real(dp) function monomial_derivative( &
            point, powers, direction) result(value)
        real(dp), intent(in) :: point(3)
        integer, intent(in) :: powers(3), direction

        integer :: reduced(3)

        if (powers(direction) == 0) then
            value = 0.0_dp
            return
        end if
        reduced = powers
        reduced(direction) = reduced(direction) - 1
        value = real(powers(direction), dp)*monomial(point, reduced)
    end function monomial_derivative

    pure subroutine reference_face(face, u, v, point, area_normal)
        integer, intent(in) :: face
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: point(3), area_normal(3)

        select case (face)
        case (1)
            point = [0.0_dp, u, v]
            area_normal = [-1.0_dp, 0.0_dp, 0.0_dp]
        case (2)
            point = [u, 0.0_dp, v]
            area_normal = [0.0_dp, -1.0_dp, 0.0_dp]
        case (3)
            point = [v, u, 0.0_dp]
            area_normal = [0.0_dp, 0.0_dp, -1.0_dp]
        case default
            point = [1.0_dp - u - v, u, v]
            area_normal = [1.0_dp, 1.0_dp, 1.0_dp]
        end select
    end subroutine reference_face

    subroutine assign_tetra_rt(left, right)
        type(tetra_rt_t), intent(out) :: left
        type(tetra_rt_t), intent(in) :: right

        left%degree = right%degree
        left%dof_count = right%dof_count
        if (allocated(left%coefficients)) deallocate(left%coefficients)
        if (allocated(right%coefficients)) then
            allocate(left%coefficients, source=right%coefficients)
        end if
    end subroutine assign_tetra_rt
end module fortfem_tetra_rt_arbitrary_order
