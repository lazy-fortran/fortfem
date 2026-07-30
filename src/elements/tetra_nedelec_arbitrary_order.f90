module fortfem_tetra_nedelec_arbitrary_order
    use fortfem_generated_tetra_modal_vector_identities, only: &
        evaluate_tetra_modal_vector_identities
    use fortfem_generated_tetra_nedelec_candidates_order_1, only: &
        evaluate_candidates_order_1
    use fortfem_generated_tetra_nedelec_candidates_order_2, only: &
        evaluate_candidates_order_2
    use fortfem_generated_tetra_nedelec_candidates_order_3, only: &
        evaluate_candidates_order_3
    use fortfem_generated_tetra_nedelec_candidates_order_4, only: &
        evaluate_candidates_order_4
    use fortfem_generated_tetra_nedelec_coefficients, only: &
        load_tetra_nedelec_coefficients
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: dense_solve
    use fortnum_quadrature, only: gauss_legendre_ab
    use fortnum_special_jacobi, only: tetrahedron_koornwinder, &
        tetrahedron_koornwinder_gradient, triangle_dubiner
    implicit none

    private

    public :: assignment(=)
    public :: evaluate_tetra_nedelec_first_kind
    public :: initialize_tetra_nedelec_first_kind
    public :: tetra_nedelec_dof_count
    public :: tetra_nedelec_first_kind_t

    type :: tetra_nedelec_first_kind_t
        private
        integer :: order = 0
        integer :: dof_count = 0
        real(dp), allocatable :: coefficients(:, :)
    end type tetra_nedelec_first_kind_t

    interface assignment(=)
        module procedure assign_tetra_nedelec_first_kind
    end interface

contains

    subroutine initialize_tetra_nedelec_first_kind(order, basis, status)
        integer, intent(in) :: order
        type(tetra_nedelec_first_kind_t), intent(out) :: basis
        integer, intent(out) :: status

        status = 1
        if (order < 1) return

        basis%order = order
        basis%dof_count = order * (order + 2) * (order + 3) / 2
        if (order <= 4) then
            call load_tetra_nedelec_coefficients( &
                order, basis%coefficients, status)
        else
            call build_runtime_coefficients( &
                order, basis%coefficients, status)
        end if
        if (status /= 0) return
        if (size(basis%coefficients, 1) /= basis%dof_count .or. &
            size(basis%coefficients, 2) /= basis%dof_count) then
            status = 2
            return
        end if
        status = 0
    end subroutine initialize_tetra_nedelec_first_kind

    subroutine evaluate_tetra_nedelec_first_kind( &
            basis, point, values, curls, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), curls(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: candidate_curls(:, :)
        real(dp), allocatable :: candidate_values(:, :)
        real(dp) :: tolerance

        values = 0.0_dp
        curls = 0.0_dp
        status = 1
        if (basis%order < 1 .or. basis%dof_count < 1) return
        if (size(values, 1) /= 3 .or. &
            size(values, 2) /= basis%dof_count) return
        if (size(curls, 1) /= 3 .or. &
            size(curls, 2) /= basis%dof_count) return
        tolerance = 64.0_dp * epsilon(1.0_dp)
        if (any(point < -tolerance)) return
        if (sum(point) > 1.0_dp + tolerance) return

        allocate( &
            candidate_values(3, basis%dof_count), &
            candidate_curls(3, basis%dof_count))
        call evaluate_nedelec_candidates( &
            basis, point, candidate_values, candidate_curls)
        values = matmul(candidate_values, basis%coefficients)
        curls = matmul(candidate_curls, basis%coefficients)
        status = 0
    end subroutine evaluate_tetra_nedelec_first_kind

    pure function tetra_nedelec_dof_count(basis) result(dof_count)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        integer :: dof_count

        dof_count = basis%dof_count
    end function tetra_nedelec_dof_count

    pure subroutine evaluate_nedelec_candidates( &
            basis, point, values, curls)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), curls(:, :)

        select case (basis%order)
        case (1)
            call evaluate_candidates_order_1( &
                point(1), point(2), point(3), values, curls)
        case (2)
            call evaluate_candidates_order_2( &
                point(1), point(2), point(3), values, curls)
        case (3)
            call evaluate_candidates_order_3( &
                point(1), point(2), point(3), values, curls)
        case (4)
            call evaluate_candidates_order_4( &
                point(1), point(2), point(3), values, curls)
        case default
            call evaluate_runtime_candidates( &
                basis%order, point, values, curls)
        end select
    end subroutine evaluate_nedelec_candidates

    subroutine build_runtime_coefficients(order, coefficients, status)
        integer, intent(in) :: order
        real(dp), allocatable, intent(out) :: coefficients(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: candidates(:, :), curls(:, :)
        real(dp), allocatable :: column_scale(:), row_scale(:)
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: matrix(:, :), right_hand_side(:, :)
        real(dp), allocatable :: tetra_weights(:), triangle_weights(:)
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: point(3), tangent(3), tangents(3, 2)
        integer :: component, diagonal, edge, exponent, face, info, moment
        integer :: node, total_degree, x_degree, y_degree, z_degree

        status = 1
        allocate(matrix( &
            order*(order + 2)*(order + 3)/2, &
            order*(order + 2)*(order + 3)/2))
        allocate(right_hand_side(size(matrix, 1), size(matrix, 2)))
        allocate(candidates(3, size(matrix, 2)), curls(3, size(matrix, 2)))
        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        matrix = 0.0_dp
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    call evaluate_runtime_candidates( &
                        order, point, candidates, curls)
                    matrix(moment, :) = matrix(moment, :) + &
                        edge_weights(node)* &
                        shifted_legendre(exponent, edge_nodes(node))* &
                        matmul(tangent, candidates)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2*order + 4, x, y, triangle_weights, status)
        if (status /= 0) return
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            call evaluate_runtime_candidates( &
                                order, point, candidates, curls)
                            matrix(moment, :) = matrix(moment, :) + &
                                triangle_weights(node)*face_moment_value( &
                                order, x_degree, y_degree, &
                                x(node), y(node))* &
                                matmul(tangents(:, component), candidates)
                        end do
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2*order + 4, x, y, z, tetra_weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            point = [x(node), y(node), z(node)]
                            call evaluate_runtime_candidates( &
                                order, point, candidates, curls)
                            matrix(moment, :) = matrix(moment, :) + &
                                tetra_weights(node)*volume_moment_value( &
                                order, x_degree, y_degree, z_degree, &
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
            column_scale(component) = 1.0_dp/maxval(abs(matrix(:, component)))
            matrix(:, component) = column_scale(component)*matrix(:, component)
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
            order, point, values, curls)
        integer, intent(in) :: order
        real(dp), intent(in) :: point(3)
        real(dp), intent(out) :: values(:, :), curls(:, :)

        integer :: candidate, component, total_degree
        integer :: x_degree, y_degree, z_degree

        values = 0.0_dp
        curls = 0.0_dp
        candidate = 0
        do component = 1, 3
            do total_degree = 0, order - 1
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        candidate = candidate + 1
                        if (order <= 5) then
                            call add_candidate_term( &
                                candidate, component, 1.0_dp, &
                                [x_degree, y_degree, z_degree], point, &
                                values, curls)
                        else
                            call add_modal_component_candidate( &
                                candidate, component, x_degree, y_degree, &
                                z_degree, point, values, curls)
                        end if
                    end do
                end do
            end do
        end do
        total_degree = order - 1
        do component = 4, 5
            do x_degree = 0, total_degree
                do y_degree = 0, total_degree - x_degree
                    z_degree = total_degree - x_degree - y_degree
                    candidate = candidate + 1
                    if (order > 5) then
                        call add_modal_cross_candidate( &
                            candidate, component, x_degree, y_degree, &
                            z_degree, point, values, curls)
                    else if (component == 4) then
                        call add_candidate_term( &
                            candidate, 1, -1.0_dp, &
                            [x_degree, y_degree + 1, z_degree], point, &
                            values, curls)
                        call add_candidate_term( &
                            candidate, 2, 1.0_dp, &
                            [x_degree + 1, y_degree, z_degree], point, &
                            values, curls)
                    else
                        call add_candidate_term( &
                            candidate, 1, -1.0_dp, &
                            [x_degree, y_degree, z_degree + 1], point, &
                            values, curls)
                        call add_candidate_term( &
                            candidate, 3, 1.0_dp, &
                            [x_degree + 1, y_degree, z_degree], point, &
                            values, curls)
                    end if
                end do
            end do
        end do
        do y_degree = 0, total_degree
            z_degree = total_degree - y_degree
            candidate = candidate + 1
            if (order > 5) then
                call add_modal_cross_candidate( &
                    candidate, 6, 0, y_degree, z_degree, &
                    point, values, curls)
            else
                call add_candidate_term( &
                    candidate, 2, -1.0_dp, [0, y_degree, z_degree + 1], &
                    point, values, curls)
                call add_candidate_term( &
                    candidate, 3, 1.0_dp, [0, y_degree + 1, z_degree], &
                    point, values, curls)
            end if
        end do
    end subroutine evaluate_runtime_candidates

    pure subroutine add_modal_component_candidate( &
            candidate, component, first_degree, second_degree, third_degree, &
            point, values, curls)
        integer, intent(in) :: candidate, component
        integer, intent(in) :: first_degree, second_degree, third_degree
        real(dp), intent(in) :: point(3)
        real(dp), intent(inout) :: values(:, :), curls(:, :)

        real(dp) :: component_curls(3, 3), cross_curls(3, 3)
        real(dp) :: cross_values(3, 3), gradient(3), value

        value = tetrahedron_koornwinder( &
            first_degree, second_degree, third_degree, &
            point(1), point(2), point(3))
        call tetrahedron_koornwinder_gradient( &
            first_degree, second_degree, third_degree, &
            point(1), point(2), point(3), gradient)
        call evaluate_tetra_modal_vector_identities( &
            point(1), point(2), point(3), value, &
            gradient(1), gradient(2), gradient(3), &
            component_curls, cross_values, cross_curls)
        values(component, candidate) = value
        curls(:, candidate) = component_curls(:, component)
    end subroutine add_modal_component_candidate

    pure subroutine add_modal_cross_candidate( &
            candidate, family, first_degree, second_degree, third_degree, &
            point, values, curls)
        integer, intent(in) :: candidate, family
        integer, intent(in) :: first_degree, second_degree, third_degree
        real(dp), intent(in) :: point(3)
        real(dp), intent(inout) :: values(:, :), curls(:, :)

        real(dp) :: component_curls(3, 3), cross_curls(3, 3)
        real(dp) :: cross_values(3, 3), gradient(3), value, x, y, z
        integer :: family_index

        x = point(1)
        y = point(2)
        z = point(3)
        value = tetrahedron_koornwinder( &
            first_degree, second_degree, third_degree, x, y, z)
        call tetrahedron_koornwinder_gradient( &
            first_degree, second_degree, third_degree, x, y, z, gradient)
        call evaluate_tetra_modal_vector_identities( &
            x, y, z, value, gradient(1), gradient(2), gradient(3), &
            component_curls, cross_values, cross_curls)
        family_index = family - 3
        values(:, candidate) = cross_values(:, family_index)
        curls(:, candidate) = cross_curls(:, family_index)
    end subroutine add_modal_cross_candidate

    pure real(dp) function face_moment_value( &
            order, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: order, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (order <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function face_moment_value

    pure real(dp) function volume_moment_value( &
            order, first_degree, second_degree, third_degree, point) &
            result(value)
        integer, intent(in) :: order
        integer, intent(in) :: first_degree, second_degree, third_degree
        real(dp), intent(in) :: point(3)

        if (order <= 5) then
            value = point(1)**first_degree*point(2)**second_degree* &
                point(3)**third_degree
        else
            value = tetrahedron_koornwinder( &
                first_degree, second_degree, third_degree, &
                point(1), point(2), point(3))
        end if
    end function volume_moment_value

    pure subroutine add_candidate_term( &
            candidate, component, coefficient, powers, point, values, curls)
        integer, intent(in) :: candidate, component, powers(3)
        real(dp), intent(in) :: coefficient, point(3)
        real(dp), intent(inout) :: values(:, :), curls(:, :)

        real(dp) :: derivative(3), value
        integer :: direction

        value = coefficient*monomial(point, powers)
        do direction = 1, 3
            derivative(direction) = coefficient* &
                monomial_derivative(point, powers, direction)
        end do
        values(component, candidate) = values(component, candidate) + value
        select case (component)
        case (1)
            curls(2, candidate) = curls(2, candidate) + derivative(3)
            curls(3, candidate) = curls(3, candidate) - derivative(2)
        case (2)
            curls(1, candidate) = curls(1, candidate) - derivative(3)
            curls(3, candidate) = curls(3, candidate) + derivative(1)
        case (3)
            curls(1, candidate) = curls(1, candidate) + derivative(2)
            curls(2, candidate) = curls(2, candidate) - derivative(1)
        end select
    end subroutine add_candidate_term

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

    pure function shifted_legendre(degree, parameter) result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: parameter
        real(dp) :: value
        real(dp) :: coordinate, current, previous, next
        integer :: polynomial_degree

        coordinate = 2.0_dp*parameter - 1.0_dp
        if (degree == 0) then
            value = 1.0_dp
            return
        end if
        previous = 1.0_dp
        current = coordinate
        do polynomial_degree = 2, degree
            next = ( &
                real(2*polynomial_degree - 1, dp)*coordinate*current - &
                real(polynomial_degree - 1, dp)*previous)/ &
                real(polynomial_degree, dp)
            previous = current
            current = next
        end do
        value = current
    end function shifted_legendre

    pure subroutine reference_edge(edge, parameter, point, tangent)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(3), tangent(3)
        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6)

        call reference_topology(vertices, edge_vertices)
        tangent = vertices(:, edge_vertices(2, edge)) - &
            vertices(:, edge_vertices(1, edge))
        point = vertices(:, edge_vertices(1, edge)) + parameter*tangent
    end subroutine reference_edge

    pure subroutine reference_face(face, u, v, point, tangents)
        integer, intent(in) :: face
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: point(3), tangents(3, 2)
        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), face_vertices(3, 4)

        call reference_topology(vertices, edge_vertices)
        face_vertices(:, 1) = [1, 2, 3]
        face_vertices(:, 2) = [1, 2, 4]
        face_vertices(:, 3) = [1, 3, 4]
        face_vertices(:, 4) = [2, 3, 4]
        tangents(:, 1) = vertices(:, face_vertices(2, face)) - &
            vertices(:, face_vertices(1, face))
        tangents(:, 2) = vertices(:, face_vertices(3, face)) - &
            vertices(:, face_vertices(1, face))
        point = vertices(:, face_vertices(1, face)) + &
            u*tangents(:, 1) + v*tangents(:, 2)
    end subroutine reference_face

    pure subroutine reference_topology(vertices, edge_vertices)
        real(dp), intent(out) :: vertices(3, 4)
        integer, intent(out) :: edge_vertices(2, 6)

        vertices(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
        vertices(:, 3) = [0.0_dp, 1.0_dp, 0.0_dp]
        vertices(:, 4) = [0.0_dp, 0.0_dp, 1.0_dp]
        edge_vertices(:, 1) = [1, 2]
        edge_vertices(:, 2) = [1, 3]
        edge_vertices(:, 3) = [1, 4]
        edge_vertices(:, 4) = [2, 3]
        edge_vertices(:, 5) = [2, 4]
        edge_vertices(:, 6) = [3, 4]
    end subroutine reference_topology

    subroutine assign_tetra_nedelec_first_kind(left, right)
        type(tetra_nedelec_first_kind_t), intent(out) :: left
        type(tetra_nedelec_first_kind_t), intent(in) :: right

        left%order = right%order
        left%dof_count = right%dof_count
        if (allocated(right%coefficients)) then
            allocate(left%coefficients( &
                size(right%coefficients, 1), size(right%coefficients, 2)))
            left%coefficients = right%coefficients
        end if
    end subroutine assign_tetra_nedelec_first_kind

end module fortfem_tetra_nedelec_arbitrary_order
