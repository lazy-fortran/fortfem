module fortfem_tetra_nedelec_arbitrary_order
    use fortfem_generated_tetra_nedelec_candidates_order_1, only: &
        evaluate_candidates_order_1
    use fortfem_generated_tetra_nedelec_candidates_order_2, only: &
        evaluate_candidates_order_2
    use fortfem_generated_tetra_nedelec_candidates_order_3, only: &
        evaluate_candidates_order_3
    use fortfem_generated_tetra_nedelec_candidates_order_4, only: &
        evaluate_candidates_order_4
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_quadrature, only: gauss_legendre_ab
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

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dgesv
    end interface

contains

    subroutine initialize_tetra_nedelec_first_kind(order, basis, status)
        integer, intent(in) :: order
        type(tetra_nedelec_first_kind_t), intent(out) :: basis
        integer, intent(out) :: status

        real(dp), allocatable :: moment_matrix(:, :)
        integer, allocatable :: pivots(:)
        integer :: candidate, info

        status = 1
        if (order < 1 .or. order > 4) return

        basis%order = order
        basis%dof_count = order * (order + 2) * (order + 3) / 2
        allocate(basis%coefficients(basis%dof_count, basis%dof_count))
        allocate(moment_matrix(basis%dof_count, basis%dof_count))
        allocate(pivots(basis%dof_count))

        call build_nedelec_moment_matrix(basis, moment_matrix, status)
        if (status /= 0) return
        basis%coefficients = 0.0_dp
        do candidate = 1, basis%dof_count
            basis%coefficients(candidate, candidate) = 1.0_dp
        end do
        call dgesv( &
            basis%dof_count, basis%dof_count, moment_matrix, basis%dof_count, &
            pivots, basis%coefficients, basis%dof_count, info)
        if (info /= 0) then
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

    subroutine build_nedelec_moment_matrix(basis, matrix, status)
        type(tetra_nedelec_first_kind_t), intent(in) :: basis
        real(dp), intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: candidate_curls(:, :)
        real(dp), allocatable :: candidate_values(:, :)
        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: face_weights(:), tetra_weights(:)
        real(dp), allocatable :: u(:), v(:), x(:), y(:), z(:)
        real(dp) :: point(3), tangent(3), tangents(3, 2)
        integer :: candidate, component, edge, exponent, face, moment, node
        integer :: total_degree, x_degree, y_degree, z_degree

        matrix = 0.0_dp
        status = 1
        allocate( &
            candidate_values(3, basis%dof_count), &
            candidate_curls(3, basis%dof_count))
        allocate( &
            edge_nodes(basis%order + 1), edge_weights(basis%order + 1))
        call gauss_legendre_ab( &
            basis%order + 1, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)

        moment = 0
        do edge = 1, 6
            do exponent = 0, basis%order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    call evaluate_nedelec_candidates( &
                        basis, point, candidate_values, candidate_curls)
                    do candidate = 1, basis%dof_count
                        matrix(moment, candidate) = &
                            matrix(moment, candidate) + &
                            edge_weights(node) * &
                            shifted_legendre(exponent, edge_nodes(node)) * &
                            dot_product(candidate_values(:, candidate), tangent)
                    end do
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * basis%order, u, v, face_weights, status)
        if (status /= 0) return
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, basis%order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(u)
                            call reference_face( &
                                face, u(node), v(node), point, tangents)
                            call evaluate_nedelec_candidates( &
                                basis, point, candidate_values, candidate_curls)
                            do candidate = 1, basis%dof_count
                                matrix(moment, candidate) = &
                                    matrix(moment, candidate) + &
                                    face_weights(node) * u(node)**x_degree * &
                                    v(node)**y_degree * dot_product( &
                                    candidate_values(:, candidate), &
                                    tangents(:, component))
                            end do
                        end do
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2 * basis%order, x, y, z, tetra_weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total_degree = 0, basis%order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            point = [x(node), y(node), z(node)]
                            call evaluate_nedelec_candidates( &
                                basis, point, candidate_values, candidate_curls)
                            matrix(moment, :) = matrix(moment, :) + &
                                tetra_weights(node) * x(node)**x_degree * &
                                y(node)**y_degree * z(node)**z_degree * &
                                candidate_values(component, :)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= basis%dof_count) then
            status = 2
            return
        end if
        status = 0
    end subroutine build_nedelec_moment_matrix

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
            values = 0.0_dp
            curls = 0.0_dp
        end select
    end subroutine evaluate_nedelec_candidates

    pure subroutine reference_edge(edge, parameter, point, tangent)
        integer, intent(in) :: edge
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: point(3), tangent(3)

        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), first, second

        call reference_topology(vertices, edge_vertices)
        first = edge_vertices(1, edge)
        second = edge_vertices(2, edge)
        point = (1.0_dp - parameter) * vertices(:, first) + &
            parameter * vertices(:, second)
        tangent = vertices(:, second) - vertices(:, first)
    end subroutine reference_edge

    pure subroutine reference_face(face, u, v, point, tangents)
        integer, intent(in) :: face
        real(dp), intent(in) :: u, v
        real(dp), intent(out) :: point(3), tangents(3, 2)

        real(dp) :: vertices(3, 4)
        integer :: edge_vertices(2, 6), face_vertices(3, 4)
        integer :: first, second, third

        call reference_topology(vertices, edge_vertices)
        face_vertices(:, 1) = [1, 2, 3]
        face_vertices(:, 2) = [1, 2, 4]
        face_vertices(:, 3) = [1, 3, 4]
        face_vertices(:, 4) = [2, 3, 4]
        first = face_vertices(1, face)
        second = face_vertices(2, face)
        third = face_vertices(3, face)
        tangents(:, 1) = vertices(:, second) - vertices(:, first)
        tangents(:, 2) = vertices(:, third) - vertices(:, first)
        point = vertices(:, first) + u * tangents(:, 1) + &
            v * tangents(:, 2)
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

    pure function shifted_legendre(degree, parameter) result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: parameter
        real(dp) :: value

        real(dp) :: current, previous, coordinate
        integer :: polynomial_degree

        coordinate = 2.0_dp * parameter - 1.0_dp
        if (degree == 0) then
            value = 1.0_dp
            return
        end if
        previous = 1.0_dp
        current = coordinate
        do polynomial_degree = 1, degree - 1
            value = ( &
                real(2 * polynomial_degree + 1, dp) * coordinate * current - &
                real(polynomial_degree, dp) * previous) / &
                real(polynomial_degree + 1, dp)
            previous = current
            current = value
        end do
        value = current
    end function shifted_legendre

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
