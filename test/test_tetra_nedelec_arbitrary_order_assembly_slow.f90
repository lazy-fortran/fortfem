program test_tetra_nedelec_arbitrary_order_assembly_slow
    use check, only: check_condition, check_summary
    use fortfem_api, only: assemble_tetra_nedelec_curl_mass_element, &
        tetra_duffy_quadrature, triangle_duffy_quadrature
    use fortfem_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    real(dp), parameter :: curl_coefficient = 1.7_dp
    real(dp), parameter :: mass_coefficient = 0.8_dp
    real(dp), allocatable :: dofs(:), matrix(:, :)
    real(dp) :: energy, exact_energy, vertices(3, 4)
    integer :: dof_count, order, status
    logical :: all_passed

    all_passed = .true.
    vertices(:, 1) = [0.2_dp, -0.3_dp, 0.4_dp]
    vertices(:, 2) = vertices(:, 1) + [2.0_dp, 0.0_dp, 0.0_dp]
    vertices(:, 3) = vertices(:, 1) + [0.0_dp, 3.0_dp, 0.0_dp]
    vertices(:, 4) = vertices(:, 1) + [0.0_dp, 0.0_dp, 4.0_dp]
    do order = 1, 4
        dof_count = order * (order + 2) * (order + 3) / 2
        allocate(dofs(dof_count))
        call reference_field_dofs(order, dofs)
        call assemble_tetra_nedelec_curl_mass_element( &
            vertices, order, 2 * order + 2, matrix, status, &
            curl_coefficient, mass_coefficient)
        energy = dot_product(dofs, matmul(matrix, dofs))
        exact_energy = exact_mapped_energy(order)
        call record_condition(status == 0 .and. &
            abs(energy - exact_energy) < 3.0e-8_dp, &
            "Higher-order tetrahedral assembly has exact polynomial energy")
        call record_condition(maxval(abs(matrix - transpose(matrix))) < &
            3.0e-12_dp, &
            "Higher-order tetrahedral curl-mass matrix is symmetric")
        deallocate(dofs, matrix)
    end do

    call assemble_tetra_nedelec_curl_mass_element( &
        vertices, 0, 2, matrix, status)
    call record_condition(status /= 0 .and. .not. allocated(matrix), &
        "Higher-order tetrahedral assembly rejects order zero")

    call check_summary("Higher-order tetrahedral Nedelec assembly")
    if (.not. all_passed) error stop 1

contains

    subroutine reference_field_dofs(order, dofs)
        integer, intent(in) :: order
        real(dp), intent(out) :: dofs(:)

        real(dp), allocatable :: edge_nodes(:), edge_weights(:)
        real(dp), allocatable :: tetra_weights(:), triangle_weights(:)
        real(dp), allocatable :: x(:), y(:), z(:)
        real(dp) :: field(3), point(3), tangent(3), tangents(3, 2)
        integer :: component, edge, exponent, face, moment, node, status
        integer :: total_degree, x_degree, y_degree, z_degree

        allocate(edge_nodes(order + 2), edge_weights(order + 2))
        call gauss_legendre_ab( &
            order + 2, 0.0_dp, 1.0_dp, edge_nodes, edge_weights)
        dofs = 0.0_dp
        moment = 0
        do edge = 1, 6
            do exponent = 0, order - 1
                moment = moment + 1
                do node = 1, size(edge_nodes)
                    call reference_edge( &
                        edge, edge_nodes(node), point, tangent)
                    field = reference_field(order, point)
                    dofs(moment) = dofs(moment) + edge_weights(node) * &
                        shifted_legendre(exponent, edge_nodes(node)) * &
                        dot_product(field, tangent)
                end do
            end do
        end do

        call triangle_duffy_quadrature( &
            2 * order + 2, x, y, triangle_weights, status)
        do face = 1, 4
            do component = 1, 2
                do total_degree = 0, order - 2
                    do x_degree = 0, total_degree
                        y_degree = total_degree - x_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            call reference_face( &
                                face, x(node), y(node), point, tangents)
                            field = reference_field(order, point)
                            dofs(moment) = dofs(moment) + &
                                triangle_weights(node) * &
                                x(node)**x_degree * y(node)**y_degree * &
                                dot_product(field, tangents(:, component))
                        end do
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2 * order + 2, x, y, z, tetra_weights, status)
        do component = 1, 3
            do total_degree = 0, order - 3
                do x_degree = 0, total_degree
                    do y_degree = 0, total_degree - x_degree
                        z_degree = total_degree - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(x)
                            point = [x(node), y(node), z(node)]
                            field = reference_field(order, point)
                            dofs(moment) = dofs(moment) + &
                                tetra_weights(node) * x(node)**x_degree * &
                                y(node)**y_degree * z(node)**z_degree * &
                                field(component)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine reference_field_dofs

    pure function exact_mapped_energy(order) result(value)
        integer, intent(in) :: order
        real(dp) :: value

        value = mass_coefficient * (24.0_dp / 9.0_dp) * &
            monomial_integral(2 * order - 2)
        if (order > 1) then
            value = value + curl_coefficient * (2.0_dp / 3.0_dp) * &
                real((order - 1)**2, dp) * &
                monomial_integral(2 * order - 4)
        end if
    end function exact_mapped_energy

    pure function monomial_integral(degree) result(value)
        integer, intent(in) :: degree
        real(dp) :: value

        value = 1.0_dp / real( &
            (degree + 1) * (degree + 2) * (degree + 3), dp)
    end function monomial_integral

    pure function reference_field(order, point) result(field)
        integer, intent(in) :: order
        real(dp), intent(in) :: point(3)
        real(dp) :: field(3)

        field = [0.0_dp, point(1)**(order - 1), 0.0_dp]
    end function reference_field

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

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_tetra_nedelec_arbitrary_order_assembly_slow
