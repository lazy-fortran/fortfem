module fortfem_tetra_rt_interpolation
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_rt_arbitrary_order, only: &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortnum_linalg, only: det3, inv3
    use fortnum_special_jacobi, only: tetrahedron_koornwinder, &
        triangle_dubiner
    implicit none

    private

    public :: interpolate_reference_tetra_rt
    public :: interpolate_physical_tetra_rt

    abstract interface
        pure subroutine physical_vector_field_3d(x, y, z, value)
            import :: dp
            real(dp), intent(in) :: x, y, z
            real(dp), intent(out) :: value(3)
        end subroutine physical_vector_field_3d

        subroutine reference_vector_field_3d(point, value)
            import :: dp
            real(dp), intent(in) :: point(3)
            real(dp), intent(out) :: value(3)
        end subroutine reference_vector_field_3d
    end interface

contains

    subroutine interpolate_physical_tetra_rt( &
            basis, vertices, field, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        procedure(physical_vector_field_3d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp) :: determinant, inverse_jacobian(3, 3), jacobian(3, 3)
        integer :: inverse_status

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        status = 1
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call inv3(jacobian, inverse_jacobian, inverse_status)
        if (inverse_status /= 0) return
        call interpolate_reference_tetra_rt( &
            basis, pulled_back_field, dofs, status)

    contains

        subroutine pulled_back_field(reference_point, value)
            real(dp), intent(in) :: reference_point(3)
            real(dp), intent(out) :: value(3)

            real(dp) :: physical_point(3), physical_value(3)

            physical_point = vertices(:, 1) + &
                matmul(jacobian, reference_point)
            call field( &
                physical_point(1), physical_point(2), physical_point(3), &
                physical_value)
            value = determinant*matmul(inverse_jacobian, physical_value)
        end subroutine pulled_back_field

    end subroutine interpolate_physical_tetra_rt

    subroutine interpolate_reference_tetra_rt(basis, field, dofs, status)
        type(tetra_rt_t), intent(in) :: basis
        procedure(reference_vector_field_3d) :: field
        real(dp), allocatable, intent(out) :: dofs(:)
        integer, intent(out) :: status

        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: area_normal(3), field_value(3), point(3)
        integer :: component, degree, face, moment, node, total
        integer :: x_degree, y_degree, z_degree

        status = 1
        degree = degree_from_dof_count(tetra_rt_dof_count(basis))
        if (degree < 0) return
        allocate(dofs(tetra_rt_dof_count(basis)))
        dofs = 0.0_dp
        moment = 0
        call triangle_duffy_quadrature( &
            2*degree + 4, x, y, weights, status)
        if (status /= 0) return
        do face = 1, 4
            do total = 0, degree
                do x_degree = 0, total
                    y_degree = total - x_degree
                    moment = moment + 1
                    do node = 1, size(weights)
                        call reference_face( &
                            face, x(node), y(node), point, area_normal)
                        call field(point, field_value)
                        dofs(moment) = dofs(moment) + weights(node)* &
                            dot_product(field_value, area_normal)* &
                            face_moment_value( &
                            degree, x_degree, y_degree, &
                            x(node), y(node))
                    end do
                end do
            end do
        end do

        call tetra_duffy_quadrature( &
            2*degree + 4, x, y, z, weights, status)
        if (status /= 0) return
        do component = 1, 3
            do total = 0, degree - 1
                do x_degree = 0, total
                    do y_degree = 0, total - x_degree
                        z_degree = total - x_degree - y_degree
                        moment = moment + 1
                        do node = 1, size(weights)
                            point = [x(node), y(node), z(node)]
                            call field(point, field_value)
                            dofs(moment) = dofs(moment) + weights(node)* &
                                field_value(component)*volume_moment_value( &
                                degree, x_degree, y_degree, z_degree, &
                                point)
                        end do
                    end do
                end do
            end do
        end do
        if (moment /= size(dofs)) return
        status = 0
    end subroutine interpolate_reference_tetra_rt

    pure function degree_from_dof_count(dof_count) result(degree)
        integer, intent(in) :: dof_count
        integer :: degree

        do degree = 0, 1000
            if (dof_count == &
                (degree + 1)*(degree + 2)*(degree + 4)/2) return
            if ((degree + 1)*(degree + 2)*(degree + 4)/2 > dof_count) exit
        end do
        degree = -1
    end function degree_from_dof_count

    pure real(dp) function face_moment_value( &
            degree, first_degree, second_degree, x, y) result(value)
        integer, intent(in) :: degree, first_degree, second_degree
        real(dp), intent(in) :: x, y

        if (degree <= 5) then
            value = x**first_degree*y**second_degree
        else
            value = triangle_dubiner(first_degree, second_degree, x, y)
        end if
    end function face_moment_value

    pure real(dp) function volume_moment_value( &
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
    end function volume_moment_value

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

end module fortfem_tetra_rt_interpolation
