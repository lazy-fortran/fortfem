module fortfem_tetra_rt_interpolation
    use fortfem_kinds, only: dp
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortfem_tetra_rt_arbitrary_order, only: &
        tetra_rt_dof_count, tetra_rt_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    private

    public :: interpolate_reference_tetra_rt

    abstract interface
        subroutine reference_vector_field_3d(point, value)
            import :: dp
            real(dp), intent(in) :: point(3)
            real(dp), intent(out) :: value(3)
        end subroutine reference_vector_field_3d
    end interface

contains

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
                            x(node)**x_degree*y(node)**y_degree
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
                                field_value(component)*x(node)**x_degree* &
                                y(node)**y_degree*z(node)**z_degree
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

        do degree = 0, 4
            if (dof_count == &
                (degree + 1)*(degree + 2)*(degree + 4)/2) return
        end do
        degree = -1
    end function degree_from_dof_count

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
