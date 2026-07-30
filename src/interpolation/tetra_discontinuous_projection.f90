module fortfem_tetra_discontinuous_projection
    use fortfem_kinds, only: dp
    use fortfem_tetra_discontinuous_arbitrary_order, only: &
        evaluate_tetra_discontinuous, tetra_discontinuous_dof_count, &
        tetra_discontinuous_t
    use fortfem_tetra_duffy_quadrature, only: tetra_duffy_quadrature
    use fortnum_linalg, only: dense_solve, det3
    implicit none

    private

    public :: project_physical_tetra_discontinuous

    abstract interface
        pure real(dp) function physical_scalar_field_3d(x, y, z) result(value)
            import :: dp
            real(dp), intent(in) :: x, y, z
        end function physical_scalar_field_3d
    end interface

contains

    subroutine project_physical_tetra_discontinuous( &
            basis, vertices, field, coefficients, status)
        type(tetra_discontinuous_t), intent(in) :: basis
        real(dp), intent(in) :: vertices(3, 4)
        procedure(physical_scalar_field_3d) :: field
        real(dp), allocatable, intent(out) :: coefficients(:)
        integer, intent(out) :: status

        real(dp), allocatable :: mass(:, :), right_hand_side(:), values(:)
        real(dp), allocatable :: weights(:), x(:), y(:), z(:)
        real(dp) :: determinant, jacobian(3, 3), physical_point(3)
        integer :: column, degree, dof_count, node, row

        status = 1
        degree = degree_from_dof_count(tetra_discontinuous_dof_count(basis))
        if (degree < 0) return
        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
        determinant = det3(jacobian)
        if (determinant <= 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)) return
        call tetra_duffy_quadrature( &
            2*degree + 8, x, y, z, weights, status)
        if (status /= 0) return
        dof_count = tetra_discontinuous_dof_count(basis)
        allocate ( &
            mass(dof_count, dof_count), right_hand_side(dof_count), &
            values(dof_count), coefficients(dof_count))
        mass = 0.0_dp
        right_hand_side = 0.0_dp
        do node = 1, size(weights)
            call evaluate_tetra_discontinuous( &
                basis, x(node), y(node), z(node), values, status)
            if (status /= 0) return
            physical_point = vertices(:, 1) + &
                matmul(jacobian, [x(node), y(node), z(node)])
            do column = 1, dof_count
                right_hand_side(column) = right_hand_side(column) + &
                    weights(node)*values(column)*field( &
                    physical_point(1), physical_point(2), physical_point(3))
                do row = 1, dof_count
                    mass(row, column) = mass(row, column) + &
                        weights(node)*values(row)*values(column)
                end do
            end do
        end do
        call dense_solve(mass, right_hand_side, coefficients, status)
    end subroutine project_physical_tetra_discontinuous

    pure integer function degree_from_dof_count(dof_count) result(degree)
        integer, intent(in) :: dof_count

        do degree = 0, 4
            if (dof_count == (degree + 1)*(degree + 2)*(degree + 3)/6) return
        end do
        degree = -1
    end function degree_from_dof_count

end module fortfem_tetra_discontinuous_projection
