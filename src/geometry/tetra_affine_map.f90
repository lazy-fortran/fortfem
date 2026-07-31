module fortfem_tetra_affine_map
    !! Analytical physical-to-reference products for an affine tetrahedron.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3, inv3, inv3_jvp, inv3_vjp
    implicit none

    private

    public :: invert_tetra_affine_map
    public :: invert_tetra_affine_map_jvp
    public :: invert_tetra_affine_map_vjp

contains

    pure subroutine invert_tetra_affine_map( &
            vertices, point, reference, status)
        real(dp), intent(in) :: vertices(3, 4), point(3)
        real(dp), intent(out) :: reference(3)
        integer, intent(out) :: status

        real(dp) :: inverse(3, 3), jacobian(3, 3)

        reference = 0.0_dp
        call tetrahedron_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call inv3(jacobian, inverse, status)
        if (status /= 0) return
        reference = matmul(inverse, point - vertices(:, 1))
    end subroutine invert_tetra_affine_map

    pure subroutine invert_tetra_affine_map_jvp( &
            vertices, point, vertices_dot, point_dot, reference_dot, status)
        real(dp), intent(in) :: vertices(3, 4), point(3)
        real(dp), intent(in) :: vertices_dot(3, 4), point_dot(3)
        real(dp), intent(out) :: reference_dot(3)
        integer, intent(out) :: status

        real(dp) :: inverse(3, 3), inverse_dot(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_dot(3, 3)
        real(dp) :: relative(3), relative_dot(3)

        reference_dot = 0.0_dp
        call tetrahedron_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call tetrahedron_jacobian(vertices_dot, jacobian_dot)
        call inv3_jvp( &
            jacobian, jacobian_dot, inverse, inverse_dot, status)
        if (status /= 0) return
        relative = point - vertices(:, 1)
        relative_dot = point_dot - vertices_dot(:, 1)
        reference_dot = matmul(inverse_dot, relative) + &
            matmul(inverse, relative_dot)
    end subroutine invert_tetra_affine_map_jvp

    pure subroutine invert_tetra_affine_map_vjp( &
            vertices, point, reference_bar, vertices_bar, point_bar, status)
        real(dp), intent(in) :: vertices(3, 4), point(3), reference_bar(3)
        real(dp), intent(out) :: vertices_bar(3, 4), point_bar(3)
        integer, intent(out) :: status

        real(dp) :: inverse(3, 3), inverse_bar(3, 3)
        real(dp) :: jacobian(3, 3), jacobian_bar(3, 3)
        real(dp) :: relative(3)

        vertices_bar = 0.0_dp
        point_bar = 0.0_dp
        call tetrahedron_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        relative = point - vertices(:, 1)
        inverse_bar = spread(reference_bar, 2, 3)*spread(relative, 1, 3)
        call inv3_vjp( &
            jacobian, inverse_bar, inverse, jacobian_bar, status)
        if (status /= 0) return
        point_bar = matmul(transpose(inverse), reference_bar)
        vertices_bar(:, 1) = -point_bar - jacobian_bar(:, 1) - &
            jacobian_bar(:, 2) - jacobian_bar(:, 3)
        vertices_bar(:, 2) = jacobian_bar(:, 1)
        vertices_bar(:, 3) = jacobian_bar(:, 2)
        vertices_bar(:, 4) = jacobian_bar(:, 3)
    end subroutine invert_tetra_affine_map_vjp

    pure subroutine tetrahedron_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(3, 4)
        real(dp), intent(out) :: jacobian(3, 3)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
        jacobian(:, 3) = vertices(:, 4) - vertices(:, 1)
    end subroutine tetrahedron_jacobian

    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp) :: determinant

        determinant = det3(jacobian)
        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
    end function valid_jacobian

end module fortfem_tetra_affine_map
