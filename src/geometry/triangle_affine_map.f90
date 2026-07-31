module fortfem_triangle_affine_map
    !! Analytical physical-to-reference products for an affine triangle.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det2, inv2, inv2_jvp, inv2_vjp
    implicit none

    private

    public :: invert_triangle_affine_map
    public :: invert_triangle_affine_map_jvp
    public :: invert_triangle_affine_map_vjp

contains

    pure subroutine invert_triangle_affine_map( &
            vertices, point, reference, status)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(out) :: reference(2)
        integer, intent(out) :: status

        real(dp) :: inverse(2, 2), jacobian(2, 2)

        reference = 0.0_dp
        call triangle_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call inv2(jacobian, inverse, status)
        if (status /= 0) return
        reference = matmul(inverse, point - vertices(:, 1))
    end subroutine invert_triangle_affine_map

    pure subroutine invert_triangle_affine_map_jvp( &
            vertices, point, vertices_dot, point_dot, reference_dot, status)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(in) :: vertices_dot(2, 3), point_dot(2)
        real(dp), intent(out) :: reference_dot(2)
        integer, intent(out) :: status

        real(dp) :: inverse(2, 2), inverse_dot(2, 2)
        real(dp) :: jacobian(2, 2), jacobian_dot(2, 2)
        real(dp) :: relative(2), relative_dot(2)

        reference_dot = 0.0_dp
        call triangle_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call triangle_jacobian(vertices_dot, jacobian_dot)
        call inv2_jvp( &
            jacobian, jacobian_dot, inverse, inverse_dot, status)
        if (status /= 0) return
        relative = point - vertices(:, 1)
        relative_dot = point_dot - vertices_dot(:, 1)
        reference_dot = matmul(inverse_dot, relative) + &
            matmul(inverse, relative_dot)
    end subroutine invert_triangle_affine_map_jvp

    pure subroutine invert_triangle_affine_map_vjp( &
            vertices, point, reference_bar, vertices_bar, point_bar, status)
        real(dp), intent(in) :: vertices(2, 3), point(2), reference_bar(2)
        real(dp), intent(out) :: vertices_bar(2, 3), point_bar(2)
        integer, intent(out) :: status

        real(dp) :: inverse(2, 2), inverse_bar(2, 2)
        real(dp) :: jacobian(2, 2), jacobian_bar(2, 2)
        real(dp) :: relative(2)

        vertices_bar = 0.0_dp
        point_bar = 0.0_dp
        call triangle_jacobian(vertices, jacobian)
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        relative = point - vertices(:, 1)
        inverse_bar = spread(reference_bar, 2, 2)*spread(relative, 1, 2)
        call inv2_vjp( &
            jacobian, inverse_bar, inverse, jacobian_bar, status)
        if (status /= 0) return
        point_bar = matmul(transpose(inverse), reference_bar)
        vertices_bar(:, 1) = -point_bar - &
            jacobian_bar(:, 1) - jacobian_bar(:, 2)
        vertices_bar(:, 2) = jacobian_bar(:, 1)
        vertices_bar(:, 3) = jacobian_bar(:, 2)
    end subroutine invert_triangle_affine_map_vjp

    pure subroutine triangle_jacobian(vertices, jacobian)
        real(dp), intent(in) :: vertices(2, 3)
        real(dp), intent(out) :: jacobian(2, 2)

        jacobian(:, 1) = vertices(:, 2) - vertices(:, 1)
        jacobian(:, 2) = vertices(:, 3) - vertices(:, 1)
    end subroutine triangle_jacobian

    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp) :: determinant

        determinant = det2(jacobian)
        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**2)
    end function valid_jacobian

end module fortfem_triangle_affine_map
