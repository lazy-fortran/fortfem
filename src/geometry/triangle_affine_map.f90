module fortfem_triangle_affine_map
    !! Analytical physical-to-reference products for an affine triangle.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det2
    use fortfem_generated_triangle_affine, only: &
        generated_triangle_affine
    use fortfem_generated_triangle_affine_jvp, only: &
        generated_triangle_affine_jvp
    use fortfem_generated_triangle_affine_vjp, only: &
        generated_triangle_affine_vjp
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

        real(dp) :: jacobian(2, 2)
        real(dp) :: relative(2)

        reference = 0.0_dp
        include "../generated/fortfem_triangle_geometry.inc"
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call generated_triangle_affine(jacobian, relative, reference)
        status = 0
    end subroutine invert_triangle_affine_map

    pure subroutine invert_triangle_affine_map_jvp( &
            vertices, point, vertices_dot, point_dot, reference_dot, status)
        real(dp), intent(in) :: vertices(2, 3), point(2)
        real(dp), intent(in) :: vertices_dot(2, 3), point_dot(2)
        real(dp), intent(out) :: reference_dot(2)
        integer, intent(out) :: status

        real(dp) :: jacobian(2, 2)
        real(dp) :: jacobian_dot(2, 2)
        real(dp) :: relative(2), relative_dot(2)

        reference_dot = 0.0_dp
        include "../generated/fortfem_triangle_geometry.inc"
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        include "../generated/fortfem_triangle_geometry_jvp.inc"
        call generated_triangle_affine_jvp( &
            jacobian, relative, jacobian_dot, relative_dot, reference_dot)
        status = 0
    end subroutine invert_triangle_affine_map_jvp

    pure subroutine invert_triangle_affine_map_vjp( &
            vertices, point, reference_bar, vertices_bar, point_bar, status)
        real(dp), intent(in) :: vertices(2, 3), point(2), reference_bar(2)
        real(dp), intent(out) :: vertices_bar(2, 3), point_bar(2)
        integer, intent(out) :: status

        real(dp) :: jacobian(2, 2)
        real(dp) :: jacobian_bar(2, 2)
        real(dp) :: relative(2), relative_bar(2)

        vertices_bar = 0.0_dp
        point_bar = 0.0_dp
        include "../generated/fortfem_triangle_geometry.inc"
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call generated_triangle_affine_vjp( &
            jacobian, relative, reference_bar, jacobian_bar, relative_bar)
        include "../generated/fortfem_triangle_geometry_vjp.inc"
        status = 0
    end subroutine invert_triangle_affine_map_vjp


    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(2, 2)
        real(dp) :: determinant

        determinant = det2(jacobian)
        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**2)
    end function valid_jacobian

end module fortfem_triangle_affine_map
