module fortfem_tetra_affine_map
    !! Analytical physical-to-reference products for an affine tetrahedron.
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: det3
    use fortfem_generated_tetra_affine, only: &
        generated_tetra_affine
    use fortfem_generated_tetra_affine_jvp, only: &
        generated_tetra_affine_jvp
    use fortfem_generated_tetra_affine_vjp, only: &
        generated_tetra_affine_vjp
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

        real(dp) :: jacobian(3, 3)
        real(dp) :: relative(3)

        reference = 0.0_dp
        include "../generated/fortfem_tetra_geometry.inc"
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call generated_tetra_affine(jacobian, relative, reference)
        status = 0
    end subroutine invert_tetra_affine_map

    pure subroutine invert_tetra_affine_map_jvp( &
            vertices, point, vertices_dot, point_dot, reference_dot, status)
        real(dp), intent(in) :: vertices(3, 4), point(3)
        real(dp), intent(in) :: vertices_dot(3, 4), point_dot(3)
        real(dp), intent(out) :: reference_dot(3)
        integer, intent(out) :: status

        real(dp) :: jacobian(3, 3)
        real(dp) :: jacobian_dot(3, 3)
        real(dp) :: relative(3), relative_dot(3)

        reference_dot = 0.0_dp
        include "../generated/fortfem_tetra_geometry.inc"
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        include "../generated/fortfem_tetra_geometry_jvp.inc"
        call generated_tetra_affine_jvp( &
            jacobian, relative, jacobian_dot, relative_dot, reference_dot)
        status = 0
    end subroutine invert_tetra_affine_map_jvp

    pure subroutine invert_tetra_affine_map_vjp( &
            vertices, point, reference_bar, vertices_bar, point_bar, status)
        real(dp), intent(in) :: vertices(3, 4), point(3), reference_bar(3)
        real(dp), intent(out) :: vertices_bar(3, 4), point_bar(3)
        integer, intent(out) :: status

        real(dp) :: jacobian(3, 3)
        real(dp) :: jacobian_bar(3, 3)
        real(dp) :: relative(3), relative_bar(3)

        vertices_bar = 0.0_dp
        point_bar = 0.0_dp
        include "../generated/fortfem_tetra_geometry.inc"
        if (.not. valid_jacobian(jacobian)) then
            status = 1
            return
        end if
        call generated_tetra_affine_vjp( &
            jacobian, relative, reference_bar, jacobian_bar, relative_bar)
        include "../generated/fortfem_tetra_geometry_vjp.inc"
        status = 0
    end subroutine invert_tetra_affine_map_vjp


    pure logical function valid_jacobian(jacobian) result(valid)
        real(dp), intent(in) :: jacobian(3, 3)
        real(dp) :: determinant

        determinant = det3(jacobian)
        valid = determinant > 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, maxval(abs(jacobian))**3)
    end function valid_jacobian

end module fortfem_tetra_affine_map
