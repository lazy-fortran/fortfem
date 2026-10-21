module parametric_bem_test_geometry
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    use fortfem_sphere_curved_panel, only: evaluate_sphere_curved_panel
    use fortfem_torus_curved_panel, only: evaluate_torus_curved_panel
    implicit none
    type, extends(bem_surface_t) :: test_sphere_t
        real(dp), allocatable :: vertices(:, :)
        integer, allocatable :: faces(:, :)
    contains
        procedure :: evaluate => sphere_evaluate
    end type
    type, extends(bem_surface_t) :: test_torus_t
        real(dp), allocatable :: parameters(:, :)
        integer, allocatable :: faces(:, :)
        real(dp) :: deformation(3,3) = reshape( &
            [1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp], &
            [3,3])
    contains
        procedure :: evaluate => torus_evaluate
    end type
    interface assignment(=)
        module procedure assign_sphere, assign_torus
    end interface
contains
    subroutine assign_sphere(lhs,rhs)
        type(test_sphere_t), intent(out) :: lhs
        type(test_sphere_t), intent(in) :: rhs
        if (allocated(rhs%vertices)) lhs%vertices = rhs%vertices
        if (allocated(rhs%faces)) lhs%faces = rhs%faces
    end subroutine
    subroutine assign_torus(lhs,rhs)
        type(test_torus_t), intent(out) :: lhs
        type(test_torus_t), intent(in) :: rhs
        if (allocated(rhs%parameters)) lhs%parameters = rhs%parameters
        if (allocated(rhs%faces)) lhs%faces = rhs%faces
        lhs%deformation = rhs%deformation
    end subroutine
    pure subroutine torus_evaluate(self, panel, uv, x, a, b, status)
        class(test_torus_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        real(dp) :: jacobian, p(3), ta(3), tb(3), local(2,3)
        local = self%parameters(:, self%faces(:,panel))
        call evaluate_torus_curved_panel(local, &
            2.0_dp, 0.6_dp, uv(1), uv(2), p, ta, tb, jacobian, status)
        x = matmul(self%deformation, p)
        a = matmul(self%deformation, ta)
        b = matmul(self%deformation, tb)
    end subroutine

    pure subroutine sphere_evaluate(self, panel, uv, x, a, b, status)
        class(test_sphere_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3), a(3), b(3)
        integer, intent(out) :: status
        real(dp) :: jacobian, local(3,3)
        local = self%vertices(:, self%faces(:,panel))
        call evaluate_sphere_curved_panel(local, &
            1.0_dp, uv(1), uv(2), x, a, b, jacobian, status)
    end subroutine
end module

