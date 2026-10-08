module fortfem_basis_rt_2d
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_generated_rt_basis, only: generated_rt_jet
    implicit none
    private

    public :: evaluate_rt_basis_2d
    public :: evaluate_rt_basis_div_2d
    public :: evaluate_rt_basis_2d_piola

contains

    pure subroutine evaluate_rt_basis_2d(xi, eta, triangle_area, values)
        real(dp), intent(in) :: xi,eta,triangle_area
        real(dp), intent(out) :: values(2,3)
        real(dp) :: divergence,curl
        integer :: i
        if (triangle_area <= 0.0_dp) error stop "triangle area must be positive"
        do i=1,3
            call generated_rt_jet(i,xi,eta,values(:,i),divergence,curl)
        end do
    end subroutine evaluate_rt_basis_2d

    pure subroutine evaluate_rt_basis_div_2d(xi, eta, triangle_area, divergences)
        real(dp), intent(in) :: xi,eta,triangle_area
        real(dp), intent(out) :: divergences(3)
        real(dp) :: v(2),divergence,curl
        integer :: i
        if (triangle_area <= 0.0_dp) error stop "triangle area must be positive"
        do i=1,3
            call generated_rt_jet(i,xi,eta,v,divergence,curl)
            divergences(i) = divergence/(2*triangle_area)
        end do
    end subroutine evaluate_rt_basis_div_2d

    subroutine evaluate_rt_basis_2d_piola( &
            mesh, triangle_idx, xi, eta, values)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3)

        real(dp) :: ref_values(2, 3)
        real(dp) :: jacobian(2, 2), det_jacobian
        real(dp) :: x1, y1, x2, y2, x3, y3
        integer :: basis_id

        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))

        jacobian(:, 1) = [x2 - x1, y2 - y1]
        jacobian(:, 2) = [x3 - x1, y3 - y1]
        det_jacobian = jacobian(1, 1) * jacobian(2, 2) - &
            jacobian(1, 2) * jacobian(2, 1)
        if (det_jacobian <= 0.0_dp) then
            error stop &
                "evaluate_rt_basis_2d_piola: triangle must be counter-clockwise"
        end if

        call evaluate_rt_basis_2d(xi, eta, 0.5_dp, ref_values)

        ! Contravariant Piola transformation: phi_phys = J*phi_ref/det(J).
        do basis_id = 1, 3
            values(:, basis_id) = &
                matmul(jacobian, ref_values(:, basis_id)) / det_jacobian
        end do
    end subroutine evaluate_rt_basis_2d_piola

end module fortfem_basis_rt_2d
