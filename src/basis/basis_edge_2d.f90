module fortfem_basis_edge_2d
    use fortfem_kinds
    use fortfem_mesh_2d
    implicit none
    private

    public :: edge_basis_2d_t
    public :: evaluate_edge_basis_2d
    public :: evaluate_edge_basis_curl_2d
    public :: evaluate_edge_basis_div_2d
    public :: evaluate_edge_basis_2d_piola

    ! Edge element basis functions (Nédélec elements)
    type :: edge_basis_2d_t
        integer :: n_edges = 0
        real(dp), allocatable :: edge_vectors(:,:) ! Direction vectors
        real(dp), allocatable :: edge_lengths(:) ! Edge lengths
    contains
        procedure :: init => edge_basis_init
        procedure :: destroy => edge_basis_destroy
        procedure :: n_dofs => edge_basis_n_dofs
    end type edge_basis_2d_t

contains

    subroutine edge_basis_init(this, mesh)
        class(edge_basis_2d_t), intent(inout) :: this
        type(mesh_2d_t), intent(in) :: mesh
        integer :: i, j, k, edge_count
        real(dp) :: dx, dy

        ! Count edges (3 per triangle)
        this%n_edges = 3 * mesh%n_triangles

        allocate(this%edge_vectors(2, this%n_edges))
        allocate(this%edge_lengths(this%n_edges))

        edge_count = 0

        ! For each triangle, compute edge vectors
        do i = 1, mesh%n_triangles
            do j = 1, 3
                edge_count = edge_count + 1

                ! Next vertex (cyclic)
                k = mod(j, 3) + 1

                ! Edge vector
                dx = mesh%vertices(1, mesh%triangles(k, i)) - &
                    mesh%vertices(1, mesh%triangles(j, i))
                dy = mesh%vertices(2, mesh%triangles(k, i)) - &
                    mesh%vertices(2, mesh%triangles(j, i))

                this%edge_vectors(1, edge_count) = dx
                this%edge_vectors(2, edge_count) = dy
                this%edge_lengths(edge_count) = sqrt(dx**2 + dy**2)
            end do
        end do
    end subroutine edge_basis_init

    subroutine edge_basis_destroy(this)
        class(edge_basis_2d_t), intent(inout) :: this
        if (allocated(this%edge_vectors)) deallocate(this%edge_vectors)
        if (allocated(this%edge_lengths)) deallocate(this%edge_lengths)
        this%n_edges = 0
    end subroutine edge_basis_destroy

    function edge_basis_n_dofs(this) result(n)
        class(edge_basis_2d_t), intent(in) :: this
        integer :: n
        n = this%n_edges
    end function edge_basis_n_dofs

    ! Evaluate edge basis functions at reference coordinates
    subroutine evaluate_edge_basis_2d(xi, eta, triangle_area, values)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: values(2, 3) ! 2D vectors, 3 edges

        ! RT0Ortho/Nédélec elements on reference triangle (0,0)-(1,0)-(0,1)
        ! Based on FreeFEM's RT0Ortho implementation which is RT0 rotated by 90°
        ! RT0 basis functions rotated: (x,y) → (-y,x)
        ! This gives H(curl) conforming edge elements

        real(dp) :: rt0_values(2, 3)

        ! First compute RT0 basis functions
        ! RT0 edge 0: φ₀ = (1-eta, 0)/(2*area)
        ! RT0 edge 1: φ₁ = (xi, eta-1)/(2*area)
        ! RT0 edge 2: φ₂ = (-xi, 1-eta)/(2*area)

        ! Normalize by triangle area for proper scaling
        real(dp) :: scale_factor
        scale_factor = 1.0_dp / (2.0_dp * triangle_area)

        rt0_values(1, 1) = (1.0_dp - eta) * scale_factor
        rt0_values(2, 1) = 0.0_dp

        rt0_values(1, 2) = xi * scale_factor
        rt0_values(2, 2) = (eta - 1.0_dp) * scale_factor

        rt0_values(1, 3) = -xi * scale_factor
        rt0_values(2, 3) = (1.0_dp - eta) * scale_factor

        ! Apply 90° rotation: (x,y) → (-y,x) to get RT0Ortho
        values(1, 1) = -rt0_values(2, 1) ! -0 = 0
        values(2, 1) = rt0_values(1, 1) ! (1-eta)/2A

        values(1, 2) = -rt0_values(2, 2) ! -(eta-1)/2A = (1-eta)/2A
        values(2, 2) = rt0_values(1, 2) ! xi/2A

        values(1, 3) = -rt0_values(2, 3) ! -(1-eta)/2A = (eta-1)/2A
        values(2, 3) = rt0_values(1, 3) ! -xi/2A
    end subroutine evaluate_edge_basis_2d

    ! Evaluate curl of edge basis functions
    subroutine evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: curls(3) ! Scalar curl in 2D

        ! For 2D vector field φ = (φˣ, φʸ), curl(φ) = ∂φʸ/∂ξ - ∂φˣ/∂η
        ! Then transform from reference to physical: curl_phys = curl_ref / jacobian_det
        ! where jacobian_det = 2 * triangle_area for linear triangular mapping

        real(dp) :: jacobian_det, scale_factor
        jacobian_det = 2.0_dp * triangle_area
        scale_factor = 1.0_dp / (2.0_dp * triangle_area)

        ! For RT0Ortho basis functions (RT0 rotated by 90°):
        ! φ₀ = (0, (1-eta)/(2A)): curl = ∂((1-eta)/(2A))/∂ξ - ∂(0)/∂η = 0 - 0 = 0
        ! φ₁ = ((1-eta)/(2A), xi/(2A)): curl = ∂(xi/(2A))/∂ξ - ∂((1-eta)/(2A))/∂η
        !      = 1/(2A) - (-1/(2A)) = 2/(2A) = 1/A
        ! φ₂ = ((eta-1)/(2A), -xi/(2A)): curl = ∂(-xi/(2A))/∂ξ - ∂((eta-1)/(2A))/∂η
        !      = -1/(2A) - 1/(2A) = -2/(2A) = -1/A

        ! Transform to physical element (already in physical coordinates)
        curls(1) = 0.0_dp
        curls(2) = 1.0_dp / triangle_area
        curls(3) = -1.0_dp / triangle_area
    end subroutine evaluate_edge_basis_curl_2d

    ! Evaluate edge basis functions with Piola transformation
    subroutine evaluate_edge_basis_2d_piola(mesh, triangle_idx, xi, eta, values)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: values(2, 3) ! 2D vectors, 3 edges

        real(dp) :: ref_values(2, 3) ! Reference basis values
        real(dp) :: jacobian(2, 2), inv_jacobian(2, 2), det_jacobian
        real(dp) :: triangle_area
        integer :: i, j, k
        real(dp) :: x1, y1, x2, y2, x3, y3

        ! Get triangle vertices
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))

        ! Compute Jacobian matrix: F = [∂x/∂ξ, ∂x/∂η; ∂y/∂ξ, ∂y/∂η]
        jacobian(1, 1) = x2 - x1 ! ∂x/∂ξ
        jacobian(1, 2) = x3 - x1 ! ∂x/∂η
        jacobian(2, 1) = y2 - y1 ! ∂y/∂ξ
        jacobian(2, 2) = y3 - y1 ! ∂y/∂η

        det_jacobian = jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1)
        triangle_area = 0.5_dp * abs(det_jacobian)

        ! Evaluate reference basis functions
        call evaluate_edge_basis_2d(xi, eta, triangle_area, ref_values)

        ! Apply Piola transformation: φ_phys = (1/J) * F * φ_ref
        do i = 1, 3
            values(1, i) = (jacobian(1, 1) * ref_values(1, i) + jacobian(1, 2) * ref_values(2, i)) / det_jacobian
            values(2, i) = (jacobian(2, 1) * ref_values(1, i) + jacobian(2, 2) * ref_values(2, i)) / det_jacobian
        end do
    end subroutine evaluate_edge_basis_2d_piola

    ! Evaluate divergence of edge basis functions
    subroutine evaluate_edge_basis_div_2d(xi, eta, triangle_area, divs)
        real(dp), intent(in) :: xi, eta, triangle_area
        real(dp), intent(out) :: divs(3) ! Divergence values

        ! Divergence of Nédélec elements (should be zero for H(curl))
        divs(1) = 0.0_dp
        divs(2) = 0.0_dp
        divs(3) = 0.0_dp
    end subroutine evaluate_edge_basis_div_2d

end module fortfem_basis_edge_2d
