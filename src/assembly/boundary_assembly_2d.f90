module boundary_assembly_2d_module
    use fortfem_kinds
    use basis_p1_2d_module
    use basis_p2_2d_module
    use fortfem_mesh_2d
    use fortfem_sparse_matrix, only: triplet_matrix_t
    implicit none
    private
    
    public :: boundary_assembly_2d_t
    public :: assemble_neumann_bc
    
    type :: boundary_assembly_2d_t
        integer :: element_order = 1  ! 1 for P1, 2 for P2
        integer :: boundary_marker = 0
    contains
        procedure :: init
        procedure :: destroy
        procedure :: element_boundary_vector_p1
        procedure :: element_boundary_vector_p2
        procedure :: element_boundary_matrix_p1
        procedure :: element_boundary_matrix_p2
    end type boundary_assembly_2d_t
    
    ! 1D Gauss quadrature for edge integrals
    integer, parameter :: n_gauss_1d_order2 = 2
    real(dp), parameter :: gauss_points_1d_order2(2) = [ &
        0.5_dp - sqrt(3.0_dp)/6.0_dp, &
        0.5_dp + sqrt(3.0_dp)/6.0_dp  &
    ]
    real(dp), parameter :: gauss_weights_1d_order2(2) = [0.5_dp, 0.5_dp]
    
    integer, parameter :: n_gauss_1d_order4 = 3
    real(dp), parameter :: gauss_points_1d_order4(3) = [ &
        0.5_dp - sqrt(15.0_dp)/10.0_dp, &
        0.5_dp, &
        0.5_dp + sqrt(15.0_dp)/10.0_dp  &
    ]
    real(dp), parameter :: gauss_weights_1d_order4(3) = [ &
        5.0_dp/18.0_dp, &
        8.0_dp/18.0_dp, &
        5.0_dp/18.0_dp  &
    ]
    
contains

    subroutine init(this, element_order, boundary_marker)
        class(boundary_assembly_2d_t), intent(inout) :: this
        integer, intent(in) :: element_order
        integer, intent(in), optional :: boundary_marker
        
        this%element_order = element_order
        if (present(boundary_marker)) then
            this%boundary_marker = boundary_marker
        else
            this%boundary_marker = 0  ! All boundaries
        end if
        
    end subroutine init
    
    subroutine destroy(this)
        class(boundary_assembly_2d_t), intent(inout) :: this
        
        this%element_order = 1
        this%boundary_marker = 0
        
    end subroutine destroy
    
    subroutine element_boundary_vector_p1(this, edge_vertices, g_func, be)
        class(boundary_assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: edge_vertices(2,2)  ! 2D coordinates of edge endpoints
        interface
            pure function g_func(x, y) result(g)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: g
            end function g_func
        end interface
        real(dp), intent(out) :: be(2)  ! Element boundary vector for 2 edge nodes
        
        real(dp) :: x1, y1, x2, y2, edge_length
        real(dp) :: x_phys, y_phys, g_val, s
        integer :: g, i
        type(basis_p1_2d_t) :: basis
        
        x1 = edge_vertices(1,1)
        y1 = edge_vertices(2,1)
        x2 = edge_vertices(1,2)
        y2 = edge_vertices(2,2)
        
        edge_length = sqrt((x2-x1)**2 + (y2-y1)**2)
        
        be = 0.0_dp
        
        ! Gauss quadrature along edge
        do g = 1, n_gauss_1d_order2
            s = gauss_points_1d_order2(g)
            
            ! Physical coordinates
            x_phys = (1.0_dp - s) * x1 + s * x2
            y_phys = (1.0_dp - s) * y1 + s * y2
            
            ! Neumann data
            g_val = g_func(x_phys, y_phys)
            
            ! P1 basis functions on edge: 1-s and s
            be(1) = be(1) + gauss_weights_1d_order2(g) * edge_length * g_val * (1.0_dp - s)
            be(2) = be(2) + gauss_weights_1d_order2(g) * edge_length * g_val * s
        end do
        
    end subroutine element_boundary_vector_p1
    
    subroutine element_boundary_vector_p2(this, edge_vertices, g_func, be)
        class(boundary_assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: edge_vertices(2,2)  ! 2D coordinates of edge endpoints
        interface
            pure function g_func(x, y) result(g)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: g
            end function g_func
        end interface
        real(dp), intent(out) :: be(3)  ! Element boundary vector for 3 edge nodes (P2)
        
        real(dp) :: x1, y1, x2, y2, edge_length
        real(dp) :: x_phys, y_phys, g_val, s
        real(dp) :: phi(3)
        integer :: g, i
        
        x1 = edge_vertices(1,1)
        y1 = edge_vertices(2,1)
        x2 = edge_vertices(1,2)
        y2 = edge_vertices(2,2)
        
        edge_length = sqrt((x2-x1)**2 + (y2-y1)**2)
        
        be = 0.0_dp
        
        ! Gauss quadrature along edge
        do g = 1, n_gauss_1d_order4
            s = gauss_points_1d_order4(g)
            
            ! Physical coordinates
            x_phys = (1.0_dp - s) * x1 + s * x2
            y_phys = (1.0_dp - s) * y1 + s * y2
            
            ! Neumann data
            g_val = g_func(x_phys, y_phys)
            
            ! P2 basis functions on edge
            phi(1) = (1.0_dp - s) * (1.0_dp - 2.0_dp * s)  ! Node 1
            phi(2) = s * (2.0_dp * s - 1.0_dp)             ! Node 2
            phi(3) = 4.0_dp * s * (1.0_dp - s)             ! Midpoint
            
            do i = 1, 3
                be(i) = be(i) + gauss_weights_1d_order4(g) * edge_length * g_val * phi(i)
            end do
        end do
        
    end subroutine element_boundary_vector_p2
    
    subroutine element_boundary_matrix_p1(this, edge_vertices, alpha, Ae)
        class(boundary_assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: edge_vertices(2,2)
        real(dp), intent(in) :: alpha  ! Robin BC coefficient
        real(dp), intent(out) :: Ae(2,2)  ! Element boundary matrix
        
        real(dp) :: edge_length
        real(dp) :: x1, y1, x2, y2
        
        x1 = edge_vertices(1,1)
        y1 = edge_vertices(2,1)
        x2 = edge_vertices(1,2)
        y2 = edge_vertices(2,2)
        
        edge_length = sqrt((x2-x1)**2 + (y2-y1)**2)
        
        ! Mass matrix on edge for Robin BC: alpha * ∫ u*v ds
        Ae(1,1) = alpha * edge_length / 3.0_dp
        Ae(1,2) = alpha * edge_length / 6.0_dp
        Ae(2,1) = alpha * edge_length / 6.0_dp
        Ae(2,2) = alpha * edge_length / 3.0_dp
        
    end subroutine element_boundary_matrix_p1
    
    subroutine element_boundary_matrix_p2(this, edge_vertices, alpha, Ae)
        class(boundary_assembly_2d_t), intent(in) :: this
        real(dp), intent(in) :: edge_vertices(2,2)
        real(dp), intent(in) :: alpha  ! Robin BC coefficient
        real(dp), intent(out) :: Ae(3,3)  ! Element boundary matrix for P2
        
        real(dp) :: edge_length, s
        real(dp) :: phi(3)
        integer :: i, j, g
        
        edge_length = sqrt(sum((edge_vertices(:,2) - edge_vertices(:,1))**2))
        
        Ae = 0.0_dp
        
        ! Gauss quadrature
        do g = 1, n_gauss_1d_order4
            s = gauss_points_1d_order4(g)
            
            ! P2 basis functions on edge
            phi(1) = (1.0_dp - s) * (1.0_dp - 2.0_dp * s)
            phi(2) = s * (2.0_dp * s - 1.0_dp)
            phi(3) = 4.0_dp * s * (1.0_dp - s)
            
            do i = 1, 3
                do j = 1, 3
                    Ae(i,j) = Ae(i,j) + alpha * gauss_weights_1d_order4(g) * &
                              edge_length * phi(i) * phi(j)
                end do
            end do
        end do
        
    end subroutine element_boundary_matrix_p2
    
    ! High-level assembly routine for Neumann BC
    subroutine assemble_neumann_bc(mesh, g_func, rhs, element_order, boundary_marker)
        type(mesh_2d_t), intent(in) :: mesh
        interface
            pure function g_func(x, y) result(g)
                import :: dp
                real(dp), intent(in) :: x, y
                real(dp) :: g
            end function g_func
        end interface
        real(dp), intent(inout) :: rhs(:)
        integer, intent(in), optional :: element_order
        integer, intent(in), optional :: boundary_marker
        
        type(boundary_assembly_2d_t) :: assembler
        real(dp) :: edge_vertices(2,2), be_p1(2), be_p2(3)
        integer :: e, v1, v2, order, marker
        integer :: edge_in_triangle, t, local_edge
        integer :: v3, e_mid
        
        ! Set defaults
        order = 1
        marker = 0
        if (present(element_order)) order = element_order
        if (present(boundary_marker)) marker = boundary_marker
        
        call assembler%init(order, marker)
        
        ! Loop over boundary edges
        do e = 1, mesh%n_edges
            if (.not. mesh%is_boundary_edge(e)) cycle
            
            ! Get edge vertices
            v1 = mesh%edges(1, e)
            v2 = mesh%edges(2, e)
            edge_vertices(:,1) = mesh%vertices(:, v1)
            edge_vertices(:,2) = mesh%vertices(:, v2)
            
            if (order == 1) then
                ! P1 elements
                call assembler%element_boundary_vector_p1(edge_vertices, g_func, be_p1)
                rhs(v1) = rhs(v1) + be_p1(1)
                rhs(v2) = rhs(v2) + be_p1(2)
                
            else if (order == 2) then
                ! P2 elements
                call assembler%element_boundary_vector_p2(edge_vertices, g_func, be_p2)
                
                ! Add contributions to vertices
                rhs(v1) = rhs(v1) + be_p2(1)
                rhs(v2) = rhs(v2) + be_p2(2)
                
                ! Add contribution to edge midpoint DOF
                e_mid = mesh%n_vertices + e
                rhs(e_mid) = rhs(e_mid) + be_p2(3)
            end if
        end do
        
        call assembler%destroy()
        
    end subroutine assemble_neumann_bc

end module boundary_assembly_2d_module