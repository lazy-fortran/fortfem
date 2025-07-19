program test_debug_matrix_properties
    ! Debug matrix properties for the failing case
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    integer :: n, nx, ny, n_dofs, i, j, n_boundary
    real(dp), allocatable :: A(:,:), b(:)
    
    ! Test with the failing mesh size: n=4 (level=1)
    n = 4
    nx = n + 1
    ny = n + 1
    
    call mesh%create_rectangular(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()
    
    call space%init(mesh)
    n_dofs = space%get_n_dofs()
    
    allocate(A(n_dofs, n_dofs), b(n_dofs))
    A = 0.0_dp
    b = 0.0_dp
    
    print *, "Testing matrix assembly for n=", n, "DOFs=", n_dofs
    
    ! Assemble curl-curl and mass terms
    call assemble_curl_curl_bilinear_form(mesh, space, A)
    call assemble_mass_bilinear_form(mesh, space, A)
    call assemble_rhs_linear_form(mesh, space, b)
    
    print *, "Before boundary conditions:"
    print *, "  Matrix diagonal min/max:", minval([(A(i,i), i=1,n_dofs)]), maxval([(A(i,i), i=1,n_dofs)])
    print *, "  Matrix non-zeros:", count(abs(A) > 1e-14)
    print *, "  RHS non-zeros:", count(abs(b) > 1e-14)
    
    ! Apply boundary conditions
    call apply_essential_boundary_conditions(mesh, space, A, b)
    
    print *, "After boundary conditions:"
    print *, "  Matrix diagonal min/max:", minval([(A(i,i), i=1,n_dofs)]), maxval([(A(i,i), i=1,n_dofs)])
    print *, "  Matrix non-zeros:", count(abs(A) > 1e-14)
    print *, "  RHS non-zeros:", count(abs(b) > 1e-14)
    
    ! Count boundary DOFs
    n_boundary = 0
    do i = 1, n_dofs
        if (mesh%is_boundary_edge(i)) n_boundary = n_boundary + 1
    end do
    print *, "  Boundary DOFs:", n_boundary, "Interior DOFs:", n_dofs - n_boundary
    
    call space%destroy()
    call mesh%destroy()
    deallocate(A, b)

contains

    ! Include the same assembly routines from the curl-curl solver
    subroutine assemble_curl_curl_bilinear_form(mesh, space, A)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: A(:,:)
        
        integer :: t, i, j, dof_i, dof_j, triangle_dofs(3)
        real(dp) :: curl_i(3), curl_j(3), triangle_area
        real(dp) :: xi, eta, contribution
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curl_i)
            curl_j = curl_i
            call space%get_triangle_dofs(t, triangle_dofs)
            
            do i = 1, 3
                dof_i = triangle_dofs(i) + 1
                do j = 1, 3
                    dof_j = triangle_dofs(j) + 1
                    contribution = curl_i(i) * curl_j(j) * triangle_area
                    A(dof_i, dof_j) = A(dof_i, dof_j) + contribution
                end do
            end do
        end do
    end subroutine
    
    subroutine assemble_mass_bilinear_form(mesh, space, A)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: A(:,:)
        
        integer :: t, i, j, dof_i, dof_j, triangle_dofs(3)
        real(dp) :: values_i(2,3), values_j(2,3), triangle_area
        real(dp) :: xi, eta, contribution
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            call evaluate_edge_basis_2d(xi, eta, triangle_area, values_i)
            values_j = values_i
            call space%get_triangle_dofs(t, triangle_dofs)
            
            do i = 1, 3
                dof_i = triangle_dofs(i) + 1
                do j = 1, 3
                    dof_j = triangle_dofs(j) + 1
                    contribution = (values_i(1,i) * values_j(1,j) + &
                                   values_i(2,i) * values_j(2,j)) * triangle_area
                    A(dof_i, dof_j) = A(dof_i, dof_j) + contribution
                end do
            end do
        end do
    end subroutine
    
    subroutine assemble_rhs_linear_form(mesh, space, b)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: b(:)
        
        integer :: t, i, dof_i, triangle_dofs(3)
        real(dp) :: values_i(2,3), triangle_area
        real(dp) :: xi, eta, x_phys, y_phys, contribution
        real(dp) :: f_exact(2)
        
        b = 0.0_dp
        
        do t = 1, mesh%n_triangles
            triangle_area = compute_triangle_area(mesh, t)
            xi = 1.0_dp/3.0_dp
            eta = 1.0_dp/3.0_dp
            
            call map_to_physical(mesh, t, xi, eta, x_phys, y_phys)
            f_exact(1) = x_phys * y_phys
            f_exact(2) = x_phys**2
            
            call evaluate_edge_basis_2d(xi, eta, triangle_area, values_i)
            call space%get_triangle_dofs(t, triangle_dofs)
            
            do i = 1, 3
                dof_i = triangle_dofs(i) + 1
                contribution = (f_exact(1) * values_i(1,i) + &
                               f_exact(2) * values_i(2,i)) * triangle_area
                b(dof_i) = b(dof_i) + contribution
            end do
        end do
    end subroutine
    
    subroutine apply_essential_boundary_conditions(mesh, space, A, b)
        type(mesh_2d_t), intent(in) :: mesh
        type(hcurl_space_t), intent(in) :: space
        real(dp), intent(inout) :: A(:,:), b(:)
        
        integer :: edge_idx, n_dofs
        
        n_dofs = space%get_n_dofs()
        
        do edge_idx = 1, n_dofs
            if (mesh%is_boundary_edge(edge_idx)) then
                A(edge_idx, :) = 0.0_dp
                A(:, edge_idx) = 0.0_dp
                A(edge_idx, edge_idx) = 1.0_dp
                b(edge_idx) = 0.0_dp
            end if
        end do
    end subroutine
    
    function compute_triangle_area(mesh, triangle_idx) result(area)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp) :: area
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        area = 0.5_dp * abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3))
    end function
    
    subroutine map_to_physical(mesh, triangle_idx, xi, eta, x_phys, y_phys)
        type(mesh_2d_t), intent(in) :: mesh
        integer, intent(in) :: triangle_idx
        real(dp), intent(in) :: xi, eta
        real(dp), intent(out) :: x_phys, y_phys
        
        real(dp) :: x1, y1, x2, y2, x3, y3
        
        x1 = mesh%vertices(1, mesh%triangles(1, triangle_idx))
        y1 = mesh%vertices(2, mesh%triangles(1, triangle_idx))
        x2 = mesh%vertices(1, mesh%triangles(2, triangle_idx))
        y2 = mesh%vertices(2, mesh%triangles(2, triangle_idx))
        x3 = mesh%vertices(1, mesh%triangles(3, triangle_idx))
        y3 = mesh%vertices(2, mesh%triangles(3, triangle_idx))
        
        x_phys = x1 + (x2 - x1) * xi + (x3 - x1) * eta
        y_phys = y1 + (y2 - y1) * xi + (y3 - y1) * eta
    end subroutine

end program test_debug_matrix_properties