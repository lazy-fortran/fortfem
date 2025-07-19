program test_debug_basis_verification
    ! Verify our basis functions match FreeFEM's RT0Ortho implementation
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d
    use fortfem_hcurl_space
    use fortfem_basis_edge_2d
    implicit none
    
    type(mesh_2d_t) :: mesh
    type(hcurl_space_t) :: space
    integer :: n, nx, ny
    real(dp) :: xi, eta, triangle_area
    real(dp) :: values(2,3), curls(3)
    integer :: t
    
    ! Create simple 2x2 mesh (same as FreeFEM square(4,4))
    n = 4
    nx = n + 1  ! 5x5 grid
    ny = n + 1
    
    call mesh%create_rectangular(nx, ny, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp)
    call mesh%build_edge_connectivity()
    call mesh%build_edge_dof_numbering()
    
    call space%init(mesh)
    
    print *, "Basis Function Verification for RT0Ortho Elements"
    print *, "================================================"
    print *, "Mesh: 5x5 vertices, ", mesh%n_triangles, "triangles"
    print *, ""
    
    ! Test first triangle
    t = 1
    triangle_area = compute_triangle_area(mesh, t)
    
    print *, "Triangle 1 area:", triangle_area
    print *, "Triangle vertices:"
    print *, "  V1:", mesh%vertices(:, mesh%triangles(1, t))
    print *, "  V2:", mesh%vertices(:, mesh%triangles(2, t))  
    print *, "  V3:", mesh%vertices(:, mesh%triangles(3, t))
    print *, ""
    
    ! Test basis functions at different points
    print *, "Basis function values at triangle center (1/3, 1/3):"
    xi = 1.0_dp/3.0_dp
    eta = 1.0_dp/3.0_dp
    
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
    
    print *, "  Edge 1: value = [", values(1,1), ",", values(2,1), "], curl =", curls(1)
    print *, "  Edge 2: value = [", values(1,2), ",", values(2,2), "], curl =", curls(2)
    print *, "  Edge 3: value = [", values(1,3), ",", values(2,3), "], curl =", curls(3)
    print *, ""
    
    ! Test at different points
    print *, "Basis function values at (0.5, 0.2):"
    xi = 0.5_dp
    eta = 0.2_dp
    
    call evaluate_edge_basis_2d(xi, eta, triangle_area, values)
    call evaluate_edge_basis_curl_2d(xi, eta, triangle_area, curls)
    
    print *, "  Edge 1: value = [", values(1,1), ",", values(2,1), "], curl =", curls(1)
    print *, "  Edge 2: value = [", values(1,2), ",", values(2,2), "], curl =", curls(2)
    print *, "  Edge 3: value = [", values(1,3), ",", values(2,3), "], curl =", curls(3)
    print *, ""
    
    ! Verify curl computation manually
    print *, "Expected RT0Ortho curl values:"
    print *, "  Based on our implementation:"
    print *, "    curl(φ₁) = 0.0"
    print *, "    curl(φ₂) = 1/area =", 1.0_dp/triangle_area
    print *, "    curl(φ₃) = -1/area =", -1.0_dp/triangle_area
    print *, ""
    
    ! Test constant curl field representation
    print *, "Testing if constant curl = 2 can be represented:"
    print *, "  For E = [-y, x] with curl(E) = 2"
    print *, "  We need: c₁*0 + c₂*(1/A) + c₃*(-1/A) = 2"
    print *, "  This gives: c₂ - c₃ = 2*A =", 2.0_dp*triangle_area
    print *, "  With A =", triangle_area, ", we need c₂ - c₃ =", 2.0_dp*triangle_area
    
    call space%destroy()
    call mesh%destroy()

contains

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
    end function compute_triangle_area

end program test_debug_basis_verification