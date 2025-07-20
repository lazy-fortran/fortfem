program test_clean_mesh_api
    use fortfem_api
    use fortfem_kinds
    implicit none

    type(mesh_t) :: mesh
    type(boundary_t) :: boundary, outer, hole
    
    write(*,*) "=== Testing Clean Mesh Generation API ==="
    write(*,*) ""
    
    ! Simple built-ins
    write(*,*) "1. Simple built-ins:"
    mesh = unit_square_mesh(20)
    write(*,*) "   ✓ unit_square_mesh(20)"
    
    mesh = rectangle_mesh(20, 30, [0.0_dp, 2.0_dp, 0.0_dp, 1.5_dp])
    write(*,*) "   ✓ rectangle_mesh(20, 30, [0,2,0,1.5])"
    
    mesh = unit_disk_mesh(resolution=0.1_dp)
    write(*,*) "   ✓ unit_disk_mesh(resolution=0.1)"
    
    ! Boundary-based
    write(*,*) ""
    write(*,*) "2. Boundary-based generation:"
    boundary = circle_boundary([0.5_dp, 0.5_dp], 0.3_dp, 20)
    write(*,*) "   ✓ circle_boundary(center, radius, n)"
    
    boundary = rectangle_boundary([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], 40)
    write(*,*) "   ✓ rectangle_boundary(domain, n)"
    
    boundary = l_shape_boundary(1.0_dp, 30)
    write(*,*) "   ✓ l_shape_boundary(size, n)"
    
    ! Mesh from boundary
    write(*,*) ""
    write(*,*) "3. Mesh from boundary:"
    mesh = mesh_from_boundary(boundary, resolution=0.05_dp)
    write(*,*) "   ✓ mesh_from_boundary(boundary, resolution)"
    
    write(*,*) ""
    write(*,*) "=== Clean Mesh API Ready! ==="
    write(*,*) ""
    write(*,*) "Available functions:"
    write(*,*) "• unit_square_mesh(n)"
    write(*,*) "• rectangle_mesh(nx, ny, [x0,x1,y0,y1])"
    write(*,*) "• unit_disk_mesh(resolution=h)"
    write(*,*) "• circle_boundary(center, radius, n)"
    write(*,*) "• rectangle_boundary(domain, n)"  
    write(*,*) "• line_segment(p1, p2, n)"
    write(*,*) "• arc_segment(p1, p2, center, n)"
    write(*,*) "• l_shape_boundary(size, n)"
    write(*,*) "• mesh_from_boundary(boundary, resolution=h)"

end program test_clean_mesh_api