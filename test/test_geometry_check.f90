program test_geometry_check
    ! Check geometric properties of the 3 input points
    use fortfem_kinds
    use geometric_predicates
    use delaunay_types
    implicit none

    type(point_t) :: p1, p2, p3
    integer :: orient
    real(dp) :: area
    
    write(*,*) "=== Geometry Check ==="
    
    ! Create the 3 points
    p1%x = 0.0_dp; p1%y = 0.0_dp
    p2%x = 1.0_dp; p2%y = 0.0_dp  
    p3%x = 0.5_dp; p3%y = 1.0_dp
    
    write(*,'(A,2F8.3)') "Point 1: ", p1%x, p1%y
    write(*,'(A,2F8.3)') "Point 2: ", p2%x, p2%y
    write(*,'(A,2F8.3)') "Point 3: ", p3%x, p3%y
    
    ! Check orientation
    orient = orientation(p1, p2, p3)
    write(*,*) "Orientation P1->P2->P3:", orient
    if (orient == ORIENTATION_CCW) then
        write(*,*) "  Counter-clockwise (valid triangle)"
    else if (orient == ORIENTATION_CW) then
        write(*,*) "  Clockwise (valid triangle, reversed)"
    else
        write(*,*) "  Collinear (degenerate)"
    end if
    
    ! Calculate area
    area = triangle_area(p1, p2, p3)
    write(*,'(A,F12.6)') "Triangle area: ", area
    
    if (area > 1.0e-10_dp) then
        write(*,*) "Triangle is valid (non-degenerate)"
    else
        write(*,*) "Triangle is degenerate (zero area)"
    end if

end program test_geometry_check