program test_lab_tensor_mesh
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh, &
                                   build_tensor_mesh, ogrid_tensor_mesh
    implicit none
    type(tensor_mesh_t) :: mesh, oriented
    integer :: status, d, p, e, i, j, n, expected
    real(dp) :: radius
    integer, allocatable :: corners(:, :), gather(:, :)
    real(dp), allocatable :: geom(:, :, :)
    do d = 2, 3
        do p = 1, 6
            call structured_tensor_mesh(d, [2, 2, 2], p, mesh, status, &
                                        map=curved, periodic=[.true., .false., .false.])
            if (status /= 0) error stop 'structured mesh'
            if (mesh%ndof /= (2*p)*(2*p + 1)**(d - 1)) error stop 'periodic dofs'
            ! Reverse two reference directions: same physical cells and orientation.
            n = p + 1
            geom = mesh%geometry
            corners = mesh%vertices
            do e = 1, mesh%nelem
                do i = 1, size(geom, 2)
                    j = (i - 1)/n**2*n**2 + (n - 1 - modulo((i - 1)/n, n))*n + &
                        n - modulo(i - 1, n)
                    geom(:, i, e) = mesh%geometry(:, j, e)
                end do
                do i = 1, 2**d
                    j = (i - 1)/4*4 + (1 - modulo((i - 1)/2, 2))*2 + &
                        2 - modulo(i - 1, 2)
                    corners(i, e) = mesh%vertices(j, e)
                end do
            end do
            ! Explicit gather maps also support periodic multigraph topology.
            gather = mesh%dof
            do e = 1, mesh%nelem
                do i = 1, size(gather, 1)
                    j = (i - 1)/n**2*n**2 + (n - 1 - modulo((i - 1)/n, n))*n + &
                        n - modulo(i - 1, n)
                    gather(i, e) = mesh%dof(j, e)
                end do
            end do
            call build_tensor_mesh(d, p, p, corners, geom, oriented, status, gather)
            if (status /= 0) error stop 'oriented mesh'
            if (maxval(abs(mesh%coordinates - oriented%coordinates)) > &
                1.e-13_dp) error stop 'oriented gather coordinates'
            deallocate (geom, corners, gather)
            call automatic_orientation_case(d, p)
        end do
        call ogrid_tensor_mesh(d, 2, 3, mesh, status, nphi=5)
        if (status /= 0) error stop 'O-grid mesh'
        if (.not. all(ieee_is_finite(mesh%coordinates))) error stop 'finite coordinates'
        expected = 5*(2*3)**2 + 2*(2*3) + 1
        if (d == 3) expected = expected*5*3
        if (mesh%ndof /= expected) error stop 'O-grid seam identifications'
        do i = 1, mesh%ndof
            if (d == 2) then
                radius = norm2(mesh%coordinates(:, i))
            else
           radius = norm2([norm2(mesh%coordinates(1:2, i)) - 3, mesh%coordinates(3, i)])
            end if
            if (radius < 1.e-12_dp .and. mesh%boundary(i)) &
                error stop 'axis is interior'
        end do
    end do
    print *, 'PASS: periodic, oriented and axis-crossing tensor meshes'
contains
    subroutine automatic_orientation_case(dim, degree)
        integer, intent(in) :: dim, degree
        type(tensor_mesh_t) :: original, reoriented
        real(dp), allocatable :: geometry(:, :, :)
        integer, allocatable :: vertices(:, :)
        integer :: cell, a, b, n, status
       call structured_tensor_mesh(dim, [2, 2, 2], degree, original, status, map=curved)
        if (status /= 0) error stop 'automatic orientation source'
        geometry = original%geometry
        vertices = original%vertices
        n = degree + 1
        ! Adjacent cells have different edge and face reference orientations.
        do cell = 1, original%nelem, 2
            do a = 1, size(geometry, 2)
         b = (a - 1)/n**2*n**2 + (n - 1 - modulo((a - 1)/n, n))*n + n - modulo(a - 1, n)
                geometry(:, a, cell) = original%geometry(:, b, cell)
            end do
            do a = 1, 2**dim
                b = (a - 1)/4*4 + (1 - modulo((a - 1)/2, 2))*2 + 2 - modulo(a - 1, 2)
                vertices(a, cell) = original%vertices(b, cell)
            end do
        end do
        call build_tensor_mesh(dim, degree, degree, vertices, geometry, &
                               reoriented, status)
        if (status /= 0) error stop 'automatic oriented gather'
        if (reoriented%ndof /= original%ndof) error stop 'shared oriented dof count'
        do cell = 1, original%nelem
            do a = 1, size(geometry, 2)
                if (maxval(abs(reoriented%coordinates(:, reoriented%dof(a, cell)) - &
                 geometry(:, a, cell))) > 1.e-13_dp) error stop 'oriented face location'
            end do
        end do
    end subroutine
    subroutine curved(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x = xi
        x(2) = xi(2) + 0.04_dp*sin(2*acos(-1.0_dp)*xi(1))*xi(2)
    end subroutine
end program
