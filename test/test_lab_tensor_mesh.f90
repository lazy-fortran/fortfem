program test_lab_tensor_mesh
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh, &
        build_tensor_mesh, ogrid_tensor_mesh
    implicit none
    type(tensor_mesh_t) :: mesh, oriented
    integer :: status, d, p, e, i, j, n
    integer, allocatable :: corners(:,:), gather(:,:)
    real(dp), allocatable :: geom(:,:,:)
    do d = 2, 3
        do p = 1, 6
            call structured_tensor_mesh(d, [2,2,2], p, mesh, status, &
                map=curved, periodic=[.true.,.false.,.false.])
            if (status /= 0) error stop 'structured mesh'
            if (mesh%ndof /= (2*p)*(2*p+1)**(d-1)) error stop 'periodic dofs'
            ! Reverse two reference directions: same physical cells and orientation.
            n = p+1
            geom = mesh%geometry
            corners = mesh%vertices
            do e = 1, mesh%nelem
                do i = 1, size(geom,2)
                    j = (i-1)/n**2*n**2 + (n-1-modulo((i-1)/n,n))*n + &
                        n-modulo(i-1,n)
                    geom(:,i,e) = mesh%geometry(:,j,e)
                end do
                do i = 1, 2**d
                    j = (i-1)/4*4 + (1-modulo((i-1)/2,2))*2 + &
                        2-modulo(i-1,2)
                    corners(i,e) = mesh%vertices(j,e)
                end do
            end do
            ! Explicit gather maps also support periodic multigraph topology.
            gather = mesh%dof
            do e = 1, mesh%nelem
                do i = 1, size(gather,1)
                    j = (i-1)/n**2*n**2 + (n-1-modulo((i-1)/n,n))*n + &
                        n-modulo(i-1,n)
                    gather(i,e) = mesh%dof(j,e)
                end do
            end do
            call build_tensor_mesh(d,p,p,corners,geom,oriented,status,gather)
            if (status /= 0) error stop 'oriented mesh'
            if (maxval(abs(mesh%coordinates-oriented%coordinates)) > &
                1.e-13_dp) error stop 'oriented gather coordinates'
            deallocate(geom,corners,gather)
        end do
        call ogrid_tensor_mesh(d, 2, 3, mesh, status, nphi=5)
        if (status /= 0) error stop 'O-grid mesh'
        if (minval(sum(mesh%coordinates(1:2,:)**2,dim=1)) < 0.0_dp) &
            error stop 'finite coordinates'
    end do
    print *, 'PASS: periodic, oriented and axis-crossing tensor meshes'
contains
    subroutine curved(xi,x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x = xi
        x(2) = xi(2) + 0.04_dp*sin(2*acos(-1.0_dp)*xi(1))*xi(2)
    end subroutine
end program
