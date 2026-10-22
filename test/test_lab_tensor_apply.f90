program test_lab_tensor_apply
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
             apply_tensor_operator, assemble_tensor_operator, tensor_diagonal, metric_at
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: op, stream
    type(csr_matrix_t) :: a
    real(dp), allocatable :: x(:), y(:), z(:), v(:), diagonal(:)
    real(dp) :: err, denom, jac(3, 3), kmat(3, 3), g(3, 3), mass, det, zphys(3), vref(3)
    integer :: dim, p, status, i, j
    jac = reshape([2.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, 3.0_dp, 0.0_dp, &
                   0.2_dp, 0.4_dp, 4.0_dp], [3, 3])
    kmat = reshape([2.0_dp, 0.1_dp, 0.2_dp, 0.1_dp, 3.0_dp, 0.3_dp, &
                    0.2_dp, 0.3_dp, 4.0_dp], [3, 3])
    zphys = [0.7_dp, -1.3_dp, 0.9_dp]
    vref = matmul(transpose(jac), zphys)
    call metric_at(jac, kmat, 0.6_dp, 0.8_dp, g, mass, det)
    err = abs(dot_product(vref, matmul(g, vref)) - &
              0.6_dp*24*dot_product(zphys, matmul(kmat, zphys)))
    if (err > 1.e-12_dp .or. abs(mass - 0.6_dp*24*0.8_dp) > 1.e-14_dp) &
        error stop 'independent physical-energy pullback'
    jac(:, 1) = 0
    call metric_at(jac, kmat, 0.6_dp, 0.8_dp, g, mass, det)
    if (abs(det) > tiny(1.0_dp)) error stop 'singular Jacobian rejected before division'
    do dim = 2, 3
        do p = 1, 6
            call structured_tensor_mesh(dim, [2, 2, 2], p, mesh, status, &
                               map=curved, periodic=[.true., .false., .false.], q=p + 1)
            if (status /= 0) error stop 'mesh'
            call build_tensor_operator(mesh, coefficient, op, status)
            if (status /= 0) error stop 'operator'
            call build_tensor_operator(mesh, coefficient, stream, status, &
                                       precompute=.false.)
            if (status /= 0) error stop 'stream operator'
            allocate (x(mesh%ndof), y(mesh%ndof), z(mesh%ndof), v(mesh%ndof))
            do i = 1, mesh%ndof
                x(i) = sin(real(13*i, dp))
            end do
            call apply_tensor_operator(op, x, y)
            call apply_tensor_operator(stream, x, v)
            denom = max(1.0_dp, norm2(y))
            if (norm2(y - v)/denom > 1.e-13_dp) error stop 'recompute metric'
            call assemble_tensor_operator(op, a)
            z = 0
            do i = 1, a%n
                do j = a%rowptr(i), a%rowptr(i + 1) - 1
                    z(i) = z(i) + a%val(j)*x(a%colidx(j))
                end do
            end do
            err = norm2(y - z)/denom
            print '(a,i1,a,i1,a,es10.2)', 'dim=', dim, ' p=', p, ' relative=', err
            if (err > 1.e-12_dp) error stop 'MF versus independent CSR'
            call tensor_diagonal(op, diagonal)
            do i = 1, a%n
                do j = a%rowptr(i), a%rowptr(i + 1) - 1
                    if (a%colidx(j) /= i) cycle
                    if (abs(diagonal(i) - a%val(j)) > &
                        1.e-12_dp*max(1.0_dp, abs(a%val(j)))) &
                        error stop 'exact diagonal'
                end do
            end do
            deallocate (x, y, z, v, diagonal)
        end do
    end do
contains
    subroutine curved(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x = xi
        x(2) = xi(2) + 0.06_dp*sin(2*acos(-1.0_dp)*xi(1))*xi(2)
        if (size(x) == 3) x(3) = xi(3) + 0.04_dp*xi(1)*xi(2)*(1 - xi(3))
    end subroutine
    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1 + 0.1_dp*x(i)**2
        end do
        k(1, 2) = 0.17_dp
        k(2, 1) = 0.17_dp
        c = 0.3_dp
    end subroutine
end program
