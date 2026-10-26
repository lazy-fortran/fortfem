program test_lab_tensor_setup_threads
    ! The threaded operator setup must reproduce the serial one: the geometry,
    ! coefficient and metric arrays bitwise (no reductions), the assembled matrix
    ! to roundoff, for every thread count, on curved 2D and 3D meshes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_set_num_threads
    use fortfem_element_to_csr, only: csr_matrix_t
    use fortfem_tensor_mesh, only: tensor_mesh_t, ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                                       assemble_tensor_operator
    implicit none
    integer, parameter :: threads(3) = [1, 4, 16]
    integer :: failure_mode = 0
    type(tensor_mesh_t) :: mesh
    type(tensor_operator_t) :: reference, op
    type(csr_matrix_t) :: coarse_reference, coarse
    integer :: dim, t, status
    do dim = 2, 3
        if (dim == 2) call ogrid_tensor_mesh(2, 3, 5, mesh, status, nphi=12, q=6)
        if (dim == 3) call ogrid_tensor_mesh(3, 2, 3, mesh, status, nphi=8, q=4)
        if (status /= 0) error stop 'mesh'
        call build_tensor_operator(mesh, coefficient, reference, status, &
                                   parallel=.false.)
        if (status /= 0) error stop 'serial operator'
        call assemble_tensor_operator(reference, coarse_reference, parallel=.false.)
        print '(a,i2,a,i7,a,i7)', 'dim', dim, ' elements', mesh%nelem, ' dofs', &
            mesh%ndof
        do t = 1, size(threads)
            call omp_set_num_threads(threads(t))
            call build_tensor_operator(mesh, coefficient, op, status)
            if (status /= 0) error stop 'threaded operator'
            if (.not. same(op%points, reference%points)) error stop 'points'
            if (.not. same(op%jac, reference%jac)) error stop 'jac'
            if (.not. same(op%metric, reference%metric)) error stop 'metric'
            if (.not. same(op%mass, reference%mass)) error stop 'mass'
            if (.not. same(op%coef, reference%coef)) error stop 'coef'
            if (.not. same(op%weight, reference%weight)) error stop 'weight'
            if (.not. same(op%volume, reference%volume)) error stop 'volume'
            if (.not. same(op%reaction, reference%reaction)) error stop 'reaction'
            if (any(op%colour_cells /= reference%colour_cells)) error stop 'colours'
            call assemble_tensor_operator(op, coarse)
            if (any(coarse%colidx /= coarse_reference%colidx)) error stop 'pattern'
            print '(a,i3,a,es10.2)', '  threads', threads(t), &
                ' matrix relative difference', &
                maxval(abs(coarse%val - coarse_reference%val))/ &
                maxval(abs(coarse_reference%val))
            if (maxval(abs(coarse%val - coarse_reference%val)) > &
                1.e-14_dp*maxval(abs(coarse_reference%val))) error stop 'matrix'
        end do
        ! A rejected coefficient (status 1) is reported for every thread count.
        do failure_mode = 1, 2
            call omp_set_num_threads(16)
            call build_tensor_operator(mesh, coefficient, op, status)
            if (status /= 1) error stop 'failure status'
            call build_tensor_operator(mesh, coefficient, op, status, &
                                       parallel=.false.)
            if (status /= 1) error stop 'serial failure status'
        end do
        failure_mode = 0
    end do
    print *, 'tensor setup threads: PASS'
contains
    logical function same(a, b)
        real(dp), intent(in) :: a(..), b(..)
        select rank (a)
        rank (2)
            select rank (b)
            rank (2)
                same = all(shape(a) == shape(b)) .and. all(a == b)
            end select
        rank (3)
            select rank (b)
            rank (3)
                same = all(shape(a) == shape(b)) .and. all(a == b)
            end select
        rank (4)
            select rank (b)
            rank (4)
                same = all(shape(a) == shape(b)) .and. all(a == b)
            end select
        end select
    end function

    subroutine coefficient(x, k, c)
        ! Thread-safe: reads only its arguments and failure_mode.
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 2.0_dp + sin(3.0_dp*x(1) + x(i))
        end do
        k(1, 2) = 0.3_dp*cos(x(1)*x(2)); k(2, 1) = k(1, 2)
        c = 1.0_dp + x(1)**2
        if (failure_mode == 1 .and. x(1) > 0.5_dp) k(1, 1) = -1
        if (failure_mode == 2 .and. x(1) > 0.5_dp) c = -1
    end subroutine
end program
