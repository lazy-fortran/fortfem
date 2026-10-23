module lab_tensor_seam_cases
    !! Shared cases of test_lab_tensor_seam and test_lab_tensor_seam_slow.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t, structured_tensor_mesh, &
                                   ogrid_tensor_mesh
    use fortfem_tensor_operator, only: tensor_operator_t, build_tensor_operator, &
                    apply_tensor_operator, tensor_load, tensor_l2_error, scalar_function
    use fortfem_tensor_seam, only: tensor_seam_t, close_tensor_seam, &
                                   apply_seam_operator, assemble_seam_operator, &
                                   seam_expand, SEAM_MORTAR, SEAM_INTERPOLATION
    use fortfem_tensor_solver, only: tensor_solve_info_t, solve_tensor_diffusion
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    private
    public :: run_tensor_seam
    real(dp), parameter :: pi = acos(-1.0_dp)
    !! Seam shift s(x) = shift0 + shift2*x**2 in periods of the slab y.
    real(dp) :: shift0 = 0.37_dp, shift2 = 0.21_dp
    !! Torus rotation iota(s) = twist0 + twist2*s turns per toroidal turn.
    real(dp) :: twist0 = 0.37_dp, twist2 = 0.21_dp
    real(dp), parameter :: translate(3, 4) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1, &
                                                      0, 0, -1]*1.0_dp, [3, 4])
    real(dp), parameter :: identity(3, 4) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1, &
                                                     0, 0, 0]*1.0_dp, [3, 4])
    real(dp), parameter :: lattice(3, 3) = reshape([0, 0, 0, 0, -1, 0, 0, 1, 0]* &
                                                   1.0_dp, [3, 3])
contains
    subroutine run_tensor_seam(slow)
        !! Standard: constraint algebra and bounded solves; slow: rate studies
        !! and the costlier twisted-torus degrees.
        logical, intent(in) :: slow
        if (slow) then
            call check_operator_agreement([1, 0], [3, 4])
            call check_reproduction([1, 0], [3, 4])
            call check_rigid_twist(2, 3)
            call check_slab_rates()
            call check_torus_rates()
        else
            call check_zero_shift()
            call check_operator_agreement([1, 6], [1, 2])
            call check_reproduction([1, 4], [1, 2])
            call check_bounded_solves()
            call check_rigid_twist(1, 2)
        end if
        print *, 'PASS: periodic seam constraints (criterion 3)'
    end subroutine

    subroutine sheared_slab(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        x(1) = xi(1)
        x(2) = xi(2) + (shift0 + shift2*xi(1)**2)*xi(3)
        x(3) = xi(3)
    end subroutine

    subroutine twisted_torus(xi, x)
        !! O-grid disk turned by 2 pi iota(s) phi/(2 pi), R = 3 + u, Z = v,
        !! with y = -R sin(phi) so that the cells are positively oriented.
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: s, angle, u, v, phi
        s = xi(1)**2 + xi(2)**2
        angle = 2*pi*xi(3)*(twist0 + twist2*s)
        u = cos(angle)*xi(1) - sin(angle)*xi(2)
        v = sin(angle)*xi(1) + cos(angle)*xi(2)
        phi = 2*pi*xi(3)
        x(1) = (3 + u)*cos(phi)
        x(2) = -(3 + u)*sin(phi)
        x(3) = v
    end subroutine

    subroutine isotropic(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
        c = 0
        if (.false.) k(1, 1) = x(1)
    end subroutine

    function slab_exact(x) result(u)
        real(dp), intent(in) :: x(:)
        real(dp) :: u
        u = cos(x(1))*(2 + sin(2*pi*x(2) + 4*pi*x(3)))
    end function

    function slab_forcing(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = slab_exact(x) + 20*pi**2*cos(x(1))*sin(2*pi*x(2) + 4*pi*x(3))
    end function

    function torus_exact(x) result(u)
        real(dp), intent(in) :: x(:)
        real(dp) :: u
        u = sin(x(1)/2)*cos(x(2)/3)*exp(x(3)/2)
    end function

    function torus_forcing(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = torus_exact(x)/9
    end function

    subroutine slab(cells, p, mesh, master, slave)
        integer, intent(in) :: cells(3), p
        type(tensor_mesh_t), intent(out) :: mesh
        integer, allocatable, intent(out) :: master(:, :), slave(:, :)
        integer :: status
        call structured_tensor_mesh(3, cells, p, mesh, status, map=sheared_slab, &
                                    periodic=[.false., .true., .false.], q=p + 1, &
                                    dirichlet=[.true., .false., .false.])
        if (status /= 0) error stop 'slab mesh'
        call layer_faces(cells(1)*cells(2), cells(3), master, slave)
    end subroutine

    subroutine torus(n, nphi, p, mesh, master, slave)
        integer, intent(in) :: n, nphi, p
        type(tensor_mesh_t), intent(out) :: mesh
        integer, allocatable, intent(out) :: master(:, :), slave(:, :)
        integer :: status
        call ogrid_tensor_mesh(3, n, p, mesh, status, nphi=nphi, q=p + 2, &
                               map=twisted_torus, periodic=.false.)
        if (status /= 0) error stop 'torus mesh'
        call layer_faces(5*n*n, nphi, master, slave)
    end subroutine

    subroutine layer_faces(per_layer, layers, master, slave)
        !! Third-direction faces of the first and last cell layers.
        integer, intent(in) :: per_layer, layers
        integer, allocatable, intent(out) :: master(:, :), slave(:, :)
        integer :: e
        allocate (master(3, per_layer), slave(3, per_layer))
        do e = 1, per_layer
            master(:, e) = [e, 3, 0]
            slave(:, e) = [e + per_layer*(layers - 1), 3, 1]
        end do
    end subroutine

    subroutine check_zero_shift()
        type(tensor_mesh_t) :: open, closed, periodic
        type(tensor_seam_t) :: seam
        type(tensor_operator_t) :: a, b
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp), allocatable :: x(:), ya(:), yb(:)
        integer :: p, status, i
        shift0 = 0
        shift2 = 0
        do p = 1, 4
            call slab([2, 3, 4], p, open, master, slave)
            call close_tensor_seam(open, master, slave, closed, seam, status, &
                                   motion=translate, images=lattice)
            if (status /= 0) error stop 'zero-shift seam'
            if (.not. seam%conforming .or. seam%nslave /= 0) &
                error stop 'zero shift is not an identification'
            call structured_tensor_mesh(3, [2, 3, 4], p, periodic, status, &
                          map=sheared_slab, periodic=[.false., .true., .true.], q=p + 1)
            if (status /= 0) error stop 'periodic mesh'
            if (closed%ndof /= periodic%ndof) error stop 'zero-shift ndof'
            if (any(closed%dof /= periodic%dof)) error stop 'zero-shift gather map'
         if (any(closed%vertices /= periodic%vertices)) error stop 'zero-shift vertices'
            if (any(closed%boundary .neqv. periodic%boundary)) &
                error stop 'zero-shift boundary'
            call build_tensor_operator(closed, isotropic, a, status)
            call build_tensor_operator(periodic, isotropic, b, status)
            allocate (x(closed%ndof), ya(closed%ndof), yb(closed%ndof))
            do i = 1, closed%ndof
                x(i) = sin(real(17*i, dp))
            end do
            call apply_seam_operator(a, seam, x, ya)
            call apply_tensor_operator(b, x, yb)
            ! Bit-for-bit, not merely within a tolerance.
            if (any(transfer(ya, 1_8, size(ya)) /= transfer(yb, 1_8, size(yb)))) &
                error stop 'zero shift is not bit-for-bit periodic'
            deallocate (x, ya, yb)
        end do
        shift0 = 0.37_dp
        shift2 = 0.21_dp
        print *, 'zero shift: identification equals periodic map bit-for-bit, p=1..4'
    end subroutine

    subroutine check_operator_agreement(pslab, ptorus)
        !! Degree ranges per geometry; an empty range skips the geometry.
        integer, intent(in) :: pslab(2), ptorus(2)
        type(tensor_mesh_t) :: mesh, closed
        type(tensor_seam_t) :: seam
        type(tensor_operator_t) :: op
        type(csr_matrix_t) :: matrix
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp), allocatable :: x(:), ymf(:), ycsr(:)
        real(dp) :: diff
        integer :: p, method, status, i, k, geometry
        print *, 'geometry method p slaves reduced |MF-CSR|/|CSR| (criterion 3)'
        do geometry = 1, 2
            do method = SEAM_MORTAR, SEAM_INTERPOLATION
                do p = merge(pslab(1), ptorus(1), geometry == 1), &
                    merge(pslab(2), ptorus(2), geometry == 1)
                    if (geometry == 1) then
                        call slab([2, 3, 3], p, mesh, master, slave)
                        call close_tensor_seam(mesh, master, slave, closed, seam, &
                                status, motion=translate, images=lattice, method=method)
                    else
                        call torus(1, 3, p, mesh, master, slave)
                        call close_tensor_seam(mesh, master, slave, closed, seam, &
                                               status, motion=identity, method=method)
                    end if
                    if (status /= 0) print *, 'status', status, geometry, method, p
                    if (status /= 0) error stop 'seam constraint'
                    if (seam%conforming) error stop 'twisted seam matched'
                    call build_tensor_operator(closed, isotropic, op, status)
                    if (status /= 0) error stop 'seam operator'
                    call assemble_seam_operator(op, seam, matrix)
                    allocate (x(seam%nreduced), ymf(seam%nreduced), ycsr(seam%nreduced))
                    do i = 1, seam%nreduced
                        x(i) = cos(real(29*i, dp))
                    end do
                    call apply_seam_operator(op, seam, x, ymf)
                    do i = 1, matrix%n
                        ycsr(i) = 0
                        do k = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                            ycsr(i) = ycsr(i) + matrix%val(k)*x(matrix%colidx(k))
                        end do
                    end do
                    diff = norm2(ymf - ycsr)/norm2(ycsr)
                    print '(3i7,2i8,es12.3)', geometry, method, p, seam%nslave, &
                        seam%nreduced, diff
                    if (diff > 1.e-12_dp) error stop 'MF and assembled seam operators'
                    deallocate (x, ymf, ycsr)
                end do
            end do
        end do
    end subroutine

    subroutine check_reproduction(pslab, ptorus)
        integer, intent(in) :: pslab(2), ptorus(2)
        type(tensor_mesh_t) :: mesh, closed
        type(tensor_seam_t) :: seam
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp), allocatable :: reduced(:), full(:)
        real(dp) :: err(3)
        integer :: p, method, status, f, k, geometry
        print *, 'geometry method p const linear-x quadratic-x  max slave error'
        do geometry = 1, 2
            do method = SEAM_MORTAR, SEAM_INTERPOLATION
                do p = merge(pslab(1), ptorus(1), geometry == 1), &
                    merge(pslab(2), ptorus(2), geometry == 1)
                    if (geometry == 1) then
                        call slab([2, 3, 2], p, mesh, master, slave)
                        call close_tensor_seam(mesh, master, slave, closed, seam, &
                                status, motion=translate, images=lattice, method=method)
                    else
                        call torus(1, 2, p, mesh, master, slave)
                        call close_tensor_seam(mesh, master, slave, closed, seam, &
                                               status, motion=identity, method=method)
                    end if
                    if (status /= 0) error stop 'reproduction seam'
                    allocate (reduced(seam%nreduced), full(seam%nfull))
                    err = 0
                    do f = 1, merge(3, 1, geometry == 1)
                        if (f == 3 .and. p < 2) cycle
                        do k = 1, seam%nreduced
                            reduced(k) = probe(f, mesh%coordinates(:, seam%primary(k)))
                        end do
                        call seam_expand(seam, reduced, full)
                        do k = 1, seam%nslave
                            err(f) = max(err(f), abs(full(seam%slave(k)) - &
                                          probe(f, mesh%coordinates(:, seam%slave(k)))))
                        end do
                    end do
                    print '(3i7,3es12.3)', geometry, method, p, err
                    if (maxval(err) > 1.e-11_dp) error stop 'seam reproduction'
                    deallocate (reduced, full)
                end do
            end do
        end do
    end subroutine

    real(dp) function probe(f, x)
        !! Constants everywhere; x and x^2 are in both slab trace spaces.
        integer, intent(in) :: f
        real(dp), intent(in) :: x(:)
        select case (f)
        case (1)
            probe = 1
        case (2)
            probe = 0.3_dp + 2*x(1)
        case default
            probe = x(1)**2 - 0.5_dp*x(1)
        end select
    end function

    subroutine solve_error(mesh, seam, exact, forcing, p, error, iterations)
        type(tensor_mesh_t), intent(in) :: mesh
        type(tensor_seam_t), intent(in) :: seam
        procedure(scalar_function) :: exact, forcing
        integer, intent(in) :: p
        real(dp), intent(out) :: error
        integer, intent(out) :: iterations
        type(tensor_operator_t) :: op
        type(tensor_solve_info_t) :: info
        real(dp), allocatable :: rhs(:), x(:), bc(:)
        integer :: status, i
        call build_tensor_operator(mesh, isotropic, op, status, nq=p + 3)
        if (status /= 0) error stop 'MMS operator'
        call tensor_load(op, forcing, rhs)
        allocate (x(mesh%ndof), bc(mesh%ndof))
        do i = 1, mesh%ndof
            bc(i) = exact(mesh%coordinates(:, i))
        end do
        x = 0
        call solve_tensor_diffusion(mesh, isotropic, rhs, 1.e-11_dp, x, info, bc=bc, &
                                    maxiter=200, seam=seam)
        if (info%status /= 0) error stop 'MMS seam PCG'
        iterations = info%iterations
        error = tensor_l2_error(op, x, exact)
    end subroutine

    subroutine check_slab_rates()
        type(tensor_mesh_t) :: mesh, closed
        type(tensor_seam_t) :: seam
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp) :: error(2), previous(2), rate(2)
        integer :: p, n, method, status, its(2)
        print *, 'sheared slab: p n dofs its(mortar,interp) L2(mortar,interp) rates'
        do p = 2, 3
            previous = 0
            do n = 2, 8, 2
                if (n == 6) cycle
                do method = SEAM_MORTAR, SEAM_INTERPOLATION
                    call slab([n, n, n], p, mesh, master, slave)
                    call close_tensor_seam(mesh, master, slave, closed, seam, status, &
                                        motion=translate, images=lattice, method=method)
                    if (status /= 0) error stop 'slab seam'
                    call solve_error(closed, seam, slab_exact, slab_forcing, p, &
                                     error(method), its(method))
                end do
                rate = 0
                if (n > 2) rate = log(previous/error)/log(2.0_dp)
                print '(3i7,2i5,2es12.3,2f7.2)', p, n, closed%ndof, its, error, rate
               if (n == 8 .and. rate(1) < p + 0.7_dp) error stop 'mortar slab rate loss'
                if (maxval(its) > 60) error stop 'seam p-MG iterations unbounded'
                previous = error
            end do
        end do
    end subroutine

    subroutine check_torus_rates()
        type(tensor_mesh_t) :: mesh, closed
        type(tensor_seam_t) :: seam
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp) :: error, previous, rate
        integer :: p, n, status, its
        print *, 'twisted O-grid torus, mortar: p n dofs its L2 rate area-defect'
        do p = 2, 3
            previous = 0
            do n = 1, 4
                if (n == 3) cycle
                call torus(n, 6*n, p, mesh, master, slave)
                call close_tensor_seam(mesh, master, slave, closed, seam, status, &
                                       motion=identity)
                if (status /= 0) error stop 'torus seam'
               call solve_error(closed, seam, torus_exact, torus_forcing, p, error, its)
                rate = 0
                if (n > 1) rate = log(previous/error)/log(2.0_dp)
                print '(4i7,es12.3,f7.2,es12.3)', p, n, closed%ndof, its, error, rate, &
                    seam%area_defect
                if (n == 4 .and. rate < p + 0.5_dp) error stop 'mortar torus rate loss'
                if (its > 60) error stop 'torus seam p-MG iterations unbounded'
                previous = error
            end do
        end do
    end subroutine

    subroutine check_rigid_twist(n, p)
        !! A shear-free twist keeps unsplit, nearly affine slave faces as large
        !! as curved master cells; its overlaps must still be clipped exactly.
        integer, intent(in) :: n, p
        type(tensor_mesh_t) :: mesh, closed
        type(tensor_seam_t) :: seam
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp), allocatable :: reduced(:), full(:)
        real(dp) :: error, matched
        integer :: status, its
        twist2 = 0
        twist0 = 0.25_dp
        ! A quarter turn maps the O-grid onto itself: an identified seam.
        call torus(n, 6*n, p, mesh, master, slave)
        call close_tensor_seam(mesh, master, slave, closed, seam, status, &
                               motion=identity)
        if (status /= 0) error stop 'quarter-turn seam'
        if (.not. seam%conforming) error stop 'quarter turn not identified'
        call solve_error(closed, seam, torus_exact, torus_forcing, p, matched, its)
        twist0 = 0.37_dp
        call torus(n, 6*n, p, mesh, master, slave)
        call close_tensor_seam(mesh, master, slave, closed, seam, status, &
                               motion=identity)
        if (status /= 0) error stop 'rigid twist seam'
        allocate (reduced(seam%nreduced), full(seam%nfull))
        reduced = 1
        call seam_expand(seam, reduced, full)
        call solve_error(closed, seam, torus_exact, torus_forcing, p, error, its)
        print '(a,2i3,4es12.3)', ' rigid twist n p defect const L2(0.37,0.25)', n, p, &
            seam%area_defect, maxval(abs(full - 1)), error, matched
        twist0 = 0.37_dp
        twist2 = 0.21_dp
        if (seam%area_defect > 1.e-10_dp) error stop 'rigid twist overlap area'
        if (maxval(abs(full - 1)) > 1.e-11_dp) error stop 'rigid twist constants'
        if (error > 1.5_dp*matched) error stop 'rigid twist mortar error'
    end subroutine

    subroutine check_bounded_solves()
        !! p-MG on the constrained spaces converges in a bounded number of steps.
        type(tensor_mesh_t) :: mesh, closed
        type(tensor_seam_t) :: seam
        integer, allocatable :: master(:, :), slave(:, :)
        real(dp) :: error
        integer :: method, status, its
        print *, 'bounded p-MG: geometry method p n its L2'
        do method = SEAM_MORTAR, SEAM_INTERPOLATION
            call slab([4, 4, 4], 2, mesh, master, slave)
            call close_tensor_seam(mesh, master, slave, closed, seam, status, &
                                   motion=translate, images=lattice, method=method)
            if (status /= 0) error stop 'bounded slab seam'
            call solve_error(closed, seam, slab_exact, slab_forcing, 2, error, its)
            print '(4i7,i5,es12.3)', 1, method, 2, 4, its, error
            if (its > 30) error stop 'slab seam p-MG iterations unbounded'
        end do
        call torus(1, 6, 2, mesh, master, slave)
        call close_tensor_seam(mesh, master, slave, closed, seam, status, &
                               motion=identity)
        if (status /= 0) error stop 'bounded torus seam'
        call solve_error(closed, seam, torus_exact, torus_forcing, 2, error, its)
        print '(4i7,i5,es12.3)', 2, SEAM_MORTAR, 2, 1, its, error
        if (its > 30) error stop 'torus seam p-MG iterations unbounded'
    end subroutine
end module
