program test_lab_line_hierarchy
    ! Nested line hierarchies through the shared MF multigrid: B-spline O-grid
    ! h-coarsening (rate p+1, bounded iterations), GLL degree lowering, spline x
    ! Fourier mode coarsening on a torus through the axis and aligned anisotropy
    ! with semi-coarsening.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_line_tables, only: LINE_GLL, LINE_BSPLINE
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_line_hierarchy, only: ogrid_line_hierarchy, box_line_hierarchy, &
                                      PHI_NONE, PHI_FOURIER
    use fortfem_line_transfer, only: gll_line_prolongation
    use fortfem_tensor_operator, only: tensor_load, tensor_boundary_projection, &
                                       tensor_l2_error
    use fortfem_tensor_solver, only: tensor_multigrid_t, tensor_solve_info_t, &
                                     build_line_multigrid, solve_tensor_pcg
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: chi = 1
    logical :: aniso = .false.
    call gll_prolongation_unity()
    call spline_ogrid_rate()
    call gll_ogrid_degree_lowering()
    call spline_fourier_torus()
    call anisotropic_semi_coarsening()
    print '(a)', 'PASS: nested line hierarchies through the shared MF multigrid'
contains
    subroutine gll_prolongation_unity()
        real(dp), allocatable :: p(:, :)
        integer :: status, k
        logical :: wrap
        do k = 1, 2
            wrap = k == 2
            call gll_line_prolongation(4, 2, 3, p, status, periodic=wrap)
            if (status /= 0) error stop 'GLL prolongation'
            if (size(p, 1) /= 12 + merge(0, 1, wrap)) error stop 'GLL fine size'
            if (size(p, 2) /= 6 + merge(0, 1, wrap)) error stop 'GLL coarse size'
            if (maxval(abs(sum(p, 2) - 1)) > 1.e-14_dp) &
                error stop 'GLL prolongation partition of unity'
        end do
    end subroutine

    subroutine solve(meshes, transfers, error, iterations, smoothing)
        type(tensor_mesh_t), intent(in) :: meshes(:)
        type(csr_matrix_t), intent(in) :: transfers(:)
        real(dp), intent(out) :: error
        integer, intent(out) :: iterations
        integer, intent(in), optional :: smoothing
        type(tensor_multigrid_t) :: mg
        type(tensor_solve_info_t) :: info
        real(dp), allocatable :: rhs(:), bc(:), x(:)
        integer :: status
        call build_line_multigrid(meshes, transfers, coefficient, mg, status, &
                                  smoothing=smoothing)
        if (status /= 0) error stop 'line multigrid'
        call tensor_load(mg%level(1)%op, forcing, rhs)
        call tensor_boundary_projection(mg%level(1)%op, exact, bc, status)
        if (status /= 0) error stop 'boundary projection'
        allocate (x(size(rhs)))
        x = 0
        call solve_tensor_pcg(mg, rhs, 1.e-10_dp, x, info, bc=bc, maxiter=200)
        if (info%status /= 0) error stop 'multigrid PCG'
        iterations = info%iterations
        error = tensor_l2_error(mg%level(1)%op, x, exact)
    end subroutine

    subroutine spline_ogrid_rate()
        type(tensor_mesh_t), allocatable :: meshes(:)
        type(csr_matrix_t), allocatable :: transfers(:)
        real(dp) :: error(2)
        integer :: status, k, its(2)
        do k = 1, 2
            call ogrid_line_hierarchy(2, LINE_BSPLINE, 2, 4*k, PHI_NONE, 0, meshes, &
                                      transfers, status)
            if (status /= 0) error stop 'spline O-grid hierarchy'
            if (size(meshes) /= 2 + k) error stop 'spline O-grid level count'
            call solve(meshes, transfers, error(k), its(k))
        end do
        print '(a,2es10.2,a,2i3)', 'spline O-grid p=2 L2:', error, ' its:', its
        if (error(1)/error(2) < 2**2.5_dp) error stop 'spline O-grid rate'
        if (maxval(its) > 10) error stop 'spline O-grid iterations'
    end subroutine

    subroutine gll_ogrid_degree_lowering()
        type(tensor_mesh_t), allocatable :: meshes(:)
        type(csr_matrix_t), allocatable :: transfers(:)
        real(dp) :: error
        integer :: status, its
        call ogrid_line_hierarchy(2, LINE_GLL, 4, 2, PHI_NONE, 0, meshes, &
                                  transfers, status)
        if (status /= 0) error stop 'GLL O-grid hierarchy'
        ! Degrees 4 -> 2 -> 1 on the same cells.
        if (size(meshes) /= 3) error stop 'GLL level count'
        call solve(meshes, transfers, error, its)
        print '(a,es10.2,a,i3)', 'GLL O-grid p=4 L2:', error, ' its:', its
        if (error > 1.e-6_dp) error stop 'GLL O-grid accuracy'
        if (its > 12) error stop 'GLL O-grid iterations'
    end subroutine

    subroutine spline_fourier_torus()
        type(tensor_mesh_t), allocatable :: meshes(:)
        type(csr_matrix_t), allocatable :: transfers(:)
        real(dp) :: error
        integer :: status, its
        call ogrid_line_hierarchy(3, LINE_BSPLINE, 2, 2, PHI_FOURIER, 4, meshes, &
                                  transfers, status, map=torus)
        if (status /= 0) error stop 'spline x Fourier hierarchy'
        call solve(meshes, transfers, error, its)
        print '(a,i3,a,es10.2,a,i3)', 'spline x Fourier torus levels:', &
            size(meshes), ' L2:', error, ' its:', its
        if (its > 15) error stop 'spline x Fourier iterations'
    end subroutine

    subroutine anisotropic_semi_coarsening()
        type(tensor_mesh_t), allocatable :: meshes(:)
        type(csr_matrix_t), allocatable :: transfers(:)
        real(dp) :: error
        integer :: status, its, k
        aniso = .true.
        do k = 1, 2
            chi = merge(1.e3_dp, 1.e6_dp, k == 1)
            call box_line_hierarchy(2, 2, [16, 16, 1], [.true., .false., .false.], &
                                    meshes, transfers, status)
            if (status /= 0) error stop 'box hierarchy'
            call solve(meshes, transfers, error, its)
            print '(a,es9.1,a,es10.2,a,i3)', 'aligned chi=', chi, ' L2:', error, &
                ' its:', its
            if (its > 10) error stop 'semi-coarsening iterations'
        end do
        aniso = .false.
    end subroutine

    subroutine torus(xi, x)
        real(dp), intent(in) :: xi(:)
        real(dp), intent(out) :: x(:)
        real(dp) :: r
        r = 3 + xi(2)
        x = [r*cos(xi(3)), r*sin(xi(3)), xi(1)]
    end subroutine

    function exact(x) result(u)
        real(dp), intent(in) :: x(:)
        real(dp) :: u
        if (aniso) then
            u = sin(pi*x(1))*sin(pi*x(2))
        else if (size(x) == 3) then
            u = sin(x(1)) + cos(x(2)) + sin(x(3))
        else
            u = exp(0.1_dp*x(1) + 0.2_dp*x(2))
        end if
    end function

    function forcing(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        if (aniso) then
            f = (chi + 1)*pi**2*exact(x)
        else if (size(x) == 3) then
            f = exact(x)
        else
            f = -0.05_dp*exact(x)
        end if
    end function

    subroutine coefficient(x, k, c)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k(:, :), c
        integer :: i
        k = 0
        do i = 1, size(x)
            k(i, i) = 1
        end do
        if (aniso) k(1, 1) = chi
        c = 0
    end subroutine
end program
