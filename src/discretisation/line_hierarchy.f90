module fortfem_line_hierarchy
    ! Nested hierarchies for the line-table MF operator. O-grid: five conforming
    ! patches through the axis, times a periodic phi line of the same family or
    ! a Fourier phi line. B-spline levels halve the patch spans and phi spans
    ! (periodic splines keep at least degree+1 spans) and Fourier levels halve
    ! the modes; GLL levels lower the degree p -> (p+1)/2 on the same cells, as
    ! MF1's p-multigrid. Box: one structured B-spline patch with per-direction
    ! (semi-)coarsening. Geometry degree and the phi rule stay fixed per level.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_line_tables, only: tensor_line_tables_t, LINE_GLL, &
                                          LINE_BSPLINE
    use fortfem_tensor_mesh, only: tensor_mesh_t, geometry_map
    use fortfem_bspline_lines, only: build_bspline_line
    use fortfem_gll_line, only: build_gll_line
    use fortfem_fourier_lines, only: build_fourier_line
    use fortfem_line_mesh, only: line_ogrid_mesh, line_tensor_mesh, &
                                 line_mesh_prolongation, line_matrix_t
    use fortfem_line_transfer, only: bspline_line_prolongation, &
                                     fourier_line_prolongation, gll_line_prolongation
    use fortfem_element_to_csr, only: csr_matrix_t
    implicit none
    private
    public :: ogrid_line_hierarchy, box_line_hierarchy
    public :: PHI_NONE, PHI_LOCAL, PHI_FOURIER
    ! PHI_LOCAL: periodic phi line of the patch family (B-spline or GLL).
    integer, parameter :: PHI_NONE = 0, PHI_LOCAL = 1, PHI_FOURIER = 2
    real(dp), parameter :: pi = acos(-1.0_dp)
contains
    subroutine ogrid_line_hierarchy(dim, family, degree, n, phi, nphi, meshes, &
                                    transfers, status, map, continuity, q, nq, &
                                    nq_phi, levels)
        ! n: patch spans per direction (power of two for B-splines). nphi: phi
        ! cells (PHI_LOCAL) or Fourier modes (PHI_FOURIER). map(xi, x) receives
        ! (rx, ry, phi) in 3D and (rx, ry) in 2D.
        integer, intent(in) :: dim, family, degree, n, phi, nphi
        type(tensor_mesh_t), allocatable, intent(out) :: meshes(:)
        type(csr_matrix_t), allocatable, intent(out) :: transfers(:)
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        integer, intent(in), optional :: continuity, q, nq, nq_phi, levels
        type(tensor_line_tables_t) :: lines(3)
        type(line_matrix_t) :: p1d(3)
        integer :: nlev, l, nb(32), np(32), pd(32), ratio(3), k, i, rule, qg
        status = 1
        if (dim == 2 .and. phi /= PHI_NONE) return
        if (dim == 3 .and. phi == PHI_NONE) return
        if (family /= LINE_GLL .and. family /= LINE_BSPLINE) return
        if (n < 1) return
        if (family == LINE_BSPLINE .and. iand(n, n - 1) /= 0) return
        k = degree - 1
        if (present(continuity)) k = continuity
        qg = degree
        if (present(q)) qg = q
        rule = 2*nphi + 1
        if (phi == PHI_LOCAL) rule = degree + 1
        if (present(nq_phi)) rule = nq_phi
        nb(1) = n; np(1) = nphi; pd(1) = degree
        nlev = 1
        do while (nlev < 32)
            if (present(levels)) then
                if (nlev >= levels) exit
            end if
            nb(nlev + 1) = nb(nlev); np(nlev + 1) = np(nlev); pd(nlev + 1) = pd(nlev)
            if (family == LINE_GLL) then
                pd(nlev + 1) = (pd(nlev) + 1)/2
                if (phi == PHI_FOURIER) np(nlev + 1) = np(nlev)/2
            else
                nb(nlev + 1) = max(1, nb(nlev)/2)
                if (phi == PHI_LOCAL) then
                    if (np(nlev)/2 >= degree + 1 .and. modulo(np(nlev), 2) == 0) &
                        np(nlev + 1) = np(nlev)/2
                else if (phi == PHI_FOURIER) then
                    np(nlev + 1) = np(nlev)/2
                end if
            end if
            if (all([nb(nlev + 1), np(nlev + 1), pd(nlev + 1)] == &
                    [nb(nlev), np(nlev), pd(nlev)])) exit
            nlev = nlev + 1
        end do
        allocate (meshes(nlev), transfers(nlev - 1))
        do l = 1, nlev
            if (family == LINE_GLL) then
                call build_gll_line(pd(l), nb(l), lines(1), status, x0=0.0_dp, &
                                    x1=1.0_dp)
                if (status /= 0) return
                if (phi == PHI_LOCAL) then
                    call build_gll_line(pd(l), np(l), lines(3), status, &
                                        periodic=.true., x0=0.0_dp, x1=2*pi)
                    if (status /= 0) return
                end if
            else
                call build_bspline_line(degree, lines(1), status, nelem=nb(l), &
                                        continuity=k, nq=nq)
                if (status /= 0) return
                if (phi == PHI_LOCAL) then
                    call build_bspline_line(degree, lines(3), status, periodic=.true., &
                                            knots=[(2*pi*i/np(l), i=0, np(l))], &
                                            continuity=k, nq=rule)
                    if (status /= 0) return
                end if
            end if
            if (phi == PHI_FOURIER) then
                call build_fourier_line(np(l), lines(3), status, nq=rule)
                if (status /= 0) return
            end if
            call line_ogrid_mesh(dim, lines, meshes(l), status, map=map, q=qg)
            if (status /= 0) return
        end do
        do l = 1, nlev - 1
            ratio = 1
            if (family == LINE_GLL) then
                call gll_line_prolongation(pd(l), pd(l + 1), nb(l), p1d(1)%a, status)
                if (status /= 0) return
                if (phi == PHI_LOCAL) then
                    call gll_line_prolongation(pd(l), pd(l + 1), np(l), p1d(3)%a, &
                                               status, periodic=.true.)
                    if (status /= 0) return
                end if
            else
                if (nb(l) /= nb(l + 1)) then
                    ratio(1:2) = 2
                    call bspline_line_prolongation(degree, nb(l + 1), p1d(1)%a, &
                                                   status, continuity=k)
                    if (status /= 0) return
                else
                    call identity(meshes(l)%line(1)%ndof, p1d(1)%a)
                end if
                if (phi == PHI_LOCAL) then
                    if (np(l) /= np(l + 1)) then
                        ratio(3) = 2
                        call bspline_line_prolongation(degree, np(l + 1), p1d(3)%a, &
                                                       status, continuity=k, &
                                                       periodic=.true.)
                        if (status /= 0) return
                    else
                        call identity(meshes(l)%line(3)%ndof, p1d(3)%a)
                    end if
                end if
            end if
            p1d(2)%a = p1d(1)%a
            if (phi == PHI_FOURIER) then
                call fourier_line_prolongation(np(l), np(l + 1), p1d(3)%a, status)
                if (status /= 0) return
            else if (phi == PHI_NONE) then
                call identity(1, p1d(3)%a)
            end if
            call line_mesh_prolongation(meshes(l), meshes(l + 1), p1d, ratio, &
                                        transfers(l), status)
            if (status /= 0) return
        end do
        status = 0
    end subroutine

    subroutine box_line_hierarchy(dim, degree, n, coarsen, meshes, transfers, status, &
                                  map, continuity, levels)
        ! One B-spline patch with n(d) spans on [0, 1]^dim. Directions with
        ! coarsen(d) halve their spans per level (semi-coarsening otherwise).
        integer, intent(in) :: dim, degree, n(3)
        logical, intent(in) :: coarsen(3)
        type(tensor_mesh_t), allocatable, intent(out) :: meshes(:)
        type(csr_matrix_t), allocatable, intent(out) :: transfers(:)
        integer, intent(out) :: status
        procedure(geometry_map), optional :: map
        integer, intent(in), optional :: continuity, levels
        type(tensor_line_tables_t) :: lines(3)
        type(line_matrix_t) :: p1d(3)
        integer :: nlev, l, d, k, nc(3, 32), ratio(3)
        status = 1
        k = degree - 1
        if (present(continuity)) k = continuity
        nc(:, 1) = n
        nlev = 1
        do while (nlev < 32)
            if (present(levels)) then
                if (nlev >= levels) exit
            end if
            nc(:, nlev + 1) = nc(:, nlev)
            do d = 1, dim
                if (coarsen(d) .and. modulo(nc(d, nlev), 2) == 0) &
                    nc(d, nlev + 1) = nc(d, nlev)/2
            end do
            if (all(nc(:, nlev + 1) == nc(:, nlev))) exit
            nlev = nlev + 1
        end do
        allocate (meshes(nlev), transfers(nlev - 1))
        do l = 1, nlev
            do d = 1, dim
                call build_bspline_line(degree, lines(d), status, nelem=nc(d, l), &
                                        continuity=k)
                if (status /= 0) return
            end do
            call line_tensor_mesh(dim, lines, meshes(l), status, map=map)
            if (status /= 0) return
        end do
        do l = 1, nlev - 1
            ratio = 1
            do d = 1, 3
                if (d > dim) then
                    call identity(1, p1d(d)%a)
                else if (nc(d, l) /= nc(d, l + 1)) then
                    ratio(d) = 2
                    call bspline_line_prolongation(degree, nc(d, l + 1), p1d(d)%a, &
                                                   status, continuity=k)
                    if (status /= 0) return
                else
                    call identity(meshes(l)%line(d)%ndof, p1d(d)%a)
                end if
            end do
            call line_mesh_prolongation(meshes(l), meshes(l + 1), p1d, ratio, &
                                        transfers(l), status)
            if (status /= 0) return
        end do
        status = 0
    end subroutine

    subroutine identity(m, a)
        integer, intent(in) :: m
        real(dp), allocatable, intent(out) :: a(:, :)
        integer :: j
        allocate (a(m, m))
        a = 0
        do j = 1, m
            a(j, j) = 1
        end do
    end subroutine
end module
