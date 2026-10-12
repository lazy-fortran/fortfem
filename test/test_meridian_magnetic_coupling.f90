program test_meridian_magnetic_coupling
    !! Boundary chain and normal-flux pairing of the per-n FEEC complex.
    !!  1. The chain is one closed counter-clockwise loop whose panel nodes
    !!     are the boundary vertices and edge mid nodes (signed area of the
    !!     node polygon equals the mesh area up to the curved-edge segments).
    !!  2. Partition of unity of the P2 trace basis and the divergence
    !!     theorem: sum_i C2(i, j) = int_Omega div q_j dR dZ for every RT
    !!     basis field, an exact discrete identity on curved meshes.
    !!  3. Point evaluation reproduces the quadrature-point evaluation.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_qpoint, meridian_point_basis_t, meridian_point_basis_init, &
        meridian_point_basis
    use fortfem_meridian_magnetic_coupling, only: meridian_boundary_t, &
        meridian_boundary_build, meridian_flux_pairing
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use meridian_test_meshes, only: make_warped_rectangle_mesh
    implicit none

    logical :: all_passed
    integer :: k, g

    all_passed = .true.
    do k = 1, 3
        do g = 1, 2
            call run(k, g)
        end do
    end do
    call check_summary("Meridian FEEC magnetic coupling")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    subroutine run(k, g)
        integer, intent(in) :: k, g
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(meridian_boundary_t) :: bd
        type(meridian_point_basis_t) :: pb
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: vals(:), colsum(:), divint(:), phi(:, :)
        real(dp), allocatable :: phi2(:, :), xi(:), eta(:), w(:)
        real(dp) :: area, x(2), x2(2), wq, err, scale, v
        integer :: status, p, t, q, j, l, nn, nt
        character(len=40) :: tag

        write (tag, '(a,i0,a,i0)') "order ", k, ", geometry ", g
        call make_warped_rectangle_mesh(6, g, mm)
        call meridian_feec_build(mm, k, fs, status)
        call record(status == 0, "build, "//trim(tag))
        call meridian_boundary_build(fs, bd, status)
        call record(status == 0 .and. bd%n_panels == 24, &
            "closed boundary chain of 24 edges, "//trim(tag))
        area = 0.0_dp
        do p = 1, size(bd%curve%nodes, 2)
            l = mod(p, size(bd%curve%nodes, 2)) + 1
            area = area + bd%curve%nodes(1, p)*bd%curve%nodes(2, l) &
                - bd%curve%nodes(1, l)*bd%curve%nodes(2, p)
        end do
        call record(area > 0.0_dp, "counter-clockwise chain, "//trim(tag))

        ! divergence theorem
        call meridian_flux_pairing(fs, bd, rows, cols, vals, status)
        call record(status == 0, "pairing assembled, "//trim(tag))
        nn = fs%nn
        allocate (colsum(fs%sp(2)%ndof), divint(fs%sp(2)%ndof))
        colsum = 0.0_dp
        do p = 1, size(vals)
            colsum(cols(p)) = colsum(cols(p)) + vals(p)
        end do
        ! int div q_j dR dZ = int (D q_j) phys3 R dR dZ, phys3 = g / R
        divint = 0.0_dp
        nt = fs%mm%topology%n_triangles
        allocate (phi(3, fs%sp(3)%nloc))
        do j = 1, fs%div%a%ncol
            do p = fs%div%a%col_ptr(j), fs%div%a%col_ptr(j + 1) - 1
                l = fs%div%a%row_idx(p)
                t = (l - 1)/fs%nld + 1
                v = 0.0_dp
                do q = 1, fs%nq
                    call meridian_qpoint(fs, 3, t, q, x, wq, phi, status)
                    v = v + wq*phi(1, l - (t - 1)*fs%nld)
                end do
                divint(j) = divint(j) + fs%div%a%val(p)*v
            end do
        end do
        scale = maxval(abs(divint(1:nn)))
        err = maxval(abs(colsum(1:nn) - divint(1:nn)))/scale
        print '(a,a,es10.2)', "  ", trim(tag)//": divergence identity", err
        call record(err < 1.0e-11_dp, "divergence theorem, "//trim(tag))

        ! point evaluation against quadrature tables
        call meridian_point_basis_init(fs, pb, status)
        call triangle_duffy_quadrature(2*k + 2*fs%g + 2, xi, eta, w, status)
        allocate (phi2(3, fs%sp(1)%nloc))
        deallocate (phi)
        allocate (phi(3, fs%sp(1)%nloc))
        err = 0.0_dp
        do t = 1, nt, 7
            do q = 1, size(xi), 3
                call meridian_qpoint(fs, 1, t, q, x, wq, phi, status)
                call meridian_point_basis(fs, pb, 1, t, xi(q), eta(q), x2, &
                    phi2, status)
                err = max(err, maxval(abs(phi - phi2)), maxval(abs(x - x2)))
            end do
        end do
        call record(err < 1.0e-12_dp, "point evaluation, "//trim(tag))
    end subroutine run

end program test_meridian_magnetic_coupling
