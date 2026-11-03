program test_meridian_feec_second_kind
    !! Second-kind (edge_kind = 2) meridian complex
    !! P_{k+1} -> (N2_k, P_{k+1}) -> (BDM_k, DG_{k-1}) -> DG_{k-1}:
    !!  - local dimensions of the complete-polynomial spaces,
    !!  - curl_n grad_n = 0 and div_n curl_n = 0 for n = 0, 1, 3, 6, orders
    !!    1-3 and geometry orders 1-3 on curved meshes,
    !!  - the L2 projection of an analytic field into V1 converges at order
    !!    k + 1 and the curl of that projection at order k on a periodically
    !!    perturbed mesh, against the analytic field and curl.
    !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, &
        meridian_mesh_straight
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_op_complex, meridian_zcsc_matmul, meridian_error_norm
    use fortsparse, only: csc_matvec, csc_z_t
    use meridian_test_fields, only: mode_n, vec_a, curl_a, project
    use meridian_test_meshes, only: make_annulus_mesh, &
        make_warped_rectangle_mesh
    implicit none

    integer, parameter :: ns(4) = [0, 1, 3, 6]
    logical :: all_passed

    all_passed = .true.
    !$ call omp_set_num_threads(min(4, omp_get_max_threads()))
    call dimensions()
    call identities()
    call projection_orders()
    call check_summary("Meridian second-kind FEEC complex")
    if (.not. all_passed) error stop 1

contains

    subroutine record(c, d)
        logical, intent(in) :: c
        character(len=*), intent(in) :: d
        all_passed = all_passed .and. c
        call check_condition(c, d)
    end subroutine record

    function max_abs(z) result(m)
        type(csc_z_t), intent(in) :: z
        real(dp) :: m
        m = 0.0_dp
        if (z%nnz > 0) m = maxval(abs(z%val))
    end function max_abs

    subroutine dimensions()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        integer :: k, status, ne, nv, nt
        character(len=40) :: tag

        call make_warped_rectangle_mesh(3, 1, mm)
        nt = mm%topology%n_triangles
        nv = mm%topology%n_vertices
        ! Euler: nv - ne + nt = 1 on the rectangle.
        ne = nv + nt - 1
        do k = 1, 3
            call meridian_feec_build(mm, k, fs, status, edge_kind=2)
            write (tag, '(a,i0)') "order ", k
            call record(status == 0, "second-kind build, "//trim(tag))
            if (status /= 0) cycle
            call record(fs%nln == 2*(k + 1)*(k + 2)/2, &
                "N2_k has complete vector P_k locally, "//trim(tag))
            call record(fs%nl0 == (k + 2)*(k + 3)/2, &
                "V0 is Lagrange P_(k+1), "//trim(tag))
            call record(fs%nn == ne*(k + 1) + nt*(k*k - 1), &
                "k+1 moments per edge and k^2-1 per element, "//trim(tag))
            call record(fs%n0 == nv + ne*k + nt*k*(k - 1)/2, &
                "global P_(k+1) count, "//trim(tag))
        end do
        call meridian_feec_build(mm, 1, fs, status, edge_kind=3)
        call record(status /= 0, "unknown edge family is rejected")
    end subroutine dimensions

    subroutine identities()
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        type(csc_z_t) :: g, c, d, cg, dc
        character(len=80) :: tag
        integer :: geo, k, i, status, shape_id
        real(dp) :: e1, e2, s1, s2

        do shape_id = 1, 2
            do geo = 1, 3
                if (shape_id == 1) then
                    call make_annulus_mesh(2, 8, 0.4_dp, 1.0_dp, 3.0_dp, geo, mm)
                else
                    call make_warped_rectangle_mesh(3, geo, mm)
                end if
                do k = 1, 3
                    call meridian_feec_build(mm, k, fs, status, edge_kind=2)
                    write (tag, '(a,i0,a,i0,a,i0)') "shape ", shape_id, &
                        " geo ", geo, " order ", k
                    call record(status == 0, "build "//trim(tag))
                    if (status /= 0) cycle
                    e1 = 0.0_dp
                    e2 = 0.0_dp
                    s1 = 0.0_dp
                    s2 = 0.0_dp
                    do i = 1, size(ns)
                        call meridian_op_complex(fs%grad, ns(i), g, status)
                        call meridian_op_complex(fs%curl, ns(i), c, status)
                        call meridian_op_complex(fs%div, ns(i), d, status)
                        call meridian_zcsc_matmul(c, g, cg, status)
                        call meridian_zcsc_matmul(d, c, dc, status)
                        e1 = max(e1, max_abs(cg))
                        e2 = max(e2, max_abs(dc))
                        s1 = max(s1, max_abs(g))
                        s2 = max(s2, max_abs(c))
                    end do
                    call record(e1 < 1.0e-12_dp*max(1.0_dp, s1*s2), &
                        "curl_n grad_n = 0 for n = 0,1,3,6, "//trim(tag))
                    call record(e2 < 1.0e-12_dp*max(1.0_dp, s2), &
                        "div_n curl_n = 0 for n = 0,1,3,6, "//trim(tag))
                    call record(s1 > 0.5_dp .and. s2 > 0.5_dp, &
                        "operators are nontrivial, "//trim(tag))
                    print '(a,a,2es10.2)', "  residuals ", trim(tag), e1, e2
                end do
            end do
        end do
    end subroutine identities

    subroutine perturbed_mesh(nx, mm)
        !! [2, 3] x [-0.5, 0.5] with a repeated two-by-two vertex displacement
        !! pattern, so that no refinement becomes locally translation
        !! invariant.
        integer, intent(in) :: nx
        type(meridian_mesh_t), intent(out) :: mm
        type(mesh_2d_t) :: mesh
        real(dp) :: shift(2, 0:1, 0:1), h
        integer :: v, ix, iy, status

        shift(:, 0, 0) = [0.0_dp, 0.0_dp]
        shift(:, 1, 0) = [0.13_dp, -0.11_dp]
        shift(:, 0, 1) = [-0.09_dp, 0.16_dp]
        shift(:, 1, 1) = [0.18_dp, 0.07_dp]
        call mesh%create_rectangular(nx + 1, nx + 1, 2.0_dp, 3.0_dp, &
            -0.5_dp, 0.5_dp)
        h = 1.0_dp/nx
        do v = 1, mesh%n_vertices
            ix = nint((mesh%vertices(1, v) - 2.0_dp)/h)
            iy = nint((mesh%vertices(2, v) + 0.5_dp)/h)
            if (ix == 0 .or. ix == nx .or. iy == 0 .or. iy == nx) cycle
            mesh%vertices(:, v) = mesh%vertices(:, v) &
                + h*shift(:, modulo(ix, 2), modulo(iy, 2))
        end do
        call meridian_mesh_straight(mesh, 1, mm, status)
        if (status /= 0) error stop "perturbed mesh failed"
    end subroutine perturbed_mesh

    subroutine projection_errors(nx, k, edge_kind, ea, ec)
        integer, intent(in) :: nx, k, edge_kind
        real(dp), intent(out) :: ea, ec
        type(meridian_mesh_t) :: mm
        type(meridian_feec_t) :: fs
        complex(dp), allocatable :: u(:), w(:)
        logical :: ok
        integer :: status

        call perturbed_mesh(nx, mm)
        call meridian_feec_build(mm, k, fs, status, edge_kind=edge_kind)
        if (status /= 0) error stop "projection build failed"
        call project(fs, 1, vec_a, u, ok)
        if (.not. ok) error stop "projection failed"
        call meridian_error_norm(fs, 1, u, vec_a, ea, status)
        if (status /= 0) error stop "projection error failed"
        w = cmplx(csc_matvec(fs%curl%a, real(u, dp)), &
            csc_matvec(fs%curl%a, aimag(u)), dp)
        call meridian_error_norm(fs, 2, w, curl_a, ec, status)
        if (status /= 0) error stop "curl error failed"
    end subroutine projection_errors

    subroutine projection_orders()
        integer, parameter :: nx(2) = [8, 16]
        real(dp) :: ea(2), ec(2), ra, rc, fa(2), fc(2)
        character(len=40) :: tag
        integer :: k, lvl

        mode_n = 0
        do k = 1, 2
            do lvl = 1, 2
                call projection_errors(nx(lvl), k, 2, ea(lvl), ec(lvl))
                call projection_errors(nx(lvl), k, 1, fa(lvl), fc(lvl))
            end do
            ra = log(ea(1)/ea(2))/log(2.0_dp)
            rc = log(ec(1)/ec(2))/log(2.0_dp)
            write (tag, '(a,i0)') "order ", k
            print '(a,a,2f7.2,a,2f7.2)', "  projection L2/curl rates N2 ", &
                trim(tag), ra, rc, "  N1 ", log(fa(1)/fa(2))/log(2.0_dp), &
                log(fc(1)/fc(2))/log(2.0_dp)
            call record(ra > real(k + 1, dp) - 0.3_dp, &
                "L2 projection into N2_k is order k+1, "//trim(tag))
            call record(rc > real(k, dp) - 0.3_dp, &
                "curl of the N2_k projection is order k, "//trim(tag))
            call record(ea(2) < fa(2), &
                "N2_k projection is closer than N1_k, "//trim(tag))
        end do
    end subroutine projection_orders

end program test_meridian_feec_second_kind
