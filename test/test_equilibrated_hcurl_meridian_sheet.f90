program test_equilibrated_hcurl_meridian_sheet
    !! Edge (sheet) loads and the free-trace hook of the meridian per-n
    !! H(curl) equilibrated bound, by exact recovery on the straight mesh of
    !! [2, 3] x [-0.5, 0.5]: A_h is the projection of A = A_0 + chi(Z > 0) c Z
    !! e_R (A_0 a degree-2 polynomial field), which lies in V1 (k = 3) since
    !! Z = 0 is a mesh line; rot A jumps by -c across it, i.e. a sheet load
    !! int G (v . e_R) dl with G = -R c on the edges Z = 0 oriented along +R.
    !!   1. fixed trace (v x n = 0), sheet load: the bound vanishes; with the
    !!      sign of G reversed it does not;
    !!   2. free trace, boundary loads G_b = R rot A (the trace sigma = R B_phi
    !!      of the exact field) together with the sheet: the bound vanishes;
    !!      with G_b = 0 it does not; the side default (flood fill) and an
    !!      explicit side marker agree;
    !!   3. invalid requests (boundary loads without free trace, patchwise
    !!      method with loads, a non-separating sheet) are rejected.
    use check, only: check_condition, check_summary
    use fortfem_equilibrated_hcurl_meridian_n, only: HCURL_BAD_INPUT, &
        HCURL_GLOBAL, HCURL_OK, HCURL_PATCHWISE, &
        meridian_hcurl_equilibrated_bound, meridian_hcurl_estimate_t
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, &
        meridian_mesh_from_map
    use fortfem_meridian_feec_n, only: meridian_feec_build, meridian_feec_t
    use meridian_test_fields, only: project
    implicit none

    complex(dp), parameter :: ca1 = (1.0_dp, 0.5_dp), ca2 = (0.8_dp, -0.3_dp)
    complex(dp), parameter :: csheet = (0.7_dp, 0.4_dp)
    integer, parameter :: mode = 2
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs, fsd
    type(meridian_hcurl_estimate_t) :: est, est2
    complex(dp), allocatable :: u(:)
    integer, allocatable :: sheet_cell(:), sheet_edge(:), all_cell(:), all_edge(:)
    logical, allocatable :: upper(:)
    real(dp) :: sign_g, boundary_scale
    integer :: status
    logical :: ok

    call straight_rectangle(4, mm)
    call meridian_feec_build(mm, 3, fs, status, quad_degree=14)
    call check_condition(status == 0, "primal build")
    call meridian_feec_build(mm, 3, fsd, status, quad_degree=14)
    call check_condition(status == 0, "dual build")
    call project(fs, 1, piecewise_a, u, ok)
    call check_condition(ok, "projection")
    call find_edges()

    ! 1. fixed trace, interior sheet
    sign_g = 1.0_dp
    boundary_scale = 1.0_dp
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, edge_cell=sheet_cell, edge_local=sheet_edge, edge_load=load)
    print '(a, 4es11.3)', "fixed trace, sheet: bound, |B_h|, defect, mismatch ", &
        est%bound, est%field_norm, est%load_defect, est%lift_mismatch
    call check_condition(status == HCURL_OK, "sheet status")
    call check_condition(est%bound < 1.0e-9_dp*est%field_norm, &
        "sheet load: bound vanishes for the exact field")
    call check_condition(est%scaled_potential .and. .not. est%trace_fixed, &
        "sheet flags")
    call check_condition(est%load_defect < 1.0e-12_dp, &
        "linear sheet data are interpolated exactly")
    sign_g = -1.0_dp
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, edge_cell=sheet_cell, edge_local=sheet_edge, edge_load=load)
    print '(a, es11.3)', "fixed trace, reversed sheet: bound ", est%bound
    call check_condition(status == HCURL_OK .and. est%bound > &
        1.0e-3_dp*est%field_norm, "a reversed sheet load is detected")
    sign_g = 1.0_dp

    ! 2. free trace: sheet plus boundary traces
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, edge_cell=all_cell, edge_local=all_edge, edge_load=load, &
        free_trace=.true.)
    print '(a, 3es11.3)', "free trace: bound, |B_h|, mismatch ", est%bound, &
        est%field_norm, est%lift_mismatch
    call check_condition(status == HCURL_OK, "free trace status")
    call check_condition(est%bound < 1.0e-9_dp*est%field_norm, &
        "free trace: bound vanishes for the exact field and traces")
    call check_condition(est%trace_fixed, "free trace flag")
    call check_condition(est%lift_mismatch < 1.0e-12_dp, &
        "consistent traces at the sheet end points")
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est2, &
        status, edge_cell=all_cell, edge_local=all_edge, edge_load=load, &
        free_trace=.true., side=upper)
    call check_condition(status == HCURL_OK .and. abs(est2%bound - est%bound) &
        <= 1.0e-12_dp*est%field_norm, "explicit side marker agrees")
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est2, &
        status, edge_cell=all_cell, edge_local=all_edge, edge_load=load, &
        free_trace=.true., side=.not. upper)
    call check_condition(status == HCURL_OK .and. abs(est2%bound - est%bound) &
        <= 1.0e-12_dp*est%field_norm, "complementary side marker agrees")
    boundary_scale = 0.0_dp
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, edge_cell=all_cell, edge_local=all_edge, edge_load=load, &
        free_trace=.true.)
    print '(a, es11.3)', "free trace, zero boundary loads: bound ", est%bound
    call check_condition(status == HCURL_OK .and. est%bound > &
        1.0e-3_dp*est%field_norm, "missing boundary loads are detected")
    boundary_scale = 1.0_dp

    ! 3. rejected requests
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, edge_cell=all_cell, edge_local=all_edge, edge_load=load)
    call check_condition(status == HCURL_BAD_INPUT, &
        "boundary loads need the free trace")
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, HCURL_PATCHWISE, sheet_cell, sheet_edge, load)
    call check_condition(status == HCURL_BAD_INPUT, "loads need the global method")
    call meridian_hcurl_equilibrated_bound(fs, fsd, mode, u, piecewise_j, est, &
        status, HCURL_GLOBAL, sheet_cell(2:), sheet_edge(2:), load)
    call check_condition(status == HCURL_BAD_INPUT, &
        "a sheet that does not separate the mesh is rejected")
    call check_summary("Meridian H(curl) bound with sheet and trace loads")

contains

    subroutine straight_rectangle(nx, mesh)
        integer, intent(in) :: nx
        type(meridian_mesh_t), intent(out) :: mesh
        type(mesh_2d_t) :: topo
        real(dp), allocatable :: pv(:, :, :)
        integer :: t, st

        call topo%create_rectangular(nx + 1, nx + 1, 2.0_dp, 3.0_dp, -0.5_dp, &
            0.5_dp)
        allocate (pv(2, 3, topo%n_triangles))
        do t = 1, topo%n_triangles
            pv(:, :, t) = topo%vertices(:, topo%triangles(:, t))
        end do
        call meridian_mesh_from_map(topo, pv, 1, mesh, st)
        if (st /= 0) error stop "rectangle mesh"
    end subroutine straight_rectangle

    subroutine find_edges()
        !! Sheet edges Z = 0 seen from the upper cell (oriented along +R) and
        !! all boundary edges.
        real(dp) :: p(2, 3), centre(2)
        integer :: t, e, nt
        integer, allocatable :: bc(:), be(:)

        nt = mm%topology%n_triangles
        allocate (sheet_cell(0), sheet_edge(0), bc(0), be(0), upper(nt))
        do t = 1, nt
            p = mm%topology%vertices(:, mm%topology%triangles(1:3, t))
            centre = sum(p, dim=2)/3
            upper(t) = centre(2) > 0
            do e = 1, 3
                associate (a => p(:, e), b => p(:, mod(e, 3) + 1))
                    if (abs(a(2)) < 1.0e-12_dp .and. abs(b(2)) < 1.0e-12_dp &
                        .and. upper(t)) then
                        sheet_cell = [sheet_cell, t]
                        sheet_edge = [sheet_edge, e]
                    end if
                    if (same_side(a, b)) then
                        bc = [bc, t]
                        be = [be, e]
                    end if
                end associate
            end do
        end do
        call check_condition(size(sheet_cell) == 4, "four sheet edges")
        call check_condition(size(bc) == 16, "sixteen boundary edges")
        all_cell = [sheet_cell, bc]
        all_edge = [sheet_edge, be]
    end subroutine find_edges

    logical function same_side(a, b)
        !! Both points on the same side of the rectangle boundary.
        real(dp), intent(in) :: a(2), b(2)
        real(dp), parameter :: tol = 1.0e-12_dp

        same_side = (abs(a(1) - 2) < tol .and. abs(b(1) - 2) < tol) .or. &
            (abs(a(1) - 3) < tol .and. abs(b(1) - 3) < tol) .or. &
            (abs(a(2) + 0.5_dp) < tol .and. abs(b(2) + 0.5_dp) < tol) .or. &
            (abs(a(2) - 0.5_dp) < tol .and. abs(b(2) - 0.5_dp) < tol)
    end function same_side

    subroutine piecewise_a(x, v)
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        real(dp) :: r, z

        r = x(1)
        z = x(2)
        v(1) = ca1*(1.0_dp + 0.5_dp*z + 0.3_dp*r + 0.2_dp*r*z - 0.4_dp*z*z)
        v(2) = (0.0_dp, 0.0_dp)
        v(3) = ca2*(-0.7_dp + 0.4_dp*r + 0.6_dp*z + 0.3_dp*r*r - 0.5_dp*r*z)
        if (z > 0) v(1) = v(1) + csheet*z
    end subroutine piecewise_a

    complex(dp) function rot_a(x) result(w)
        !! rot A = d_R A_Z - d_Z A_R of piecewise_a (side by Z, upper for Z = 0
        !! taken from the caller).
        real(dp), intent(in) :: x(2)

        w = ca2*(0.4_dp + 0.6_dp*x(1) - 0.5_dp*x(2)) &
            - ca1*(0.5_dp + 0.2_dp*x(1) - 0.8_dp*x(2))
    end function rot_a

    subroutine piecewise_j(x, v)
        !! mu0 J_pol = (1/R) (-rot90 grad(R w) + n^2 A / R), w = rot A, on each
        !! side of Z = 0 (the extra field c Z e_R has w = -c).
        real(dp), intent(in) :: x(2)
        complex(dp), intent(out) :: v(3)
        complex(dp) :: w, wr, wz, a(3)
        real(dp) :: r, nn

        r = x(1)
        nn = real(mode, dp)
        call piecewise_a(x, a)
        w = rot_a(x)
        wr = ca2*(0.6_dp) - ca1*(0.2_dp)
        wz = ca2*(-0.5_dp) + ca1*(0.8_dp)
        if (x(2) > 0) w = w - csheet
        v(1) = (r*wz + nn*nn*a(1)/r)/r
        v(2) = (0.0_dp, 0.0_dp)
        v(3) = (-(w + r*wr) + nn*nn*a(3)/r)/r
    end subroutine piecewise_j

    subroutine load(cell, edge, tau, g)
        !! Sheet: G = -R c on Z = 0 (oriented along +R, left = upper);
        !! boundary: G_b = R rot A of the cell's side.
        integer, intent(in) :: cell, edge
        real(dp), intent(in) :: tau
        complex(dp), intent(out) :: g
        real(dp) :: p(2, 3), x(2)
        complex(dp) :: w

        p = mm%topology%vertices(:, mm%topology%triangles(1:3, cell))
        x = (1.0_dp - tau)*p(:, edge) + tau*p(:, mod(edge, 3) + 1)
        if (abs(p(2, edge)) < 1.0e-12_dp .and. &
            abs(p(2, mod(edge, 3) + 1)) < 1.0e-12_dp) then
            g = -sign_g*x(1)*csheet
        else
            w = rot_a(x)
            if (upper(cell)) w = w - csheet
            g = boundary_scale*x(1)*w
        end if
    end subroutine load

end program test_equilibrated_hcurl_meridian_sheet
