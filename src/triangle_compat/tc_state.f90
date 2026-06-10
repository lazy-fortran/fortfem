module tc_state
    !! Mesh state and handle primitives for the Triangle-compatible mesher.
    !! Data layout mirrors the oriented-triangle/subsegment structure of
    !! Shewchuk's published mesh generator so that traversal and allocation
    !! order reproduce its behavior exactly.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    implicit none
    private

    public :: tc_state_t, tc_init_state, tc_random
    public :: t_slot, t_orient, t_make, s_slot, s_so, s_make
    public :: t_sym, t_lnext, t_lprev, t_onext, t_oprev, t_dnext, t_dprev
    public :: t_org, t_dest, t_apex, set_org, set_dest, set_apex
    public :: t_bond, t_dissolve, ts_pivot, ts_bond, ts_dissolve
    public :: st_pivot, st_dissolve, s_sym, s_pivot, s_next
    public :: s_org, s_dest, set_s_org, set_s_dest
    public :: seg_org, seg_dest, set_seg_org, set_seg_dest
    public :: s_bond, s_dissolve
    public :: V_INPUT, V_SEGMENT, V_FREE, V_DEAD, V_UNDEAD
    public :: TRIPERBLOCK, SAMPLEFACTOR, SQRT2_TC

    integer, parameter :: V_INPUT = 0, V_SEGMENT = 1, V_FREE = 2
    integer, parameter :: V_DEAD = 3, V_UNDEAD = 4
    integer, parameter :: TRIPERBLOCK = 4092
    integer, parameter :: SAMPLEFACTOR = 11
    real(dp), parameter :: SQRT2_TC = 1.4142135623730951_dp
    integer, parameter :: P1(0:2) = [1, 2, 0]
    integer, parameter :: M1(0:2) = [2, 0, 1]

    type :: tc_state_t
        ! vertices: slot 1.. (0 = null vertex)
        real(dp), allocatable :: vx(:), vy(:)
        integer, allocatable :: vmark(:), vtype(:), vtri(:)
        integer :: nv_slots = 0, nv_live = 0
        integer, allocatable :: vfree(:)
        integer :: nv_free = 0
        ! triangles: slot 0 = outer space (dummy)
        integer, allocatable :: tnbr(:, :), tver(:, :), tsub(:, :)
        logical, allocatable :: tdead(:), tinf(:)
        integer :: nt_slots = 0, nt_live = 0
        integer, allocatable :: tfree(:)
        integer :: nt_free = 0
        ! subsegments: slot 0 = omnipresent dummy
        integer, allocatable :: snbr(:, :), sver(:, :), ssegv(:, :), stri(:, :)
        integer, allocatable :: smark(:)
        logical, allocatable :: sdead(:)
        integer :: ns_slots = 0, ns_live = 0
        integer, allocatable :: sfree(:)
        integer :: ns_free = 0
        ! behavior
        real(dp) :: minangle = 0.0_dp, goodangle = 0.0_dp, offconstant = 0.0_dp
        real(dp) :: maxarea = -1.0_dp
        integer :: nobisect = 0
        logical :: fixedarea = .false., quality = .false.
        integer :: steinerleft = -1
        ! mesh status
        integer :: invertices = 0, insegments = 0
        integer(int64) :: hullsize = 0
        logical :: checksegments = .false., checkquality = .false.
        integer :: recent = -1
        integer :: samples = 1
        integer(int64) :: rngseed = 1
        real(dp) :: xmin = 0.0_dp, xmax = 0.0_dp, ymin = 0.0_dp, ymax = 0.0_dp
        integer :: undeads = 0
        ! flip journal (one vertex insertion)
        integer, allocatable :: jtri(:), jkind(:)
        integer :: jn = 0
        ! virus pool
        integer, allocatable :: vir(:)
        integer :: vir_n = 0
        ! bad triangle priority queues (4096 FIFO buckets)
        integer, allocatable :: btq_front(:), btq_tail(:), btq_nextq(:)
        integer :: btq_first = -1
        integer, allocatable :: bt_tri(:), bt_org(:), bt_dest(:), bt_apex(:)
        integer, allocatable :: bt_next(:)
        real(dp), allocatable :: bt_key(:)
        integer :: bt_slots = 0, bt_items = 0
        integer, allocatable :: bt_free(:)
        integer :: bt_nfree = 0
        ! bad (encroached) subsegment pool
        integer, allocatable :: bs_sub(:), bs_org(:), bs_dest(:)
        logical, allocatable :: bs_dead(:)
        integer :: bs_slots = 0, bs_items = 0
        integer, allocatable :: bs_free(:)
        integer :: bs_nfree = 0
    end type tc_state_t

contains

    subroutine tc_init_state(s, npoints)
        type(tc_state_t), intent(out) :: s
        integer, intent(in) :: npoints
        integer :: cap

        cap = max(4*npoints, 64)
        allocate (s%vx(cap), s%vy(cap), s%vmark(cap), s%vtype(cap), s%vtri(cap))
        allocate (s%vfree(cap))
        s%vtri = -1
        cap = max(8*npoints, 64)
        allocate (s%tnbr(0:2, 0:cap), s%tver(0:2, 0:cap), s%tsub(0:2, 0:cap))
        allocate (s%tdead(0:cap), s%tinf(0:cap), s%tfree(cap))
        s%tnbr(:, 0) = 0
        s%tver(:, 0) = 0
        s%tsub(:, 0) = 0
        s%tdead(0) = .false.
        s%tinf(0) = .false.
        cap = max(2*npoints, 64)
        allocate (s%snbr(0:1, 0:cap), s%sver(0:1, 0:cap), s%ssegv(0:1, 0:cap))
        allocate (s%stri(0:1, 0:cap), s%smark(0:cap), s%sdead(0:cap), s%sfree(cap))
        s%snbr(:, 0) = 0
        s%sver(:, 0) = 0
        s%ssegv(:, 0) = 0
        s%stri(:, 0) = 0
        s%smark(0) = 0
        s%sdead(0) = .false.
        allocate (s%jtri(64), s%jkind(64))
        allocate (s%vir(64))
        allocate (s%btq_front(0:4095), s%btq_tail(0:4095), s%btq_nextq(0:4095))
        s%btq_front = -1
        cap = 64
        allocate (s%bt_tri(cap), s%bt_org(cap), s%bt_dest(cap), s%bt_apex(cap))
        allocate (s%bt_next(cap), s%bt_key(cap), s%bt_free(cap))
        allocate (s%bs_sub(cap), s%bs_org(cap), s%bs_dest(cap))
        allocate (s%bs_dead(cap), s%bs_free(cap))
        s%rngseed = 1
        s%samples = 1
    end subroutine tc_init_state

    function tc_random(s, choices) result(r)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: choices
        integer :: r

        s%rngseed = mod(s%rngseed*1366_int64 + 150889_int64, 714025_int64)
        r = int(s%rngseed/(714025_int64/int(choices, int64) + 1_int64))
    end function tc_random

    pure function t_slot(ot) result(t)
        integer, intent(in) :: ot
        integer :: t
        t = shiftr(ot, 2)
    end function t_slot

    pure function t_orient(ot) result(o)
        integer, intent(in) :: ot
        integer :: o
        o = iand(ot, 3)
    end function t_orient

    pure function t_make(slot, orient) result(ot)
        integer, intent(in) :: slot, orient
        integer :: ot
        ot = ior(shiftl(slot, 2), orient)
    end function t_make

    pure function s_slot(os) result(t)
        integer, intent(in) :: os
        integer :: t
        t = shiftr(os, 1)
    end function s_slot

    pure function s_so(os) result(o)
        integer, intent(in) :: os
        integer :: o
        o = iand(os, 1)
    end function s_so

    pure function s_make(slot, so) result(os)
        integer, intent(in) :: slot, so
        integer :: os
        os = ior(shiftl(slot, 1), so)
    end function s_make

    pure function t_sym(s, ot) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: r
        r = s%tnbr(iand(ot, 3), shiftr(ot, 2))
    end function t_sym

    pure function t_lnext(ot) result(r)
        integer, intent(in) :: ot
        integer :: r
        r = ior(iand(ot, not(3)), P1(iand(ot, 3)))
    end function t_lnext

    pure function t_lprev(ot) result(r)
        integer, intent(in) :: ot
        integer :: r
        r = ior(iand(ot, not(3)), M1(iand(ot, 3)))
    end function t_lprev

    pure function t_onext(s, ot) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: r
        r = t_sym(s, t_lprev(ot))
    end function t_onext

    pure function t_oprev(s, ot) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: r
        r = t_lnext(t_sym(s, ot))
    end function t_oprev

    pure function t_dnext(s, ot) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: r
        r = t_lprev(t_sym(s, ot))
    end function t_dnext

    pure function t_dprev(s, ot) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: r
        r = t_sym(s, t_lnext(ot))
    end function t_dprev

    pure function t_org(s, ot) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: v
        v = s%tver(P1(iand(ot, 3)), shiftr(ot, 2))
    end function t_org

    pure function t_dest(s, ot) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: v
        v = s%tver(M1(iand(ot, 3)), shiftr(ot, 2))
    end function t_dest

    pure function t_apex(s, ot) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: v
        v = s%tver(iand(ot, 3), shiftr(ot, 2))
    end function t_apex

    subroutine set_org(s, ot, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot, v
        s%tver(P1(iand(ot, 3)), shiftr(ot, 2)) = v
    end subroutine set_org

    subroutine set_dest(s, ot, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot, v
        s%tver(M1(iand(ot, 3)), shiftr(ot, 2)) = v
    end subroutine set_dest

    subroutine set_apex(s, ot, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot, v
        s%tver(iand(ot, 3), shiftr(ot, 2)) = v
    end subroutine set_apex

    subroutine t_bond(s, ot1, ot2)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot1, ot2
        s%tnbr(iand(ot1, 3), shiftr(ot1, 2)) = ot2
        s%tnbr(iand(ot2, 3), shiftr(ot2, 2)) = ot1
    end subroutine t_bond

    subroutine t_dissolve(s, ot)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot
        s%tnbr(iand(ot, 3), shiftr(ot, 2)) = 0
    end subroutine t_dissolve

    pure function ts_pivot(s, ot) result(os)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: ot
        integer :: os
        os = s%tsub(iand(ot, 3), shiftr(ot, 2))
    end function ts_pivot

    subroutine ts_bond(s, ot, os)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot, os
        s%tsub(iand(ot, 3), shiftr(ot, 2)) = os
        s%stri(iand(os, 1), shiftr(os, 1)) = ot
    end subroutine ts_bond

    subroutine ts_dissolve(s, ot)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot
        s%tsub(iand(ot, 3), shiftr(ot, 2)) = 0
    end subroutine ts_dissolve

    pure function st_pivot(s, os) result(ot)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: ot
        ot = s%stri(iand(os, 1), shiftr(os, 1))
    end function st_pivot

    subroutine st_dissolve(s, os)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os
        s%stri(iand(os, 1), shiftr(os, 1)) = 0
    end subroutine st_dissolve

    pure function s_sym(os) result(r)
        integer, intent(in) :: os
        integer :: r
        r = ieor(os, 1)
    end function s_sym

    pure function s_pivot(s, os) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: r
        r = s%snbr(iand(os, 1), shiftr(os, 1))
    end function s_pivot

    pure function s_next(s, os) result(r)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: r
        r = s%snbr(1 - iand(os, 1), shiftr(os, 1))
    end function s_next

    pure function s_org(s, os) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: v
        v = s%sver(iand(os, 1), shiftr(os, 1))
    end function s_org

    pure function s_dest(s, os) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: v
        v = s%sver(1 - iand(os, 1), shiftr(os, 1))
    end function s_dest

    subroutine set_s_org(s, os, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os, v
        s%sver(iand(os, 1), shiftr(os, 1)) = v
    end subroutine set_s_org

    subroutine set_s_dest(s, os, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os, v
        s%sver(1 - iand(os, 1), shiftr(os, 1)) = v
    end subroutine set_s_dest

    pure function seg_org(s, os) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: v
        v = s%ssegv(iand(os, 1), shiftr(os, 1))
    end function seg_org

    pure function seg_dest(s, os) result(v)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: os
        integer :: v
        v = s%ssegv(1 - iand(os, 1), shiftr(os, 1))
    end function seg_dest

    subroutine set_seg_org(s, os, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os, v
        s%ssegv(iand(os, 1), shiftr(os, 1)) = v
    end subroutine set_seg_org

    subroutine set_seg_dest(s, os, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os, v
        s%ssegv(1 - iand(os, 1), shiftr(os, 1)) = v
    end subroutine set_seg_dest

    subroutine s_bond(s, os1, os2)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os1, os2
        s%snbr(iand(os1, 1), shiftr(os1, 1)) = os2
        s%snbr(iand(os2, 1), shiftr(os2, 1)) = os1
    end subroutine s_bond

    subroutine s_dissolve(s, os)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: os
        s%snbr(iand(os, 1), shiftr(os, 1)) = 0
    end subroutine s_dissolve

end module tc_state
