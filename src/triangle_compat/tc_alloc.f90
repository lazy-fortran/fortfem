module tc_alloc
    !! Slot allocation and deallocation for vertices, triangles, and
    !! subsegments of the Triangle-compatible mesh state.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state, only: tc_state_t, t_make, s_make, V_INPUT, V_DEAD
    implicit none
    private

    public :: make_triangle, make_subseg, triangle_dealloc, subseg_dealloc
    public :: vertex_alloc, vertex_dealloc

contains

    function make_triangle(s) result(ot)
        type(tc_state_t), intent(inout) :: s
        integer :: ot
        integer :: slot, oldcap
        integer, allocatable :: itmp(:, :), i1(:)
        logical, allocatable :: ltmp(:)

        if (s%nt_free > 0) then
            slot = s%tfree(s%nt_free)
            s%nt_free = s%nt_free - 1
        else
            slot = s%nt_slots + 1
            oldcap = ubound(s%tdead, 1)
            if (slot > oldcap) then
                allocate (itmp(0:2, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%tnbr
                call move_alloc(itmp, s%tnbr)
                allocate (itmp(0:2, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%tver
                call move_alloc(itmp, s%tver)
                allocate (itmp(0:2, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%tsub
                call move_alloc(itmp, s%tsub)
                allocate (ltmp(0:2*oldcap))
                ltmp(0:oldcap) = s%tdead
                call move_alloc(ltmp, s%tdead)
                allocate (ltmp(0:2*oldcap))
                ltmp(0:oldcap) = s%tinf
                call move_alloc(ltmp, s%tinf)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%tfree
                call move_alloc(i1, s%tfree)
            end if
            s%nt_slots = slot
        end if
        s%nt_live = s%nt_live + 1
        s%tnbr(:, slot) = 0
        s%tver(:, slot) = 0
        s%tsub(:, slot) = 0
        s%tdead(slot) = .false.
        s%tinf(slot) = .false.
        ot = t_make(slot, 0)
    end function make_triangle

    subroutine triangle_dealloc(s, slot)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: slot
        s%tdead(slot) = .true.
        s%tver(:, slot) = 0
        s%nt_live = s%nt_live - 1
        s%nt_free = s%nt_free + 1
        s%tfree(s%nt_free) = slot
    end subroutine triangle_dealloc

    function make_subseg(s) result(os)
        type(tc_state_t), intent(inout) :: s
        integer :: os
        integer :: slot, oldcap
        integer, allocatable :: itmp(:, :), i1(:)
        logical, allocatable :: ltmp(:)

        if (s%ns_free > 0) then
            slot = s%sfree(s%ns_free)
            s%ns_free = s%ns_free - 1
        else
            slot = s%ns_slots + 1
            oldcap = ubound(s%sdead, 1)
            if (slot > oldcap) then
                allocate (itmp(0:1, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%snbr
                call move_alloc(itmp, s%snbr)
                allocate (itmp(0:1, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%sver
                call move_alloc(itmp, s%sver)
                allocate (itmp(0:1, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%ssegv
                call move_alloc(itmp, s%ssegv)
                allocate (itmp(0:1, 0:2*oldcap))
                itmp(:, 0:oldcap) = s%stri
                call move_alloc(itmp, s%stri)
                allocate (i1(0:2*oldcap))
                i1(0:oldcap) = s%smark
                call move_alloc(i1, s%smark)
                allocate (ltmp(0:2*oldcap))
                ltmp(0:oldcap) = s%sdead
                call move_alloc(ltmp, s%sdead)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%sfree
                call move_alloc(i1, s%sfree)
            end if
            s%ns_slots = slot
        end if
        s%ns_live = s%ns_live + 1
        s%snbr(:, slot) = 0
        s%sver(:, slot) = 0
        s%ssegv(:, slot) = 0
        s%stri(:, slot) = 0
        s%smark(slot) = 0
        s%sdead(slot) = .false.
        os = s_make(slot, 0)
    end function make_subseg

    subroutine subseg_dealloc(s, slot)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: slot
        s%sdead(slot) = .true.
        s%sver(:, slot) = 0
        s%ns_live = s%ns_live - 1
        s%ns_free = s%ns_free + 1
        s%sfree(s%ns_free) = slot
    end subroutine subseg_dealloc

    function vertex_alloc(s) result(v)
        type(tc_state_t), intent(inout) :: s
        integer :: v
        integer :: oldcap
        integer, allocatable :: i1(:)
        real(dp), allocatable :: r1(:)

        if (s%nv_free > 0) then
            v = s%vfree(s%nv_free)
            s%nv_free = s%nv_free - 1
        else
            v = s%nv_slots + 1
            oldcap = size(s%vx)
            if (v > oldcap) then
                allocate (r1(2*oldcap))
                r1(1:oldcap) = s%vx
                call move_alloc(r1, s%vx)
                allocate (r1(2*oldcap))
                r1(1:oldcap) = s%vy
                call move_alloc(r1, s%vy)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%vmark
                call move_alloc(i1, s%vmark)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%vtype
                call move_alloc(i1, s%vtype)
                allocate (i1(2*oldcap))
                i1 = -1
                i1(1:oldcap) = s%vtri
                call move_alloc(i1, s%vtri)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%vfree
                call move_alloc(i1, s%vfree)
            end if
            s%nv_slots = v
        end if
        s%nv_live = s%nv_live + 1
        s%vtri(v) = -1
        s%vmark(v) = 0
        s%vtype(v) = V_INPUT
    end function vertex_alloc

    subroutine vertex_dealloc(s, v)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: v
        s%vtype(v) = V_DEAD
        s%nv_live = s%nv_live - 1
        s%nv_free = s%nv_free + 1
        s%vfree(s%nv_free) = v
    end subroutine vertex_dealloc

end module tc_alloc
