module tc_insert
    !! Vertex insertion with incremental flipping, subsegment insertion,
    !! and the flip journal that allows undoing a rejected insertion.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state
    use tc_alloc
    use tc_flip, only: unflip
    use tc_geom, only: ccw_v, incircle_v
    use tc_locate, only: preciselocate, locate, LOC_INTRIANGLE, LOC_ONEDGE, &
                         LOC_ONVERTEX, LOC_OUTSIDE
    use tc_qtest, only: checkseg4encroach, testtriangle, badsubseg_alloc
    implicit none
    private

    public :: insertsubseg, insertvertex, undovertex
    public :: IV_SUCCESS, IV_ENCROACHING, IV_VIOLATING, IV_DUPLICATE
    public :: JK_TRISECT, JK_BISECT, JK_FLIP

    integer, parameter :: IV_SUCCESS = 1, IV_ENCROACHING = 2
    integer, parameter :: IV_VIOLATING = 3, IV_DUPLICATE = 4
    integer, parameter :: JK_TRISECT = 1, JK_BISECT = 2, JK_FLIP = 3

contains

    subroutine journal_push(s, ot, kind)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: ot, kind
        integer, allocatable :: i1(:)
        integer :: oldcap

        oldcap = size(s%jtri)
        if (s%jn + 1 > oldcap) then
            allocate (i1(2*oldcap))
            i1(1:oldcap) = s%jtri
            call move_alloc(i1, s%jtri)
            allocate (i1(2*oldcap))
            i1(1:oldcap) = s%jkind
            call move_alloc(i1, s%jkind)
        end if
        s%jn = s%jn + 1
        s%jtri(s%jn) = ot
        s%jkind(s%jn) = kind
    end subroutine journal_push

    subroutine insertsubseg(s, tri, subsegmark)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: tri, subsegmark
        integer :: oppotri, newsubseg, triorg, tridest

        triorg = t_org(s, tri)
        tridest = t_dest(s, tri)
        if (s%vmark(triorg) == 0) s%vmark(triorg) = subsegmark
        if (s%vmark(tridest) == 0) s%vmark(tridest) = subsegmark
        newsubseg = ts_pivot(s, tri)
        if (s_slot(newsubseg) == 0) then
            newsubseg = make_subseg(s)
            call set_s_org(s, newsubseg, tridest)
            call set_s_dest(s, newsubseg, triorg)
            call set_seg_org(s, newsubseg, tridest)
            call set_seg_dest(s, newsubseg, triorg)
            call ts_bond(s, tri, newsubseg)
            oppotri = t_sym(s, tri)
            newsubseg = s_sym(newsubseg)
            call ts_bond(s, oppotri, newsubseg)
            s%smark(s_slot(newsubseg)) = subsegmark
        else
            if (s%smark(s_slot(newsubseg)) == 0) then
                s%smark(s_slot(newsubseg)) = subsegmark
            end if
        end if
    end subroutine insertsubseg

    function insertvertex(s, newvertex, searchtri, splitseg, segmentflaws, &
                          triflaws) result(success)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: newvertex
        integer, intent(inout) :: searchtri
        integer, intent(in) :: splitseg   ! encoded osub, or -1
        logical, intent(in) :: segmentflaws, triflaws
        integer :: success
        integer :: horiz, top, botleft, botright, topleft, topright
        integer :: newbotleft, newbotright, newtopright
        integer :: botlcasing, botrcasing, toplcasing, toprcasing
        integer :: botlsubseg, botrsubseg, toplsubseg, toprsubseg
        integer :: brokensubseg, checksubseg, rightsubseg, newsubseg
        integer :: testtri
        integer :: first, leftvertex, rightvertex, botvertex, topvertex
        integer :: farvertex, segmentorg, segmentdest
        integer :: intersect, idx, splitseg_l
        logical :: mirrorflag, doflip, enq

        splitseg_l = splitseg
        if (splitseg_l < 0) then
            if (t_slot(searchtri) == 0) then
                horiz = s%tnbr(0, 0)
                intersect = locate(s, s%vx(newvertex), s%vy(newvertex), horiz)
            else
                horiz = searchtri
                intersect = preciselocate(s, s%vx(newvertex), s%vy(newvertex), &
                                          horiz, .true.)
            end if
        else
            horiz = searchtri
            intersect = LOC_ONEDGE
        end if

        if (intersect == LOC_ONVERTEX) then
            searchtri = horiz
            s%recent = horiz
            success = IV_DUPLICATE
            return
        end if
        if ((intersect == LOC_ONEDGE) .or. (intersect == LOC_OUTSIDE)) then
            if (s%checksegments .and. (splitseg_l < 0)) then
                brokensubseg = ts_pivot(s, horiz)
                if (s_slot(brokensubseg) /= 0) then
                    if (segmentflaws) then
                        enq = s%nobisect /= 2
                        if (enq .and. (s%nobisect == 1)) then
                            testtri = t_sym(s, horiz)
                            enq = t_slot(testtri) /= 0
                        end if
                        if (enq) then
                            idx = badsubseg_alloc(s)
                            s%bs_sub(idx) = brokensubseg
                            s%bs_org(idx) = s_org(s, brokensubseg)
                            s%bs_dest(idx) = s_dest(s, brokensubseg)
                        end if
                    end if
                    searchtri = horiz
                    s%recent = horiz
                    success = IV_VIOLATING
                    return
                end if
            end if

            botright = t_lprev(horiz)
            botrcasing = t_sym(s, botright)
            topright = t_sym(s, horiz)
            mirrorflag = t_slot(topright) /= 0
            if (mirrorflag) then
                topright = t_lnext(topright)
                toprcasing = t_sym(s, topright)
                newtopright = make_triangle(s)
            else
                s%hullsize = s%hullsize + 1
            end if
            newbotright = make_triangle(s)

            rightvertex = t_org(s, horiz)
            leftvertex = t_dest(s, horiz)
            botvertex = t_apex(s, horiz)
            call set_org(s, newbotright, botvertex)
            call set_dest(s, newbotright, rightvertex)
            call set_apex(s, newbotright, newvertex)
            call set_org(s, horiz, newvertex)
            if (mirrorflag) then
                topvertex = t_dest(s, topright)
                call set_org(s, newtopright, rightvertex)
                call set_dest(s, newtopright, topvertex)
                call set_apex(s, newtopright, newvertex)
                call set_org(s, topright, newvertex)
            end if

            if (s%checksegments) then
                botrsubseg = ts_pivot(s, botright)
                if (s_slot(botrsubseg) /= 0) then
                    call ts_dissolve(s, botright)
                    call ts_bond(s, newbotright, botrsubseg)
                end if
                if (mirrorflag) then
                    toprsubseg = ts_pivot(s, topright)
                    if (s_slot(toprsubseg) /= 0) then
                        call ts_dissolve(s, topright)
                        call ts_bond(s, newtopright, toprsubseg)
                    end if
                end if
            end if

            call t_bond(s, newbotright, botrcasing)
            newbotright = t_lprev(newbotright)
            call t_bond(s, newbotright, botright)
            newbotright = t_lprev(newbotright)
            if (mirrorflag) then
                call t_bond(s, newtopright, toprcasing)
                newtopright = t_lnext(newtopright)
                call t_bond(s, newtopright, topright)
                newtopright = t_lnext(newtopright)
                call t_bond(s, newtopright, newbotright)
            end if

            if (splitseg_l >= 0) then
                call set_s_dest(s, splitseg_l, newvertex)
                segmentorg = seg_org(s, splitseg_l)
                segmentdest = seg_dest(s, splitseg_l)
                splitseg_l = s_sym(splitseg_l)
                rightsubseg = s_pivot(s, splitseg_l)
                call insertsubseg(s, newbotright, s%smark(s_slot(splitseg_l)))
                newsubseg = ts_pivot(s, newbotright)
                call set_seg_org(s, newsubseg, segmentorg)
                call set_seg_dest(s, newsubseg, segmentdest)
                call s_bond(s, splitseg_l, newsubseg)
                newsubseg = s_sym(newsubseg)
                call s_bond(s, newsubseg, rightsubseg)
                splitseg_l = s_sym(splitseg_l)
                if (s%vmark(newvertex) == 0) then
                    s%vmark(newvertex) = s%smark(s_slot(splitseg_l))
                end if
            end if

            if (s%checkquality) then
                s%jn = 0
                call journal_push(s, horiz, JK_BISECT)
            end if

            horiz = t_lnext(horiz)
        else
            botleft = t_lnext(horiz)
            botright = t_lprev(horiz)
            botlcasing = t_sym(s, botleft)
            botrcasing = t_sym(s, botright)
            newbotleft = make_triangle(s)
            newbotright = make_triangle(s)

            rightvertex = t_org(s, horiz)
            leftvertex = t_dest(s, horiz)
            botvertex = t_apex(s, horiz)
            call set_org(s, newbotleft, leftvertex)
            call set_dest(s, newbotleft, botvertex)
            call set_apex(s, newbotleft, newvertex)
            call set_org(s, newbotright, botvertex)
            call set_dest(s, newbotright, rightvertex)
            call set_apex(s, newbotright, newvertex)
            call set_apex(s, horiz, newvertex)

            if (s%checksegments) then
                botlsubseg = ts_pivot(s, botleft)
                if (s_slot(botlsubseg) /= 0) then
                    call ts_dissolve(s, botleft)
                    call ts_bond(s, newbotleft, botlsubseg)
                end if
                botrsubseg = ts_pivot(s, botright)
                if (s_slot(botrsubseg) /= 0) then
                    call ts_dissolve(s, botright)
                    call ts_bond(s, newbotright, botrsubseg)
                end if
            end if

            call t_bond(s, newbotleft, botlcasing)
            call t_bond(s, newbotright, botrcasing)
            newbotleft = t_lnext(newbotleft)
            newbotright = t_lprev(newbotright)
            call t_bond(s, newbotleft, newbotright)
            newbotleft = t_lnext(newbotleft)
            call t_bond(s, botleft, newbotleft)
            newbotright = t_lprev(newbotright)
            call t_bond(s, botright, newbotright)

            if (s%checkquality) then
                s%jn = 0
                call journal_push(s, horiz, JK_TRISECT)
            end if
        end if

        success = IV_SUCCESS
        first = t_org(s, horiz)
        rightvertex = first
        leftvertex = t_dest(s, horiz)
        do
            doflip = .true.
            if (s%checksegments) then
                checksubseg = ts_pivot(s, horiz)
                if (s_slot(checksubseg) /= 0) then
                    doflip = .false.
                    if (segmentflaws) then
                        if (checkseg4encroach(s, checksubseg) /= 0) then
                            success = IV_ENCROACHING
                        end if
                    end if
                end if
            end if

            if (doflip) then
                top = t_sym(s, horiz)
                if (t_slot(top) == 0) then
                    doflip = .false.
                else
                    farvertex = t_apex(s, top)
                    doflip = incircle_v(s, leftvertex, newvertex, rightvertex, &
                                        farvertex) > 0.0_dp
                    if (doflip) then
                        topleft = t_lprev(top)
                        toplcasing = t_sym(s, topleft)
                        topright = t_lnext(top)
                        toprcasing = t_sym(s, topright)
                        botleft = t_lnext(horiz)
                        botlcasing = t_sym(s, botleft)
                        botright = t_lprev(horiz)
                        botrcasing = t_sym(s, botright)
                        call t_bond(s, topleft, botlcasing)
                        call t_bond(s, botleft, botrcasing)
                        call t_bond(s, botright, toprcasing)
                        call t_bond(s, topright, toplcasing)
                        if (s%checksegments) then
                            toplsubseg = ts_pivot(s, topleft)
                            botlsubseg = ts_pivot(s, botleft)
                            botrsubseg = ts_pivot(s, botright)
                            toprsubseg = ts_pivot(s, topright)
                            if (s_slot(toplsubseg) == 0) then
                                call ts_dissolve(s, topright)
                            else
                                call ts_bond(s, topright, toplsubseg)
                            end if
                            if (s_slot(botlsubseg) == 0) then
                                call ts_dissolve(s, topleft)
                            else
                                call ts_bond(s, topleft, botlsubseg)
                            end if
                            if (s_slot(botrsubseg) == 0) then
                                call ts_dissolve(s, botleft)
                            else
                                call ts_bond(s, botleft, botrsubseg)
                            end if
                            if (s_slot(toprsubseg) == 0) then
                                call ts_dissolve(s, botright)
                            else
                                call ts_bond(s, botright, toprsubseg)
                            end if
                        end if
                        call set_org(s, horiz, farvertex)
                        call set_dest(s, horiz, newvertex)
                        call set_apex(s, horiz, rightvertex)
                        call set_org(s, top, newvertex)
                        call set_dest(s, top, farvertex)
                        call set_apex(s, top, leftvertex)
                        if (s%checkquality) then
                            call journal_push(s, horiz, JK_FLIP)
                        end if
                        horiz = t_lprev(horiz)
                        leftvertex = farvertex
                    end if
                end if
            end if
            if (.not. doflip) then
                if (triflaws) call testtriangle(s, horiz)
                horiz = t_lnext(horiz)
                testtri = t_sym(s, horiz)
                if ((leftvertex == first) .or. (t_slot(testtri) == 0)) then
                    searchtri = t_lnext(horiz)
                    s%recent = t_lnext(horiz)
                    return
                end if
                horiz = t_lnext(testtri)
                rightvertex = leftvertex
                leftvertex = t_dest(s, horiz)
            end if
        end do
    end function insertvertex

    subroutine undovertex(s)
        type(tc_state_t), intent(inout) :: s
        integer :: fliptri, botleft, botright, topright, gluetri
        integer :: botlcasing, botrcasing, toprcasing
        integer :: botlsubseg, botrsubseg, toprsubseg
        integer :: botvertex, rightvertex
        integer :: k

        do k = s%jn, 1, -1
            fliptri = s%jtri(k)
            select case (s%jkind(k))
            case (JK_TRISECT)
                botleft = t_lnext(t_dprev(s, fliptri))
                botright = t_lprev(t_onext(s, fliptri))
                botlcasing = t_sym(s, botleft)
                botrcasing = t_sym(s, botright)
                botvertex = t_dest(s, botleft)

                call set_apex(s, fliptri, botvertex)
                fliptri = t_lnext(fliptri)
                call t_bond(s, fliptri, botlcasing)
                botlsubseg = ts_pivot(s, botleft)
                call ts_bond(s, fliptri, botlsubseg)
                fliptri = t_lnext(fliptri)
                call t_bond(s, fliptri, botrcasing)
                botrsubseg = ts_pivot(s, botright)
                call ts_bond(s, fliptri, botrsubseg)

                call triangle_dealloc(s, t_slot(botleft))
                call triangle_dealloc(s, t_slot(botright))
            case (JK_BISECT)
                gluetri = t_lprev(fliptri)
                botright = t_lnext(t_sym(s, gluetri))
                botrcasing = t_sym(s, botright)
                rightvertex = t_dest(s, botright)

                call set_org(s, fliptri, rightvertex)
                call t_bond(s, gluetri, botrcasing)
                botrsubseg = ts_pivot(s, botright)
                call ts_bond(s, gluetri, botrsubseg)

                call triangle_dealloc(s, t_slot(botright))

                gluetri = t_sym(s, fliptri)
                if (t_slot(gluetri) /= 0) then
                    gluetri = t_lnext(gluetri)
                    topright = t_dnext(s, gluetri)
                    toprcasing = t_sym(s, topright)

                    call set_org(s, gluetri, rightvertex)
                    call t_bond(s, gluetri, toprcasing)
                    toprsubseg = ts_pivot(s, topright)
                    call ts_bond(s, gluetri, toprsubseg)

                    call triangle_dealloc(s, t_slot(topright))
                end if
            case (JK_FLIP)
                call unflip(s, fliptri)
            end select
        end do
        s%jn = 0
    end subroutine undovertex

end module tc_insert
