module tc_segment
    !! PSLG segment recovery in the constrained Delaunay triangulation.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state
    use tc_alloc
    use tc_geom, only: ccw_v
    use tc_locate, only: locate, finddirection, makevertexmap, LOC_ONVERTEX, &
                         DIR_WITHIN, DIR_LEFTCOLLINEAR, DIR_RIGHTCOLLINEAR
    use tc_flip, only: flip
    use tc_insert, only: insertsubseg, insertvertex, IV_SUCCESS
    use tc_geom, only: incircle_v
    implicit none
    private

    public :: formskeleton

contains

    subroutine formskeleton(s, seglist, segmarkers)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: seglist(:, :)
        integer, intent(in) :: segmarkers(:)
        integer :: i, end1, end2, boundmarker

        s%insegments = size(seglist, 2)
        if (s%nt_live == 0) return
        if (s%insegments > 0) call makevertexmap(s)
        boundmarker = 0
        do i = 1, s%insegments
            end1 = seglist(1, i)
            end2 = seglist(2, i)
            if (size(segmarkers) > 0) boundmarker = segmarkers(i)
            if ((end1 < 1) .or. (end1 > s%invertices)) cycle
            if ((end2 < 1) .or. (end2 > s%invertices)) cycle
            if ((s%vx(end1) == s%vx(end2)) .and. &
                (s%vy(end1) == s%vy(end2))) cycle
            call insertsegment(s, end1, end2, boundmarker)
        end do
    end subroutine formskeleton

    subroutine insertsegment(s, endpoint1_in, endpoint2_in, newmark)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: endpoint1_in, endpoint2_in, newmark
        integer :: endpoint1, endpoint2
        integer :: searchtri1, searchtri2, checkvertex
        integer :: code

        endpoint1 = endpoint1_in
        endpoint2 = endpoint2_in

        checkvertex = 0
        searchtri1 = s%vtri(endpoint1)
        if (searchtri1 >= 0) checkvertex = t_org(s, searchtri1)
        if (checkvertex /= endpoint1) then
            searchtri1 = s%tnbr(0, 0)
            code = locate(s, s%vx(endpoint1), s%vy(endpoint1), searchtri1)
            if (code /= LOC_ONVERTEX) then
                error stop 'tc_segment: unable to locate PSLG vertex'
            end if
        end if
        s%recent = searchtri1
        if (scoutsegment(s, searchtri1, endpoint2, newmark)) return
        endpoint1 = t_org(s, searchtri1)

        checkvertex = 0
        searchtri2 = s%vtri(endpoint2)
        if (searchtri2 >= 0) checkvertex = t_org(s, searchtri2)
        if (checkvertex /= endpoint2) then
            searchtri2 = s%tnbr(0, 0)
            code = locate(s, s%vx(endpoint2), s%vy(endpoint2), searchtri2)
            if (code /= LOC_ONVERTEX) then
                error stop 'tc_segment: unable to locate PSLG vertex'
            end if
        end if
        s%recent = searchtri2
        if (scoutsegment(s, searchtri2, endpoint1, newmark)) return
        endpoint2 = t_org(s, searchtri2)

        call constrainededge(s, searchtri1, endpoint2, newmark)
    end subroutine insertsegment

    recursive function scoutsegment(s, searchtri, endpoint2, newmark) &
        result(done)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: searchtri
        integer, intent(in) :: endpoint2, newmark
        logical :: done
        integer :: crosstri, crosssubseg, leftvertex, rightvertex, collinear

        collinear = finddirection(s, searchtri, endpoint2)
        rightvertex = t_dest(s, searchtri)
        leftvertex = t_apex(s, searchtri)
        if (((s%vx(leftvertex) == s%vx(endpoint2)) .and. &
             (s%vy(leftvertex) == s%vy(endpoint2))) .or. &
            ((s%vx(rightvertex) == s%vx(endpoint2)) .and. &
             (s%vy(rightvertex) == s%vy(endpoint2)))) then
            if ((s%vx(leftvertex) == s%vx(endpoint2)) .and. &
                (s%vy(leftvertex) == s%vy(endpoint2))) then
                searchtri = t_lprev(searchtri)
            end if
            call insertsubseg(s, searchtri, newmark)
            done = .true.
            return
        else if (collinear == DIR_LEFTCOLLINEAR) then
            searchtri = t_lprev(searchtri)
            call insertsubseg(s, searchtri, newmark)
            done = scoutsegment(s, searchtri, endpoint2, newmark)
            return
        else if (collinear == DIR_RIGHTCOLLINEAR) then
            call insertsubseg(s, searchtri, newmark)
            searchtri = t_lnext(searchtri)
            done = scoutsegment(s, searchtri, endpoint2, newmark)
            return
        else
            crosstri = t_lnext(searchtri)
            crosssubseg = ts_pivot(s, crosstri)
            if (s_slot(crosssubseg) == 0) then
                done = .false.
                return
            else
                call segmentintersection(s, crosstri, crosssubseg, endpoint2)
                searchtri = crosstri
                call insertsubseg(s, searchtri, newmark)
                done = scoutsegment(s, searchtri, endpoint2, newmark)
                return
            end if
        end if
    end function scoutsegment

    subroutine segmentintersection(s, splittri, splitsubseg, endpoint2)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: splittri
        integer, intent(in) :: splitsubseg
        integer, intent(in) :: endpoint2
        integer :: opposubseg, endpoint1, torg, tdest
        integer :: leftvertex, rightvertex, newvertex, success, collinear
        integer :: splitsub, cursub
        real(dp) :: ex, ey, tx, ty, etx, ety, split, denom

        endpoint1 = t_apex(s, splittri)
        torg = t_org(s, splittri)
        tdest = t_dest(s, splittri)
        tx = s%vx(tdest) - s%vx(torg)
        ty = s%vy(tdest) - s%vy(torg)
        ex = s%vx(endpoint2) - s%vx(endpoint1)
        ey = s%vy(endpoint2) - s%vy(endpoint1)
        etx = s%vx(torg) - s%vx(endpoint2)
        ety = s%vy(torg) - s%vy(endpoint2)
        denom = ty*ex - tx*ey
        if (denom == 0.0_dp) then
            error stop 'tc_segment: parallel segment intersection'
        end if
        split = (ey*etx - ex*ety)/denom
        newvertex = vertex_alloc(s)
        s%vx(newvertex) = s%vx(torg) + split*(s%vx(tdest) - s%vx(torg))
        s%vy(newvertex) = s%vy(torg) + split*(s%vy(tdest) - s%vy(torg))
        s%vmark(newvertex) = s%smark(s_slot(splitsubseg))
        s%vtype(newvertex) = V_INPUT
        success = insertvertex(s, newvertex, splittri, splitsubseg, &
                               .false., .false.)
        if (success /= IV_SUCCESS) then
            error stop 'tc_segment: failure to split a segment'
        end if
        s%vtri(newvertex) = splittri
        if (s%steinerleft > 0) s%steinerleft = s%steinerleft - 1

        splitsub = s_sym(splitsubseg)
        opposubseg = s_pivot(s, splitsub)
        call s_dissolve(s, splitsub)
        call s_dissolve(s, opposubseg)
        cursub = splitsub
        do
            call set_seg_org(s, cursub, newvertex)
            cursub = s_next(s, cursub)
            if (s_slot(cursub) == 0) exit
        end do
        cursub = opposubseg
        do
            call set_seg_org(s, cursub, newvertex)
            cursub = s_next(s, cursub)
            if (s_slot(cursub) == 0) exit
        end do

        collinear = finddirection(s, splittri, endpoint1)
        rightvertex = t_dest(s, splittri)
        leftvertex = t_apex(s, splittri)
        if ((s%vx(leftvertex) == s%vx(endpoint1)) .and. &
            (s%vy(leftvertex) == s%vy(endpoint1))) then
            splittri = t_onext(s, splittri)
        else if ((s%vx(rightvertex) /= s%vx(endpoint1)) .or. &
                 (s%vy(rightvertex) /= s%vy(endpoint1))) then
            error stop 'tc_segment: topological inconsistency after split'
        end if
    end subroutine segmentintersection

    recursive subroutine delaunayfixup(s, fixuptri, leftside)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: fixuptri
        logical, intent(in) :: leftside
        integer :: neartri, fartri, faredge
        integer :: nearvertex, leftvertex, rightvertex, farvertex

        neartri = t_lnext(fixuptri)
        fartri = t_sym(s, neartri)
        if (t_slot(fartri) == 0) return
        faredge = ts_pivot(s, neartri)
        if (s_slot(faredge) /= 0) return
        nearvertex = t_apex(s, neartri)
        leftvertex = t_org(s, neartri)
        rightvertex = t_dest(s, neartri)
        farvertex = t_apex(s, fartri)
        if (leftside) then
            if (ccw_v(s, nearvertex, leftvertex, farvertex) <= 0.0_dp) return
        else
            if (ccw_v(s, farvertex, rightvertex, nearvertex) <= 0.0_dp) return
        end if
        if (ccw_v(s, rightvertex, leftvertex, farvertex) > 0.0_dp) then
            if (incircle_v(s, leftvertex, farvertex, rightvertex, nearvertex) &
                <= 0.0_dp) return
        end if
        call flip(s, neartri)
        fixuptri = t_lprev(fixuptri)
        call delaunayfixup(s, fixuptri, leftside)
        call delaunayfixup(s, fartri, leftside)
    end subroutine delaunayfixup

    recursive subroutine constrainededge(s, starttri, endpoint2, newmark)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: starttri
        integer, intent(in) :: endpoint2, newmark
        integer :: fixuptri, fixuptri2, crosssubseg
        integer :: endpoint1, farvertex
        real(dp) :: area
        logical :: collision, done

        endpoint1 = t_org(s, starttri)
        fixuptri = t_lnext(starttri)
        call flip(s, fixuptri)
        collision = .false.
        done = .false.
        do
            farvertex = t_org(s, fixuptri)
            if ((s%vx(farvertex) == s%vx(endpoint2)) .and. &
                (s%vy(farvertex) == s%vy(endpoint2))) then
                fixuptri2 = t_oprev(s, fixuptri)
                call delaunayfixup(s, fixuptri, .false.)
                call delaunayfixup(s, fixuptri2, .true.)
                done = .true.
            else
                area = ccw_v(s, endpoint1, endpoint2, farvertex)
                if (area == 0.0_dp) then
                    collision = .true.
                    fixuptri2 = t_oprev(s, fixuptri)
                    call delaunayfixup(s, fixuptri, .false.)
                    call delaunayfixup(s, fixuptri2, .true.)
                    done = .true.
                else
                    if (area > 0.0_dp) then
                        fixuptri2 = t_oprev(s, fixuptri)
                        call delaunayfixup(s, fixuptri2, .true.)
                        fixuptri = t_lprev(fixuptri)
                    else
                        call delaunayfixup(s, fixuptri, .false.)
                        fixuptri = t_oprev(s, fixuptri)
                    end if
                    crosssubseg = ts_pivot(s, fixuptri)
                    if (s_slot(crosssubseg) == 0) then
                        call flip(s, fixuptri)
                    else
                        collision = .true.
                        call segmentintersection(s, fixuptri, crosssubseg, &
                                                 endpoint2)
                        done = .true.
                    end if
                end if
            end if
            if (done) exit
        end do
        call insertsubseg(s, fixuptri, newmark)
        if (collision) then
            if (.not. scoutsegment(s, fixuptri, endpoint2, newmark)) then
                call constrainededge(s, fixuptri, endpoint2, newmark)
            end if
        end if
    end subroutine constrainededge

end module tc_segment
