module tc_enforce
    !! Ruppert/Chew quality enforcement: encroached-subsegment splitting and
    !! bad-triangle splitting at circumcenters or off-centers, with -Y
    !! (no boundary Steiner points) semantics.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state
    use tc_alloc
    use tc_geom, only: ccw_v, incircle_v
    use tc_flip, only: flip
    use tc_insert, only: insertvertex, undovertex, &
                         IV_SUCCESS, IV_ENCROACHING, IV_VIOLATING, IV_DUPLICATE
    use tc_qtest, only: checkseg4encroach, testtriangle, findcircumcenter, &
                        enqueuebadtriang, dequeuebadtriang, badtri_release, &
                        badsubseg_dealloc, tallyencs, tallyfaces
    implicit none
    private

    public :: enforcequality

contains

    subroutine enforcequality(s)
        type(tc_state_t), intent(inout) :: s
        integer :: idx

        call tallyencs(s)
        call splitencsegs(s, .false.)

        if ((s%minangle > 0.0_dp) .or. s%fixedarea) then
            s%btq_front = -1
            s%btq_first = -1
            call tallyfaces(s)
            s%checkquality = .true.
            do while ((s%bt_items > 0) .and. (s%steinerleft /= 0))
                idx = dequeuebadtriang(s)
                call splittriangle(s, idx)
                if (s%bs_items > 0) then
                    call enqueuebadtriang(s, idx)
                    call splitencsegs(s, .true.)
                else
                    call badtri_release(s, idx)
                end if
            end do
        end if
    end subroutine enforcequality

    subroutine splittriangle(s, idx)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: idx
        integer :: badotri, borg, bdest, bapex, newvertex, success
        real(dp) :: ccx, ccy, xi, eta

        badotri = s%bt_tri(idx)
        borg = t_org(s, badotri)
        bdest = t_dest(s, badotri)
        bapex = t_apex(s, badotri)

        if (s%tdead(t_slot(badotri))) return
        if ((borg /= s%bt_org(idx)) .or. (bdest /= s%bt_dest(idx)) .or. &
            (bapex /= s%bt_apex(idx))) return

        newvertex = vertex_alloc(s)
        call findcircumcenter(s, borg, bdest, bapex, ccx, ccy, xi, eta, .true.)

        if (((ccx == s%vx(borg)) .and. (ccy == s%vy(borg))) .or. &
            ((ccx == s%vx(bdest)) .and. (ccy == s%vy(bdest))) .or. &
            ((ccx == s%vx(bapex)) .and. (ccy == s%vy(bapex)))) then
            call vertex_dealloc(s, newvertex)
        else
            s%vx(newvertex) = ccx
            s%vy(newvertex) = ccy
            s%vmark(newvertex) = 0
            s%vtype(newvertex) = V_FREE
            if (eta < xi) badotri = t_lprev(badotri)
            success = insertvertex(s, newvertex, badotri, -1, .true., .true.)
            if (success == IV_SUCCESS) then
                if (s%steinerleft > 0) s%steinerleft = s%steinerleft - 1
            else if (success == IV_ENCROACHING) then
                call undovertex(s)
                call vertex_dealloc(s, newvertex)
            else if (success == IV_VIOLATING) then
                call vertex_dealloc(s, newvertex)
            else
                call vertex_dealloc(s, newvertex)
            end if
        end if
    end subroutine splittriangle

    subroutine splitencsegs(s, triflaws)
        type(tc_state_t), intent(inout) :: s
        logical, intent(in) :: triflaws
        integer :: encidx, currentenc, enctri, testtri, testsh, eorg, edest
        integer :: eapex, newvertex, success, dummy_
        logical :: acuteorg, acutedest, acuteorg2, acutedest2
        real(dp) :: segmentlength, nearestpoweroftwo, split
        real(dp) :: multiplier, divisor

        do while ((s%bs_items > 0) .and. (s%steinerleft /= 0))
            encidx = 1
            do while (encidx <= s%bs_slots)
                if (s%steinerleft == 0) exit
                if (s%bs_dead(encidx)) then
                    encidx = encidx + 1
                    cycle
                end if
                currentenc = s%bs_sub(encidx)
                eorg = s_org(s, currentenc)
                edest = s_dest(s, currentenc)
                if ((.not. s%sdead(s_slot(currentenc))) .and. &
                    (eorg == s%bs_org(encidx)) .and. &
                    (edest == s%bs_dest(encidx))) then
                    enctri = st_pivot(s, currentenc)
                    testtri = t_lnext(enctri)
                    testsh = ts_pivot(s, testtri)
                    acuteorg = s_slot(testsh) /= 0
                    testtri = t_lnext(testtri)
                    testsh = ts_pivot(s, testtri)
                    acutedest = s_slot(testsh) /= 0

                    if ((.not. acuteorg) .and. (.not. acutedest)) then
                        eapex = t_apex(s, enctri)
                        do while ((s%vtype(eapex) == V_FREE) .and. &
                                  (chew_dot(s, eorg, edest, eapex) < 0.0_dp))
                            call deletevertex(s, testtri)
                            enctri = st_pivot(s, currentenc)
                            eapex = t_apex(s, enctri)
                            testtri = t_lprev(enctri)
                        end do
                    end if

                    testtri = t_sym(s, enctri)
                    if (t_slot(testtri) /= 0) then
                        testtri = t_lnext(testtri)
                        testsh = ts_pivot(s, testtri)
                        acutedest2 = s_slot(testsh) /= 0
                        acutedest = acutedest .or. acutedest2
                        testtri = t_lnext(testtri)
                        testsh = ts_pivot(s, testtri)
                        acuteorg2 = s_slot(testsh) /= 0
                        acuteorg = acuteorg .or. acuteorg2

                        if ((.not. acuteorg2) .and. (.not. acutedest2)) then
                            eapex = t_org(s, testtri)
                            do while ((s%vtype(eapex) == V_FREE) .and. &
                                      (chew_dot(s, eorg, edest, eapex) &
                                       < 0.0_dp))
                                call deletevertex(s, testtri)
                                testtri = t_sym(s, enctri)
                                eapex = t_apex(s, testtri)
                                testtri = t_lprev(testtri)
                            end do
                        end if
                    end if

                    if (acuteorg .or. acutedest) then
                        segmentlength = sqrt((s%vx(edest) - s%vx(eorg))* &
                                             (s%vx(edest) - s%vx(eorg)) + &
                                             (s%vy(edest) - s%vy(eorg))* &
                                             (s%vy(edest) - s%vy(eorg)))
                        nearestpoweroftwo = 1.0_dp
                        do while (segmentlength > 3.0_dp*nearestpoweroftwo)
                            nearestpoweroftwo = nearestpoweroftwo*2.0_dp
                        end do
                        do while (segmentlength < 1.5_dp*nearestpoweroftwo)
                            nearestpoweroftwo = nearestpoweroftwo*0.5_dp
                        end do
                        split = nearestpoweroftwo/segmentlength
                        if (acutedest) split = 1.0_dp - split
                    else
                        split = 0.5_dp
                    end if

                    newvertex = vertex_alloc(s)
                    s%vx(newvertex) = s%vx(eorg) + &
                                      split*(s%vx(edest) - s%vx(eorg))
                    s%vy(newvertex) = s%vy(eorg) + &
                                      split*(s%vy(edest) - s%vy(eorg))

                    multiplier = ccw_v(s, eorg, edest, newvertex)
                    divisor = (s%vx(eorg) - s%vx(edest))* &
                              (s%vx(eorg) - s%vx(edest)) + &
                              (s%vy(eorg) - s%vy(edest))* &
                              (s%vy(eorg) - s%vy(edest))
                    if ((multiplier /= 0.0_dp) .and. (divisor /= 0.0_dp)) then
                        multiplier = multiplier/divisor
                        if (multiplier == multiplier) then
                            s%vx(newvertex) = s%vx(newvertex) + &
                                              multiplier* &
                                              (s%vy(edest) - s%vy(eorg))
                            s%vy(newvertex) = s%vy(newvertex) + &
                                              multiplier* &
                                              (s%vx(eorg) - s%vx(edest))
                        end if
                    end if

                    s%vmark(newvertex) = s%smark(s_slot(currentenc))
                    s%vtype(newvertex) = V_SEGMENT
                    if (((s%vx(newvertex) == s%vx(eorg)) .and. &
                         (s%vy(newvertex) == s%vy(eorg))) .or. &
                        ((s%vx(newvertex) == s%vx(edest)) .and. &
                         (s%vy(newvertex) == s%vy(edest)))) then
                        error stop 'tc_enforce: segment split out of precision'
                    end if
                    success = insertvertex(s, newvertex, enctri, currentenc, &
                                           .true., triflaws)
                    if ((success /= IV_SUCCESS) .and. &
                        (success /= IV_ENCROACHING)) then
                        error stop 'tc_enforce: failure to split a segment'
                    end if
                    if (s%steinerleft > 0) s%steinerleft = s%steinerleft - 1
                    dummy_ = checkseg4encroach(s, currentenc)
                    currentenc = s_next(s, currentenc)
                    dummy_ = checkseg4encroach(s, currentenc)
                end if
                call badsubseg_dealloc(s, encidx)
                encidx = encidx + 1
            end do
        end do
    end subroutine splitencsegs

    function chew_dot(s, eorg, edest, eapex) result(d)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: eorg, edest, eapex
        real(dp) :: d
        d = (s%vx(eorg) - s%vx(eapex))*(s%vx(edest) - s%vx(eapex)) + &
            (s%vy(eorg) - s%vy(eapex))*(s%vy(edest) - s%vy(eapex))
    end function chew_dot

    subroutine deletevertex(s, deltri)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: deltri
        integer :: countingtri, firstedge, lastedge, deltriright
        integer :: lefttri, righttri, leftcasing, rightcasing
        integer :: leftsubseg, rightsubseg, delvertex, neworg
        integer :: edgecount

        delvertex = t_org(s, deltri)
        call vertex_dealloc(s, delvertex)

        countingtri = t_onext(s, deltri)
        edgecount = 1
        do while (countingtri /= deltri)
            edgecount = edgecount + 1
            countingtri = t_onext(s, countingtri)
        end do

        if (edgecount > 3) then
            firstedge = t_onext(s, deltri)
            lastedge = t_oprev(s, deltri)
            call triangulatepolygon(s, firstedge, lastedge, edgecount, &
                                    .false., s%nobisect == 0)
        end if
        deltriright = t_lprev(deltri)
        lefttri = t_dnext(s, deltri)
        leftcasing = t_sym(s, lefttri)
        righttri = t_oprev(s, deltriright)
        rightcasing = t_sym(s, righttri)
        call t_bond(s, deltri, leftcasing)
        call t_bond(s, deltriright, rightcasing)
        leftsubseg = ts_pivot(s, lefttri)
        if (s_slot(leftsubseg) /= 0) call ts_bond(s, deltri, leftsubseg)
        rightsubseg = ts_pivot(s, righttri)
        if (s_slot(rightsubseg) /= 0) call ts_bond(s, deltriright, rightsubseg)

        neworg = t_org(s, lefttri)
        call set_org(s, deltri, neworg)
        if (s%nobisect == 0) call testtriangle(s, deltri)

        call triangle_dealloc(s, t_slot(lefttri))
        call triangle_dealloc(s, t_slot(righttri))
    end subroutine deletevertex

    recursive subroutine triangulatepolygon(s, firstedge, lastedge, &
                                            edgecount, doflip, triflaws)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: firstedge
        integer, intent(inout) :: lastedge
        integer, intent(in) :: edgecount
        logical, intent(in) :: doflip, triflaws
        integer :: testtri, besttri, tempedge
        integer :: leftbasevertex, rightbasevertex, testvertex, bestvertex
        integer :: bestnumber, i

        leftbasevertex = t_apex(s, lastedge)
        rightbasevertex = t_dest(s, firstedge)
        besttri = t_onext(s, firstedge)
        bestvertex = t_dest(s, besttri)
        testtri = besttri
        bestnumber = 1
        do i = 2, edgecount - 2
            testtri = t_onext(s, testtri)
            testvertex = t_dest(s, testtri)
            if (incircle_v(s, leftbasevertex, rightbasevertex, bestvertex, &
                           testvertex) > 0.0_dp) then
                besttri = testtri
                bestvertex = testvertex
                bestnumber = i
            end if
        end do
        if (bestnumber > 1) then
            tempedge = t_oprev(s, besttri)
            call triangulatepolygon(s, firstedge, tempedge, bestnumber + 1, &
                                    .true., triflaws)
        end if
        if (bestnumber < edgecount - 2) then
            tempedge = t_sym(s, besttri)
            call triangulatepolygon(s, besttri, lastedge, &
                                    edgecount - bestnumber, .true., triflaws)
            besttri = t_sym(s, tempedge)
        end if
        if (doflip) then
            call flip(s, besttri)
            if (triflaws) then
                testtri = t_sym(s, besttri)
                call testtriangle(s, testtri)
            end if
        end if
        lastedge = besttri
    end subroutine triangulatepolygon

end module tc_enforce
