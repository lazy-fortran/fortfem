module tc_qtest
    !! Quality testing: encroachment checks, bad-triangle classification,
    !! the 4096-bucket priority queue, and the circumcenter/off-center
    !! computation.  Floating-point evaluation order (including FMA
    !! contraction) mirrors the oracle binary so Steiner point coordinates
    !! agree bitwise.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_fma
    use tc_state
    use tc_alloc
    use tc_predicates, only: tc_counterclockwise
    implicit none
    private

    public :: checkseg4encroach, testtriangle, findcircumcenter
    public :: enqueuebadtriang, dequeuebadtriang, badtri_release
    public :: badsubseg_alloc, badsubseg_dealloc, tallyencs, tallyfaces

contains

    function badsubseg_alloc(s) result(idx)
        type(tc_state_t), intent(inout) :: s
        integer :: idx
        integer, allocatable :: i1(:)
        logical, allocatable :: l1(:)
        integer :: oldcap

        if (s%bs_nfree > 0) then
            idx = s%bs_free(s%bs_nfree)
            s%bs_nfree = s%bs_nfree - 1
        else
            idx = s%bs_slots + 1
            oldcap = size(s%bs_sub)
            if (idx > oldcap) then
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bs_sub
                call move_alloc(i1, s%bs_sub)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bs_org
                call move_alloc(i1, s%bs_org)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bs_dest
                call move_alloc(i1, s%bs_dest)
                allocate (l1(2*oldcap))
                l1(1:oldcap) = s%bs_dead
                call move_alloc(l1, s%bs_dead)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bs_free
                call move_alloc(i1, s%bs_free)
            end if
            s%bs_slots = idx
        end if
        s%bs_dead(idx) = .false.
        s%bs_items = s%bs_items + 1
    end function badsubseg_alloc

    subroutine badsubseg_dealloc(s, idx)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: idx
        s%bs_dead(idx) = .true.
        s%bs_items = s%bs_items - 1
        s%bs_nfree = s%bs_nfree + 1
        s%bs_free(s%bs_nfree) = idx
    end subroutine badsubseg_dealloc

    function checkseg4encroach(s, testsubseg) result(encroached)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: testsubseg
        integer :: encroached
        integer :: neighbortri, testsym, eorg, edest, eapex, sides, idx
        real(dp) :: dotproduct, s2, orgdist, destdist

        encroached = 0
        sides = 0
        eorg = s_org(s, testsubseg)
        edest = s_dest(s, testsubseg)
        neighbortri = st_pivot(s, testsubseg)
        if (t_slot(neighbortri) /= 0) then
            sides = sides + 1
            eapex = t_apex(s, neighbortri)
            dotproduct = lens_dot(s, eorg, edest, eapex)
            if (dotproduct < 0.0_dp) then
                s2 = 2.0_dp*s%goodangle - 1.0_dp
                orgdist = sq_dist(s, eorg, eapex)
                destdist = sq_dist(s, edest, eapex)
                if (dotproduct*dotproduct >= destdist*((s2*s2)*orgdist)) then
                    encroached = 1
                end if
            end if
        end if
        testsym = s_sym(testsubseg)
        neighbortri = st_pivot(s, testsym)
        if (t_slot(neighbortri) /= 0) then
            sides = sides + 1
            eapex = t_apex(s, neighbortri)
            dotproduct = lens_dot(s, eorg, edest, eapex)
            if (dotproduct < 0.0_dp) then
                s2 = 2.0_dp*s%goodangle - 1.0_dp
                orgdist = sq_dist(s, eorg, eapex)
                destdist = sq_dist(s, edest, eapex)
                if (dotproduct*dotproduct >= destdist*((s2*s2)*orgdist)) then
                    encroached = encroached + 2
                end if
            end if
        end if

        if (encroached /= 0 .and. (s%nobisect == 0 .or. &
                                   (s%nobisect == 1 .and. sides == 2))) then
            idx = badsubseg_alloc(s)
            if (encroached == 1) then
                s%bs_sub(idx) = testsubseg
                s%bs_org(idx) = eorg
                s%bs_dest(idx) = edest
            else
                s%bs_sub(idx) = testsym
                s%bs_org(idx) = edest
                s%bs_dest(idx) = eorg
            end if
        end if
    end function checkseg4encroach

    function lens_dot(s, eorg, edest, eapex) result(d)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: eorg, edest, eapex
        real(dp) :: d
        d = ieee_fma(s%vx(eorg) - s%vx(eapex), s%vx(edest) - s%vx(eapex), &
                     (s%vy(eorg) - s%vy(eapex))*(s%vy(edest) - s%vy(eapex)))
    end function lens_dot

    function sq_dist(s, va, vb) result(d)
        type(tc_state_t), intent(in) :: s
        integer, intent(in) :: va, vb
        real(dp) :: d
        d = ieee_fma(s%vx(va) - s%vx(vb), s%vx(va) - s%vx(vb), &
                     (s%vy(va) - s%vy(vb))*(s%vy(va) - s%vy(vb)))
    end function sq_dist

    subroutine testtriangle(s, testtri)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: testtri
        integer :: torg, tdest, tapex, base1, base2, tri1, tri2, testsub
        integer :: org1, dest1, org2, dest2, joinvertex
        real(dp) :: dxod, dyod, dxda, dyda, dxao, dyao
        real(dp), volatile :: dxod2, dyod2, dxda2, dyda2, dxao2, dyao2
        real(dp) :: apexlen, orglen, destlen, minedge, angle, area
        real(dp) :: dist1, dist2

        torg = t_org(s, testtri)
        tdest = t_dest(s, testtri)
        tapex = t_apex(s, testtri)
        dxod = s%vx(torg) - s%vx(tdest)
        dyod = s%vy(torg) - s%vy(tdest)
        dxda = s%vx(tdest) - s%vx(tapex)
        dyda = s%vy(tdest) - s%vy(tapex)
        dxao = s%vx(tapex) - s%vx(torg)
        dyao = s%vy(tapex) - s%vy(torg)
        dxod2 = dxod*dxod
        dyod2 = dyod*dyod
        dxda2 = dxda*dxda
        dyda2 = dyda*dyda
        dxao2 = dxao*dxao
        dyao2 = dyao*dyao
        apexlen = dxod2 + dyod2
        orglen = dxda2 + dyda2
        destlen = dxao2 + dyao2

        if ((apexlen < orglen) .and. (apexlen < destlen)) then
            minedge = apexlen
            angle = ieee_fma(dxda, dxao, dyda*dyao)
            angle = angle*angle/(orglen*destlen)
            base1 = torg
            base2 = tdest
            tri1 = testtri
        else if (orglen < destlen) then
            minedge = orglen
            angle = ieee_fma(dxod, dxao, dyod*dyao)
            angle = angle*angle/(apexlen*destlen)
            base1 = tdest
            base2 = tapex
            tri1 = t_lnext(testtri)
        else
            minedge = destlen
            angle = ieee_fma(dxod, dxda, dyod*dyda)
            angle = angle*angle/(apexlen*orglen)
            base1 = tapex
            base2 = torg
            tri1 = t_lprev(testtri)
        end if

        if (s%fixedarea) then
            area = 0.5_dp*(dxod*dyda - dyod*dxda)
            if (area > s%maxarea) then
                call enqueuebadtri(s, testtri, minedge, tapex, torg, tdest)
                return
            end if
        end if

        if (angle > s%goodangle) then
            if ((s%vtype(base1) == V_SEGMENT) .and. &
                (s%vtype(base2) == V_SEGMENT)) then
                testsub = ts_pivot(s, tri1)
                if (s_slot(testsub) == 0) then
                    tri2 = tri1
                    do
                        tri1 = t_oprev(s, tri1)
                        testsub = ts_pivot(s, tri1)
                        if (s_slot(testsub) /= 0) exit
                    end do
                    org1 = seg_org(s, testsub)
                    dest1 = seg_dest(s, testsub)
                    do
                        tri2 = t_dnext(s, tri2)
                        testsub = ts_pivot(s, tri2)
                        if (s_slot(testsub) /= 0) exit
                    end do
                    org2 = seg_org(s, testsub)
                    dest2 = seg_dest(s, testsub)
                    joinvertex = 0
                    if ((s%vx(dest1) == s%vx(org2)) .and. &
                        (s%vy(dest1) == s%vy(org2))) then
                        joinvertex = dest1
                    else if ((s%vx(org1) == s%vx(dest2)) .and. &
                             (s%vy(org1) == s%vy(dest2))) then
                        joinvertex = org1
                    end if
                    if (joinvertex /= 0) then
                        dist1 = sq_dist(s, base1, joinvertex)
                        dist2 = sq_dist(s, base2, joinvertex)
                        if ((dist1 < 1.001_dp*dist2) .and. &
                            (dist1 > 0.999_dp*dist2)) return
                    end if
                end if
            end if
            call enqueuebadtri(s, testtri, minedge, tapex, torg, tdest)
        end if
    end subroutine testtriangle

    subroutine enqueuebadtri(s, enqtri, minedge, enqapex, enqorg, enqdest)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: enqtri, enqapex, enqorg, enqdest
        real(dp), intent(in) :: minedge
        integer :: idx

        idx = badtri_alloc(s)
        s%bt_tri(idx) = enqtri
        s%bt_key(idx) = minedge
        s%bt_apex(idx) = enqapex
        s%bt_org(idx) = enqorg
        s%bt_dest(idx) = enqdest
        call enqueuebadtriang(s, idx)
    end subroutine enqueuebadtri

    function badtri_alloc(s) result(idx)
        type(tc_state_t), intent(inout) :: s
        integer :: idx
        integer, allocatable :: i1(:)
        real(dp), allocatable :: r1(:)
        integer :: oldcap

        if (s%bt_nfree > 0) then
            idx = s%bt_free(s%bt_nfree)
            s%bt_nfree = s%bt_nfree - 1
        else
            idx = s%bt_slots + 1
            oldcap = size(s%bt_tri)
            if (idx > oldcap) then
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bt_tri
                call move_alloc(i1, s%bt_tri)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bt_org
                call move_alloc(i1, s%bt_org)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bt_dest
                call move_alloc(i1, s%bt_dest)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bt_apex
                call move_alloc(i1, s%bt_apex)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bt_next
                call move_alloc(i1, s%bt_next)
                allocate (r1(2*oldcap))
                r1(1:oldcap) = s%bt_key
                call move_alloc(r1, s%bt_key)
                allocate (i1(2*oldcap))
                i1(1:oldcap) = s%bt_free
                call move_alloc(i1, s%bt_free)
            end if
            s%bt_slots = idx
        end if
        s%bt_items = s%bt_items + 1
    end function badtri_alloc

    subroutine badtri_release(s, idx)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: idx
        s%bt_items = s%bt_items - 1
        s%bt_nfree = s%bt_nfree + 1
        s%bt_free(s%bt_nfree) = idx
    end subroutine badtri_release

    subroutine enqueuebadtriang(s, idx)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: idx
        real(dp) :: length, multiplier
        integer :: exponent_, expincrement, queuenumber, i
        logical :: posexponent

        if (s%bt_key(idx) >= 1.0_dp) then
            length = s%bt_key(idx)
            posexponent = .true.
        else
            length = 1.0_dp/s%bt_key(idx)
            posexponent = .false.
        end if
        exponent_ = 0
        do while (length > 2.0_dp)
            expincrement = 1
            multiplier = 0.5_dp
            do while (length*multiplier*multiplier > 1.0_dp)
                expincrement = expincrement*2
                multiplier = multiplier*multiplier
            end do
            exponent_ = exponent_ + expincrement
            length = length*multiplier
        end do
        exponent_ = 2*exponent_ + merge(1, 0, length > SQRT2_TC)
        if (posexponent) then
            queuenumber = 2047 - exponent_
        else
            queuenumber = 2048 + exponent_
        end if

        if (s%btq_front(queuenumber) < 0) then
            if (queuenumber > s%btq_first) then
                s%btq_nextq(queuenumber) = s%btq_first
                s%btq_first = queuenumber
            else
                i = queuenumber + 1
                do while (s%btq_front(i) < 0)
                    i = i + 1
                end do
                s%btq_nextq(queuenumber) = s%btq_nextq(i)
                s%btq_nextq(i) = queuenumber
            end if
            s%btq_front(queuenumber) = idx
        else
            s%bt_next(s%btq_tail(queuenumber)) = idx
        end if
        s%btq_tail(queuenumber) = idx
        s%bt_next(idx) = -1
    end subroutine enqueuebadtriang

    function dequeuebadtriang(s) result(idx)
        type(tc_state_t), intent(inout) :: s
        integer :: idx

        if (s%btq_first < 0) then
            idx = -1
            return
        end if
        idx = s%btq_front(s%btq_first)
        s%btq_front(s%btq_first) = s%bt_next(idx)
        if (idx == s%btq_tail(s%btq_first)) then
            s%btq_first = s%btq_nextq(s%btq_first)
        end if
    end function dequeuebadtriang

    subroutine tallyencs(s)
        type(tc_state_t), intent(inout) :: s
        integer :: slot, dummy_

        do slot = 1, s%ns_slots
            if (s%sdead(slot)) cycle
            dummy_ = checkseg4encroach(s, s_make(slot, 0))
        end do
    end subroutine tallyencs

    subroutine tallyfaces(s)
        type(tc_state_t), intent(inout) :: s
        integer :: slot

        do slot = 1, s%nt_slots
            if (s%tdead(slot)) cycle
            call testtriangle(s, t_make(slot, 0))
        end do
    end subroutine tallyfaces

    subroutine findcircumcenter(s, torg, tdest, tapex, ccx, ccy, xi, eta, &
                                offcenter)
        type(tc_state_t), intent(inout) :: s
        integer, intent(in) :: torg, tdest, tapex
        real(dp), intent(out) :: ccx, ccy, xi, eta
        logical, intent(in) :: offcenter
        real(dp) :: xdo, ydo, xao, yao, dodist, aodist, dadist
        real(dp) :: denominator, dx, dy, dxoff, dyoff, dxa, dya

        xdo = s%vx(tdest) - s%vx(torg)
        ydo = s%vy(tdest) - s%vy(torg)
        xao = s%vx(tapex) - s%vx(torg)
        yao = s%vy(tapex) - s%vy(torg)
        dodist = ieee_fma(xdo, xdo, ydo*ydo)
        aodist = ieee_fma(xao, xao, yao*yao)
        dadist = ieee_fma(s%vx(tdest) - s%vx(tapex), s%vx(tdest) - s%vx(tapex), &
                          (s%vy(tdest) - s%vy(tapex))* &
                          (s%vy(tdest) - s%vy(tapex)))
        denominator = 0.5_dp/tc_counterclockwise(s%vx(tdest), s%vy(tdest), &
                                                 s%vx(tapex), s%vy(tapex), &
                                                 s%vx(torg), s%vy(torg))
        dx = ieee_fma(yao, dodist, -(ydo*aodist))*denominator
        dy = ieee_fma(xdo, aodist, -(xao*dodist))*denominator

        if ((dodist < aodist) .and. (dodist < dadist)) then
            if (offcenter .and. (s%offconstant > 0.0_dp)) then
                dxoff = ieee_fma(0.5_dp, xdo, -(s%offconstant*ydo))
                dyoff = ieee_fma(0.5_dp, ydo, s%offconstant*xdo)
                if (ieee_fma(dxoff, dxoff, dyoff*dyoff) < &
                    ieee_fma(dx, dx, dy*dy)) then
                    dx = dxoff
                    dy = dyoff
                end if
            end if
        else if (aodist < dadist) then
            if (offcenter .and. (s%offconstant > 0.0_dp)) then
                dxoff = ieee_fma(0.5_dp, xao, s%offconstant*yao)
                dyoff = ieee_fma(0.5_dp, yao, -(s%offconstant*xao))
                if (ieee_fma(dxoff, dxoff, dyoff*dyoff) < &
                    ieee_fma(dx, dx, dy*dy)) then
                    dx = dxoff
                    dy = dyoff
                end if
            end if
        else
            if (offcenter .and. (s%offconstant > 0.0_dp)) then
                dxa = s%vx(tapex) - s%vx(tdest)
                dya = s%vy(tapex) - s%vy(tdest)
                dxoff = ieee_fma(0.5_dp, dxa, -(s%offconstant*dya))
                dyoff = ieee_fma(0.5_dp, dya, s%offconstant*dxa)
                if (ieee_fma(dxoff, dxoff, dyoff*dyoff) < &
                    ieee_fma(dx - xdo, dx - xdo, (dy - ydo)*(dy - ydo))) then
                    dx = xdo + dxoff
                    dy = ydo + dyoff
                end if
            end if
        end if

        ccx = s%vx(torg) + dx
        ccy = s%vy(torg) + dy
        xi = (denominator + denominator)*ieee_fma(yao, dx, -(xao*dy))
        eta = (denominator + denominator)*ieee_fma(xdo, dy, -(ydo*dx))
    end subroutine findcircumcenter

end module tc_qtest
