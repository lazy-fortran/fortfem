module tc_locate
    !! Point location by random sampling plus directed walk, mirroring the
    !! reference implementation including its sampling RNG consumption.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_fma
    use tc_state
    use tc_geom, only: ccw_p
    implicit none
    private

    public :: makevertexmap, preciselocate, locate, finddirection
    public :: LOC_INTRIANGLE, LOC_ONEDGE, LOC_ONVERTEX, LOC_OUTSIDE
    public :: DIR_WITHIN, DIR_LEFTCOLLINEAR, DIR_RIGHTCOLLINEAR

    integer, parameter :: LOC_INTRIANGLE = 1, LOC_ONEDGE = 2
    integer, parameter :: LOC_ONVERTEX = 3, LOC_OUTSIDE = 4
    integer, parameter :: DIR_WITHIN = 1, DIR_LEFTCOLLINEAR = 2
    integer, parameter :: DIR_RIGHTCOLLINEAR = 3

contains

    subroutine makevertexmap(s)
        type(tc_state_t), intent(inout) :: s
        integer :: slot, orient

        do slot = 1, s%nt_slots
            if (s%tdead(slot)) cycle
            do orient = 0, 2
                s%vtri(t_org(s, t_make(slot, orient))) = t_make(slot, orient)
            end do
        end do
    end subroutine makevertexmap

    function preciselocate(s, px, py, searchtri, stopatsubsegment) result(code)
        type(tc_state_t), intent(inout) :: s
        real(dp), intent(in) :: px, py
        integer, intent(inout) :: searchtri
        logical, intent(in) :: stopatsubsegment
        integer :: code
        integer :: forg, fdest, fapex, backtracktri
        real(dp) :: orgorient, destorient
        logical :: moveleft

        forg = t_org(s, searchtri)
        fdest = t_dest(s, searchtri)
        fapex = t_apex(s, searchtri)
        do
            if ((s%vx(fapex) == px) .and. (s%vy(fapex) == py)) then
                searchtri = t_lprev(searchtri)
                code = LOC_ONVERTEX
                return
            end if
            destorient = ccw_p(s, forg, fapex, px, py)
            orgorient = ccw_p(s, fapex, fdest, px, py)
            if (destorient > 0.0_dp) then
                if (orgorient > 0.0_dp) then
                    moveleft = ieee_fma(s%vx(fapex) - px, &
                                        s%vx(fdest) - s%vx(forg), &
                                        (s%vy(fapex) - py)* &
                                        (s%vy(fdest) - s%vy(forg))) > 0.0_dp
                else
                    moveleft = .true.
                end if
            else
                if (orgorient > 0.0_dp) then
                    moveleft = .false.
                else
                    if (destorient == 0.0_dp) then
                        searchtri = t_lprev(searchtri)
                        code = LOC_ONEDGE
                        return
                    end if
                    if (orgorient == 0.0_dp) then
                        searchtri = t_lnext(searchtri)
                        code = LOC_ONEDGE
                        return
                    end if
                    code = LOC_INTRIANGLE
                    return
                end if
            end if

            if (moveleft) then
                backtracktri = t_lprev(searchtri)
                fdest = fapex
            else
                backtracktri = t_lnext(searchtri)
                forg = fapex
            end if
            searchtri = t_sym(s, backtracktri)

            if (s%checksegments .and. stopatsubsegment) then
                if (s_slot(ts_pivot(s, backtracktri)) /= 0) then
                    searchtri = backtracktri
                    code = LOC_OUTSIDE
                    return
                end if
            end if
            if (t_slot(searchtri) == 0) then
                searchtri = backtracktri
                code = LOC_OUTSIDE
                return
            end if

            fapex = t_apex(s, searchtri)
        end do
    end function preciselocate

    function locate(s, px, py, searchtri) result(code)
        type(tc_state_t), intent(inout) :: s
        real(dp), intent(in) :: px, py
        integer, intent(inout) :: searchtri
        integer :: code
        integer :: torg, tdest, sampleslot, sampleorg
        real(dp) :: searchdist, dist, ahead
        integer :: samplesperblock, totalsamplesleft, samplesleft
        integer :: population, totalpopulation, blockstart

        torg = t_org(s, searchtri)
        searchdist = ieee_fma(px - s%vx(torg), px - s%vx(torg), &
                              (py - s%vy(torg))*(py - s%vy(torg)))

        if (s%recent >= 0) then
            if (.not. s%tdead(t_slot(s%recent))) then
                torg = t_org(s, s%recent)
                if ((s%vx(torg) == px) .and. (s%vy(torg) == py)) then
                    searchtri = s%recent
                    code = LOC_ONVERTEX
                    return
                end if
                dist = ieee_fma(px - s%vx(torg), px - s%vx(torg), &
                                (py - s%vy(torg))*(py - s%vy(torg)))
                if (dist < searchdist) then
                    searchtri = s%recent
                    searchdist = dist
                end if
            end if
        end if

        do while (SAMPLEFACTOR*s%samples*s%samples*s%samples < s%nt_live)
            s%samples = s%samples + 1
        end do

        samplesperblock = (s%samples*TRIPERBLOCK - 1)/s%nt_slots + 1
        samplesleft = (s%samples*TRIPERBLOCK - 1)/s%nt_slots + 1
        totalsamplesleft = s%samples
        population = TRIPERBLOCK
        totalpopulation = s%nt_slots
        blockstart = 0
        do while (totalsamplesleft > 0)
            if (population > totalpopulation) population = totalpopulation
            do
                sampleslot = blockstart + tc_random(s, population) + 1
                if (.not. s%tdead(sampleslot)) then
                    sampleorg = t_org(s, t_make(sampleslot, 0))
                    dist = ieee_fma(px - s%vx(sampleorg), px - s%vx(sampleorg), &
                                    (py - s%vy(sampleorg))* &
                                    (py - s%vy(sampleorg)))
                    if (dist < searchdist) then
                        searchtri = t_make(sampleslot, 0)
                        searchdist = dist
                    end if
                end if
                samplesleft = samplesleft - 1
                totalsamplesleft = totalsamplesleft - 1
                if ((samplesleft <= 0) .or. (totalsamplesleft <= 0)) exit
            end do
            if (totalsamplesleft > 0) then
                blockstart = blockstart + population
                samplesleft = samplesperblock
                totalpopulation = totalpopulation - population
                population = TRIPERBLOCK
            end if
        end do

        torg = t_org(s, searchtri)
        tdest = t_dest(s, searchtri)
        if ((s%vx(torg) == px) .and. (s%vy(torg) == py)) then
            code = LOC_ONVERTEX
            return
        end if
        if ((s%vx(tdest) == px) .and. (s%vy(tdest) == py)) then
            searchtri = t_lnext(searchtri)
            code = LOC_ONVERTEX
            return
        end if
        ahead = ccw_p(s, torg, tdest, px, py)
        if (ahead < 0.0_dp) then
            searchtri = t_sym(s, searchtri)
        else if (ahead == 0.0_dp) then
            if (((s%vx(torg) < px) .eqv. (px < s%vx(tdest))) .and. &
                ((s%vy(torg) < py) .eqv. (py < s%vy(tdest)))) then
                code = LOC_ONEDGE
                return
            end if
        end if
        code = preciselocate(s, px, py, searchtri, .false.)
    end function locate

    function finddirection(s, searchtri, searchpoint) result(code)
        type(tc_state_t), intent(inout) :: s
        integer, intent(inout) :: searchtri
        integer, intent(in) :: searchpoint
        integer :: code
        integer :: startvertex, leftvertex, rightvertex, checktri
        real(dp) :: leftccw, rightccw
        logical :: leftflag, rightflag

        startvertex = t_org(s, searchtri)
        rightvertex = t_dest(s, searchtri)
        leftvertex = t_apex(s, searchtri)
        leftccw = ccw_p(s, searchpoint, startvertex, &
                        s%vx(leftvertex), s%vy(leftvertex))
        leftflag = leftccw > 0.0_dp
        rightccw = ccw_p(s, startvertex, searchpoint, &
                         s%vx(rightvertex), s%vy(rightvertex))
        rightflag = rightccw > 0.0_dp
        if (leftflag .and. rightflag) then
            checktri = t_onext(s, searchtri)
            if (t_slot(checktri) == 0) then
                leftflag = .false.
            else
                rightflag = .false.
            end if
        end if
        do while (leftflag)
            searchtri = t_onext(s, searchtri)
            if (t_slot(searchtri) == 0) then
                error stop 'tc_locate: finddirection fell off the mesh'
            end if
            leftvertex = t_apex(s, searchtri)
            rightccw = leftccw
            leftccw = ccw_p(s, searchpoint, startvertex, &
                            s%vx(leftvertex), s%vy(leftvertex))
            leftflag = leftccw > 0.0_dp
        end do
        do while (rightflag)
            searchtri = t_oprev(s, searchtri)
            if (t_slot(searchtri) == 0) then
                error stop 'tc_locate: finddirection fell off the mesh'
            end if
            rightvertex = t_dest(s, searchtri)
            leftccw = rightccw
            rightccw = ccw_p(s, startvertex, searchpoint, &
                             s%vx(rightvertex), s%vy(rightvertex))
            rightflag = rightccw > 0.0_dp
        end do
        if (leftccw == 0.0_dp) then
            code = DIR_LEFTCOLLINEAR
        else if (rightccw == 0.0_dp) then
            code = DIR_RIGHTCOLLINEAR
        else
            code = DIR_WITHIN
        end if
    end function finddirection

end module tc_locate
