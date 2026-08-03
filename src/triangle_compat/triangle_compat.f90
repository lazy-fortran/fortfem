module triangle_compat
    !! Triangle-compatible meshing pipeline: divide-and-conquer constrained
    !! Delaunay triangulation, hole carving, and quality refinement matching
    !! the behavior of the Triangle mesh generator with flags "pqY" (PSLG,
    !! quality bound, no Steiner points on boundary segments).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tc_state
    use tc_alloc, only: vertex_alloc
    use tc_divconq, only: divconq_delaunay
    use tc_segment, only: formskeleton
    use tc_carve, only: carveholes
    use tc_enforce, only: enforcequality
    implicit none
    private

    public :: tc_result_t, triangulate_compat

    real(dp), parameter :: PI_TC = 3.141592653589793_dp

    type :: tc_result_t
        integer :: npoints = 0
        integer :: ntriangles = 0
        integer :: nsegments = 0
        real(dp), allocatable :: points(:, :)
        integer, allocatable :: pointmarkers(:)
        integer, allocatable :: triangles(:, :)
        integer, allocatable :: segments(:, :)
        integer, allocatable :: segmentmarkers(:)
        integer, allocatable :: neighbors(:, :)
        real(dp) :: min_angle_deg = 0.0_dp
    end type tc_result_t

    type :: behavior_t
        real(dp) :: min_angle = 20.0_dp
        real(dp) :: max_area = -1.0_dp
        logical :: quality = .true.
        integer :: nobisect = 1
    end type behavior_t

contains

    subroutine triangulate_compat(points, segments, holes, res, stat, &
                                  min_angle, max_area, quality, nobisect, &
                                  segmentmarkers)
        real(dp), intent(in) :: points(:, :)
        integer, intent(in) :: segments(:, :)
        real(dp), intent(in) :: holes(:, :)
        type(tc_result_t), intent(out) :: res
        integer, intent(out) :: stat
        real(dp), intent(in), optional :: min_angle, max_area
        logical, intent(in), optional :: quality
        integer, intent(in), optional :: nobisect
        integer, intent(in), optional :: segmentmarkers(:)

        type(tc_state_t) :: s
        type(behavior_t) :: b
        integer, allocatable :: segmarks(:)
        integer :: i, v

        stat = 0
        if (size(points, 2) < 3) then
            stat = 1
            return
        end if
        if (present(min_angle)) b%min_angle = min_angle
        if (present(max_area)) b%max_area = max_area
        if (present(quality)) b%quality = quality
        if (present(nobisect)) b%nobisect = nobisect
        if (present(segmentmarkers)) then
            segmarks = segmentmarkers
        else
            allocate (segmarks(0))
        end if

        call tc_init_state(s, size(points, 2))
        s%invertices = size(points, 2)
        do i = 1, s%invertices
            v = vertex_alloc(s)
            s%vx(v) = points(1, i)
            s%vy(v) = points(2, i)
            s%vmark(v) = 0
            s%vtype(v) = V_INPUT
            if (i == 1) then
                s%xmin = points(1, i)
                s%xmax = points(1, i)
                s%ymin = points(2, i)
                s%ymax = points(2, i)
            else
                s%xmin = min(s%xmin, points(1, i))
                s%xmax = max(s%xmax, points(1, i))
                s%ymin = min(s%ymin, points(2, i))
                s%ymax = max(s%ymax, points(2, i))
            end if
        end do

        if (b%quality) then
            s%minangle = b%min_angle
        else
            s%minangle = 0.0_dp
        end if
        s%goodangle = cos(s%minangle*PI_TC/180.0_dp)
        if (s%goodangle == 1.0_dp) then
            s%offconstant = 0.0_dp
        else
            s%offconstant = 0.475_dp*sqrt((1.0_dp + s%goodangle)/ &
                                          (1.0_dp - s%goodangle))
        end if
        s%goodangle = s%goodangle*s%goodangle
        s%nobisect = b%nobisect
        s%fixedarea = b%max_area > 0.0_dp
        s%maxarea = b%max_area
        s%steinerleft = -1

        s%hullsize = divconq_delaunay(s)

        s%checksegments = .true.
        call formskeleton(s, segments, segmarks)

        if (s%nt_live > 0) then
            call carveholes(s, holes)
        end if

        if (b%quality .and. (s%nt_live > 0)) then
            call enforcequality(s)
        end if

        call extract_result(s, res)
    end subroutine triangulate_compat

    subroutine extract_result(s, res)
        type(tc_state_t), intent(inout) :: s
        type(tc_result_t), intent(out) :: res
        integer, allocatable :: vmap(:)
        integer :: i, n, orient, nb, os

        allocate (vmap(0:s%nv_slots))
        vmap = 0
        n = 0
        do i = 1, s%nv_slots
            if ((s%vtype(i) == V_DEAD) .or. (s%vtype(i) == V_UNDEAD)) cycle
            n = n + 1
            vmap(i) = n
        end do
        res%npoints = n
        allocate (res%points(2, n), res%pointmarkers(n))
        do i = 1, s%nv_slots
            if (vmap(i) == 0) cycle
            res%points(1, vmap(i)) = s%vx(i)
            res%points(2, vmap(i)) = s%vy(i)
            res%pointmarkers(vmap(i)) = s%vmark(i)
        end do

        n = 0
        do i = 1, s%nt_slots
            if (.not. s%tdead(i)) n = n + 1
        end do
        res%ntriangles = n
        allocate (res%triangles(3, n), res%neighbors(3, n))
        ! Triangle output numbering of live slots, for the neighbor list.
        block
            integer, allocatable :: tmap(:)
            allocate (tmap(0:s%nt_slots))
            tmap = 0
            n = 0
            do i = 1, s%nt_slots
                if (s%tdead(i)) cycle
                n = n + 1
                tmap(i) = n
            end do
            n = 0
            do i = 1, s%nt_slots
                if (s%tdead(i)) cycle
                n = n + 1
                res%triangles(1, n) = vmap(t_org(s, t_make(i, 0)))
                res%triangles(2, n) = vmap(t_dest(s, t_make(i, 0)))
                res%triangles(3, n) = vmap(t_apex(s, t_make(i, 0)))
                do orient = 0, 2
                    nb = t_sym(s, t_make(i, orient))
                    if (t_slot(nb) == 0) then
                        res%neighbors(orient + 1, n) = -1
                    else
                        res%neighbors(orient + 1, n) = tmap(t_slot(nb))
                    end if
                end do
            end do
        end block

        n = 0
        do i = 1, s%ns_slots
            if (.not. s%sdead(i)) n = n + 1
        end do
        res%nsegments = n
        allocate (res%segments(2, n), res%segmentmarkers(n))
        n = 0
        do i = 1, s%ns_slots
            if (s%sdead(i)) cycle
            n = n + 1
            os = s_make(i, 0)
            res%segments(1, n) = vmap(s_org(s, os))
            res%segments(2, n) = vmap(s_dest(s, os))
            res%segmentmarkers(n) = s%smark(i)
        end do

        res%min_angle_deg = mesh_min_angle(res)
    end subroutine extract_result

    function mesh_min_angle(res) result(amin)
        type(tc_result_t), intent(in) :: res
        real(dp) :: amin
        real(dp) :: ax, ay, bx, by, cx, cy, la, lb, lc, cosang
        integer :: i

        amin = 180.0_dp
        do i = 1, res%ntriangles
            ax = res%points(1, res%triangles(1, i))
            ay = res%points(2, res%triangles(1, i))
            bx = res%points(1, res%triangles(2, i))
            by = res%points(2, res%triangles(2, i))
            cx = res%points(1, res%triangles(3, i))
            cy = res%points(2, res%triangles(3, i))
            la = (bx - cx)**2 + (by - cy)**2
            lb = (ax - cx)**2 + (ay - cy)**2
            lc = (ax - bx)**2 + (ay - by)**2
            cosang = (lb + lc - la)/(2.0_dp*sqrt(lb*lc))
            amin = min(amin, acos(max(-1.0_dp, min(1.0_dp, cosang)))* &
                       180.0_dp/PI_TC)
            cosang = (la + lc - lb)/(2.0_dp*sqrt(la*lc))
            amin = min(amin, acos(max(-1.0_dp, min(1.0_dp, cosang)))* &
                       180.0_dp/PI_TC)
            cosang = (la + lb - lc)/(2.0_dp*sqrt(la*lb))
            amin = min(amin, acos(max(-1.0_dp, min(1.0_dp, cosang)))* &
                       180.0_dp/PI_TC)
        end do
    end function mesh_min_angle

end module triangle_compat
