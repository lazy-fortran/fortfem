module robust_predicates
    !> Robust geometric predicates using integer coordinate arithmetic.
    !
    !  This module implements exact geometric predicates by mapping floating-point
    !  coordinates to a large integer grid. This approach, used by FreeFEM/BAMG,
    !  provides topological robustness without the complexity of adaptive
    !  floating-point arithmetic.
    !
    !  The key insight is that orientation and incircle tests reduce to computing
    !  determinants. By using integer coordinates, these determinants are computed
    !  exactly (no rounding errors), giving robust topological decisions.
    !
    !  Integer coordinate strategy:
    !    1. Compute bounding box of all points
    !    2. Scale coordinates to fit in 30-bit integers (max ~10^9)
    !    3. Use 64-bit integers for determinant products (30+30=60 bits, fits in 63)
    !    4. Orientation test: 2x2 determinant (exact in 64-bit)
    !    5. Incircle test: 4x4 determinant - needs care to avoid overflow
    !
    !  Reference: FreeFEM/BAMG uses MaxICoor = 2^30 - 1 = 1073741823
    !
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use fortfem_kinds, only: dp, qp
    implicit none

    private
    public :: robust_coords_t
    public :: init_robust_coords
    public :: to_integer_coords
    public :: to_real_coords
    public :: orient2d_robust
    public :: incircle_robust
    public :: ORIENT_CCW, ORIENT_CW, ORIENT_COLLINEAR

    ! Orientation results
    integer, parameter :: ORIENT_CCW = 1        ! Counter-clockwise (positive area)
    integer, parameter :: ORIENT_CW = -1        ! Clockwise (negative area)
    integer, parameter :: ORIENT_COLLINEAR = 0  ! Collinear (zero area)

    ! Maximum integer coordinate value
    ! Increased to 10^8 to provide sufficient resolution even with large margins.
    ! Incircle determinant involves terms of order coord^4 (~10^32).
    ! This requires Quad Precision (qp, ~33 decimal digits) for exact calculation.
    integer(int64), parameter :: MAX_ICOOR = 100000000_int64

    !> Robust coordinate system parameters
    type :: robust_coords_t
        real(dp) :: xmin, ymin       ! Bounding box minimum
        real(dp) :: xmax, ymax       ! Bounding box maximum
        real(dp) :: scale            ! Scale factor: integer = scale * (real - min)
        logical :: initialized = .false.
    end type robust_coords_t

contains

    subroutine init_robust_coords(rc, points, npoints)
        !> Initialize robust coordinate system from point set.
        !
        !  Computes bounding box and scale factor to map real coordinates
        !  to integers in range [0, MAX_ICOOR].
        !
        type(robust_coords_t), intent(out) :: rc
        real(dp), intent(in) :: points(:,:)  ! (2, npoints)
        integer, intent(in) :: npoints

        real(dp) :: dx, dy, max_extent
        integer :: i

        if (npoints < 1) then
            rc%xmin = 0.0_dp
            rc%ymin = 0.0_dp
            rc%xmax = 1.0_dp
            rc%ymax = 1.0_dp
            rc%scale = real(MAX_ICOOR, dp)
            rc%initialized = .true.
            return
        end if

        ! Find bounding box
        rc%xmin = points(1, 1)
        rc%xmax = points(1, 1)
        rc%ymin = points(2, 1)
        rc%ymax = points(2, 1)

        do i = 2, npoints
            if (points(1, i) < rc%xmin) rc%xmin = points(1, i)
            if (points(1, i) > rc%xmax) rc%xmax = points(1, i)
            if (points(2, i) < rc%ymin) rc%ymin = points(2, i)
            if (points(2, i) > rc%ymax) rc%ymax = points(2, i)
        end do

        ! Add small margin to avoid edge cases
        dx = rc%xmax - rc%xmin
        dy = rc%ymax - rc%ymin
        max_extent = max(dx, dy)

        if (max_extent < 1.0e-14_dp) then
            ! Degenerate case: all points at same location
            max_extent = 1.0_dp
        end if

        ! Add large margin to accommodate super-triangle vertices
        ! which are placed far outside the data bounding box.
        rc%xmin = rc%xmin - 100.0_dp * max_extent
        rc%ymin = rc%ymin - 100.0_dp * max_extent
        rc%xmax = rc%xmax + 100.0_dp * max_extent
        rc%ymax = rc%ymax + 100.0_dp * max_extent

        ! Compute scale factor
        dx = rc%xmax - rc%xmin
        dy = rc%ymax - rc%ymin
        max_extent = max(dx, dy)

        rc%scale = real(MAX_ICOOR, dp) / max_extent
        rc%initialized = .true.
    end subroutine init_robust_coords

    subroutine to_integer_coords(rc, x, y, ix, iy)
        !> Convert real coordinates to integer coordinates.
        type(robust_coords_t), intent(in) :: rc
        real(dp), intent(in) :: x, y
        integer(int64), intent(out) :: ix, iy

        ix = int(rc%scale * (x - rc%xmin), int64)
        iy = int(rc%scale * (y - rc%ymin), int64)

        ! Clamp to valid range
        ix = max(0_int64, min(MAX_ICOOR, ix))
        iy = max(0_int64, min(MAX_ICOOR, iy))
    end subroutine to_integer_coords

    subroutine to_real_coords(rc, ix, iy, x, y)
        !> Convert integer coordinates back to real coordinates.
        type(robust_coords_t), intent(in) :: rc
        integer(int64), intent(in) :: ix, iy
        real(dp), intent(out) :: x, y

        x = real(ix, dp) / rc%scale + rc%xmin
        y = real(iy, dp) / rc%scale + rc%ymin
    end subroutine to_real_coords

    integer function orient2d_robust(rc, ax, ay, bx, by, cx, cy) result(orient)
        !> Robust orientation test using integer arithmetic.
        type(robust_coords_t), intent(in) :: rc
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy

        integer(int64) :: iax, iay, ibx, iby, icx, icy
        integer(int64) :: bax, bay, cax, cay
        integer(int64) :: det

        call to_integer_coords(rc, ax, ay, iax, iay)
        call to_integer_coords(rc, bx, by, ibx, iby)
        call to_integer_coords(rc, cx, cy, icx, icy)

        bax = ibx - iax
        bay = iby - iay
        cax = icx - iax
        cay = icy - iay

        det = bax * cay - bay * cax

        if (det > 0_int64) then
            orient = ORIENT_CCW
        else if (det < 0_int64) then
            orient = ORIENT_CW
        else
            orient = ORIENT_COLLINEAR
        end if
    end function orient2d_robust

    logical function incircle_robust(rc, ax, ay, bx, by, cx, cy, dx, dy)                &
        result(inside)
        !> Robust incircle test using integer coordinates.
        !  Uses real(qp) for final determinant to avoid int64 overflow.
        type(robust_coords_t), intent(in) :: rc
        real(dp), intent(in) :: ax, ay, bx, by, cx, cy, dx, dy

        integer(int64) :: iax, iay, ibx, iby, icx, icy, idx, idy
        integer(int64) :: adx, ady, bdx, bdy, cdx, cdy
        real(qp) :: alift, blift, clift
        real(qp) :: bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady
        real(qp) :: det

        call to_integer_coords(rc, ax, ay, iax, iay)
        call to_integer_coords(rc, bx, by, ibx, iby)
        call to_integer_coords(rc, cx, cy, icx, icy)
        call to_integer_coords(rc, dx, dy, idx, idy)

        adx = iax - idx
        ady = iay - idy
        bdx = ibx - idx
        bdy = iby - idy
        cdx = icx - idx
        cdy = icy - idy

        ! Use real(qp) for lift calculation to avoid overflow
        alift = real(adx, qp)**2 + real(ady, qp)**2
        blift = real(bdx, qp)**2 + real(bdy, qp)**2
        clift = real(cdx, qp)**2 + real(cdy, qp)**2

        ! Compute minors using real(qp)
        bdxcdy = real(bdx, qp) * real(cdy, qp)
        cdxbdy = real(cdx, qp) * real(bdy, qp)
        cdxady = real(cdx, qp) * real(ady, qp)
        adxcdy = real(adx, qp) * real(cdy, qp)
        adxbdy = real(adx, qp) * real(bdy, qp)
        bdxady = real(bdx, qp) * real(ady, qp)

        det = alift * (bdxcdy - cdxbdy)                                       &
            + blift * (cdxady - adxcdy)                                       &
            + clift * (adxbdy - bdxady)

        inside = det > 0.0_qp
    end function incircle_robust

end module robust_predicates
