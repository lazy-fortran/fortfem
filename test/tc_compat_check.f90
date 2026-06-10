module tc_compat_check
    !! Helpers shared by the Triangle-compatibility tests: readers for
    !! Triangle .node/.ele files and mesh comparison against a tc_result_t.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use triangle_compat, only: tc_result_t
    implicit none
    private

    public :: read_node_file, read_ele_file, meshes_match

contains

    subroutine read_node_file(path, xy, stat)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: xy(:, :)
        integer, intent(out) :: stat

        integer :: unit, ios, n, i, idx, base
        character(1024) :: line

        stat = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        call next_data_line(unit, line, ios)
        if (ios /= 0) return
        read (line, *) n
        allocate (xy(2, n))
        base = 1
        do i = 1, n
            call next_data_line(unit, line, ios)
            if (ios /= 0) return
            read (line, *) idx, xy(1, i), xy(2, i)
            if (i == 1) base = idx
            if (idx - base + 1 /= i) return
        end do
        close (unit)
        stat = 0
    end subroutine read_node_file

    subroutine read_ele_file(path, tris, stat)
        character(*), intent(in) :: path
        integer, allocatable, intent(out) :: tris(:, :)
        integer, intent(out) :: stat

        integer :: unit, ios, n, i, idx, base
        character(1024) :: line

        stat = 1
        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        call next_data_line(unit, line, ios)
        if (ios /= 0) return
        read (line, *) n
        allocate (tris(3, n))
        base = huge(1)
        do i = 1, n
            call next_data_line(unit, line, ios)
            if (ios /= 0) return
            read (line, *) idx, tris(1, i), tris(2, i), tris(3, i)
            base = min(base, minval(tris(:, i)))
        end do
        close (unit)
        tris = tris - base + 1
        stat = 0
    end subroutine read_ele_file

    subroutine next_data_line(unit, line, ios)
        integer, intent(in) :: unit
        character(*), intent(out) :: line
        integer, intent(out) :: ios
        character(1024) :: buf
        integer :: k

        do
            read (unit, '(a)', iostat=ios) buf
            if (ios /= 0) then
                line = ''
                return
            end if
            k = index(buf, '#')
            if (k > 0) buf = buf(1:k - 1)
            if (len_trim(buf) > 0) then
                line = buf
                return
            end if
        end do
    end subroutine next_data_line

    function meshes_match(res, xy, tris, label) result(ok)
        !! Vertices must agree to 1e-12 in direct index order; triangle
        !! connectivity must be identical as a set of vertex triples.
        type(tc_result_t), intent(in) :: res
        real(dp), intent(in) :: xy(:, :)
        integer, intent(in) :: tris(:, :)
        character(*), intent(in) :: label
        logical :: ok

        real(dp), parameter :: tol = 1.0e-12_dp
        real(dp) :: maxdiff
        integer :: i

        ok = .false.
        if (res%npoints /= size(xy, 2)) then
            print '(a, a, i0, a, i0)', label, ': vertex count ', &
                res%npoints, ' /= ', size(xy, 2)
            return
        end if
        if (res%ntriangles /= size(tris, 2)) then
            print '(a, a, i0, a, i0)', label, ': triangle count ', &
                res%ntriangles, ' /= ', size(tris, 2)
            return
        end if
        maxdiff = maxval(abs(res%points - xy))
        if (maxdiff > tol) then
            print '(a, a, es12.3)', label, ': max vertex diff ', maxdiff
            return
        end if
        do i = 1, size(tris, 2)
            if (.not. has_triangle(res%triangles, tris(:, i))) then
                print '(a, a, 3(1x, i0))', label, &
                    ': golden triangle missing:', tris(:, i)
                return
            end if
        end do
        ok = .true.
    end function meshes_match

    function has_triangle(list, tri) result(found)
        integer, intent(in) :: list(:, :), tri(3)
        logical :: found
        integer :: a(3), b(3), i

        a = sorted3(tri)
        found = .false.
        do i = 1, size(list, 2)
            b = sorted3(list(:, i))
            if (all(a == b)) then
                found = .true.
                return
            end if
        end do
    end function has_triangle

    pure function sorted3(t) result(r)
        integer, intent(in) :: t(3)
        integer :: r(3)
        r = [minval(t), sum(t) - minval(t) - maxval(t), maxval(t)]
    end function sorted3

end module tc_compat_check
