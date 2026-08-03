module tc_poly_io
    !! Triangle .poly/.node/.ele file I/O for the Triangle-compatible mesher.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use triangle_compat, only: tc_result_t
    implicit none
    private

    public :: read_poly, write_node_ele

contains

    subroutine read_poly(filename, points, segments, segmentmarkers, holes, &
                         stat)
        character(*), intent(in) :: filename
        real(dp), allocatable, intent(out) :: points(:, :), holes(:, :)
        integer, allocatable, intent(out) :: segments(:, :), segmentmarkers(:)
        integer, intent(out) :: stat

        integer :: unit, ios, nv, dim_, nattr, nmark, ns, nh, i, idx
        integer :: firstnumber
        real(dp) :: x, y
        real(dp), allocatable :: attrs(:)
        character(1024) :: line

        stat = 0
        open (newunit=unit, file=filename, status='old', action='read', &
              iostat=ios)
        if (ios /= 0) then
            stat = 1
            return
        end if

        call next_data_line(unit, line, ios)
        read (line, *) nv, dim_, nattr, nmark
        allocate (points(2, nv), attrs(nattr))
        firstnumber = 1
        do i = 1, nv
            call next_data_line(unit, line, ios)
            if (nmark > 0) then
                read (line, *) idx, x, y
            else
                read (line, *) idx, x, y
            end if
            if (i == 1) firstnumber = idx
            points(:, idx - firstnumber + 1) = [x, y]
        end do

        call next_data_line(unit, line, ios)
        read (line, *) ns, nmark
        allocate (segments(2, ns), segmentmarkers(ns))
        segmentmarkers = 0
        do i = 1, ns
            call next_data_line(unit, line, ios)
            if (nmark > 0) then
                read (line, *, iostat=ios) idx, segments(1, i), &
                    segments(2, i), segmentmarkers(i)
                if (ios /= 0) then
                    read (line, *) idx, segments(1, i), segments(2, i)
                    segmentmarkers(i) = 0
                end if
            else
                read (line, *) idx, segments(1, i), segments(2, i)
            end if
        end do
        segments = segments - firstnumber + 1

        call next_data_line(unit, line, ios)
        if (ios /= 0) then
            allocate (holes(2, 0))
        else
            read (line, *) nh
            allocate (holes(2, nh))
            do i = 1, nh
                call next_data_line(unit, line, ios)
                read (line, *) idx, holes(1, i), holes(2, i)
            end do
        end if
        close (unit)
    end subroutine read_poly

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

    subroutine write_node_ele(basename, res)
        character(*), intent(in) :: basename
        type(tc_result_t), intent(in) :: res
        integer :: unit, i

        open (newunit=unit, file=trim(basename)//'.node', action='write')
        write (unit, '(i0, a)') res%npoints, ' 2 0 1'
        do i = 1, res%npoints
            write (unit, '(i0, 2(1x, es24.16e3), 1x, i0)') i, &
                res%points(1, i), res%points(2, i), res%pointmarkers(i)
        end do
        close (unit)

        open (newunit=unit, file=trim(basename)//'.ele', action='write')
        write (unit, '(i0, a)') res%ntriangles, ' 3 0'
        do i = 1, res%ntriangles
            write (unit, '(i0, 3(1x, i0))') i, res%triangles(1, i), &
                res%triangles(2, i), res%triangles(3, i)
        end do
        close (unit)
    end subroutine write_node_ele

end module tc_poly_io
