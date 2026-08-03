program triangle_compat_cli
    !! Mesh a Triangle .poly file with the Triangle-compatible pipeline.
    !! Usage: triangle_compat_cli <input.poly> <output_basename> [min_angle]
    !! Writes <output_basename>.node and <output_basename>.ele.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use triangle_compat, only: tc_result_t, triangulate_compat
    use tc_poly_io, only: read_poly, write_node_ele
    implicit none

    character(1024) :: infile, outbase, arg
    real(dp), allocatable :: points(:, :), holes(:, :)
    integer, allocatable :: segments(:, :), segmentmarkers(:)
    type(tc_result_t) :: res
    real(dp) :: min_angle
    integer :: stat, nargs
    real(dp) :: t0, t1

    nargs = command_argument_count()
    if (nargs < 2) then
        write (error_unit, '(a)') &
            'usage: triangle_compat_cli <input.poly> <output_basename> [min_angle]'
        stop 2
    end if
    call get_command_argument(1, infile)
    call get_command_argument(2, outbase)
    min_angle = 20.0_dp
    if (nargs >= 3) then
        call get_command_argument(3, arg)
        read (arg, *) min_angle
    end if

    call read_poly(trim(infile), points, segments, segmentmarkers, holes, stat)
    if (stat /= 0) then
        write (error_unit, '(a)') 'error reading '//trim(infile)
        stop 2
    end if

    call cpu_time(t0)
    call triangulate_compat(points, segments, holes, res, stat, &
                            min_angle=min_angle, quality=min_angle > 0.0_dp, &
                            nobisect=1, segmentmarkers=segmentmarkers)
    call cpu_time(t1)
    if (stat /= 0) then
        write (error_unit, '(a, i0)') 'triangulate_compat failed, stat=', stat
        stop 1
    end if

    call write_node_ele(trim(outbase), res)
    write (*, '(a, i0, a, i0, a, f8.3, a, f10.4, a)') 'vertices=', res%npoints, &
        ' triangles=', res%ntriangles, ' min_angle=', res%min_angle_deg, &
        ' cpu=', t1 - t0, ' s'
end program triangle_compat_cli
