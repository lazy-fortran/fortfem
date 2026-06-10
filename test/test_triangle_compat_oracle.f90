program test_triangle_compat_oracle
    !! Compare the Triangle-compatible mesher against a live Triangle
    !! oracle binary on the annulus PSLG, for both plain CDT ('-pY') and
    !! quality refinement ('-pq20Y').  Skips when no oracle is available;
    !! set TRIANGLE_ORACLE to point at a triangle binary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use triangle_compat, only: tc_result_t, triangulate_compat
    use tc_poly_io, only: read_poly
    use tc_compat_check, only: read_node_file, read_ele_file, meshes_match
    use check, only: check_condition, check_summary
    implicit none

    character(*), parameter :: poly_src = 'test/data/triangle_compat/annulus.poly'
    character(*), parameter :: work = '/tmp/fortfem_tc_oracle'
    character(1024) :: oracle
    logical :: have_oracle
    integer :: stat

    call get_environment_variable('TRIANGLE_ORACLE', oracle, status=stat)
    if (stat /= 0 .or. len_trim(oracle) == 0) then
        oracle = '/tmp/triangle-oracle/triangle'
    end if
    inquire (file=trim(oracle), exist=have_oracle)

    if (.not. have_oracle) then
        print '(a)', 'SKIP: no Triangle oracle binary at '//trim(oracle)
        print '(a)', '      (set TRIANGLE_ORACLE to enable this test)'
        call check_summary('Triangle compat oracle (skipped)')
        stop
    end if

    call run_case('-pq20YQ', 'annulus q20Y')
    call run_case('-pYQ', 'annulus CDT')

    call check_summary('Triangle compat oracle')

contains

    subroutine run_case(flags, label)
        character(*), intent(in) :: flags, label

        real(dp), allocatable :: points(:, :), holes(:, :), xy(:, :)
        integer, allocatable :: segments(:, :), segmarks(:), tris(:, :)
        type(tc_result_t) :: res
        real(dp) :: min_angle
        integer :: cstat, estat

        call execute_command_line('mkdir -p '//work// &
                                  ' && cp '//poly_src//' '//work//'/in.poly', &
                                  exitstat=estat)
        call check_condition(estat == 0, label//': stage input')

        call execute_command_line(trim(oracle)//' '//flags//' '// &
                                  work//'/in.poly', &
                                  cmdstat=cstat, exitstat=estat)
        call check_condition(cstat == 0 .and. estat == 0, &
                             label//': oracle run')

        call read_node_file(work//'/in.1.node', xy, estat)
        call check_condition(estat == 0, label//': read oracle nodes')
        call read_ele_file(work//'/in.1.ele', tris, estat)
        call check_condition(estat == 0, label//': read oracle elements')

        call read_poly(poly_src, points, segments, segmarks, holes, estat)
        call check_condition(estat == 0, label//': read PSLG')

        min_angle = 0.0_dp
        if (index(flags, 'q') > 0) min_angle = 20.0_dp
        call triangulate_compat(points, segments, holes, res, estat, &
                                min_angle=min_angle, &
                                quality=min_angle > 0.0_dp, nobisect=1)
        call check_condition(estat == 0, label//': triangulate_compat')

        call check_condition(meshes_match(res, xy, tris, label), &
                             label//': mesh matches oracle')
    end subroutine run_case

end program test_triangle_compat_oracle
