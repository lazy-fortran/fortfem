program test_triangle_compat_golden
    !! Regression test for the Triangle-compatible mesher against a small
    !! golden mesh produced by the Triangle oracle ('triangle -pq20Y') and
    !! checked into test/data.  Runs without any external binary.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use triangle_compat, only: tc_result_t, triangulate_compat
    use tc_poly_io, only: read_poly
    use tc_compat_check, only: read_node_file, read_ele_file, meshes_match
    use fortfem_api, only: mesh_t, mesh_from_domain
    use check, only: check_condition, check_summary
    implicit none

    character(*), parameter :: data_dir = 'test/data/triangle_compat/'
    real(dp), allocatable :: points(:, :), holes(:, :), xy(:, :)
    integer, allocatable :: segments(:, :), segmarks(:), tris(:, :)
    type(tc_result_t) :: res
    type(mesh_t) :: mesh
    integer :: stat

    call read_poly(data_dir//'annulus.poly', points, segments, segmarks, &
                   holes, stat)
    call check_condition(stat == 0, 'read annulus.poly')

    call read_node_file(data_dir//'annulus.golden.node', xy, stat)
    call check_condition(stat == 0, 'read annulus.golden.node')
    call read_ele_file(data_dir//'annulus.golden.ele', tris, stat)
    call check_condition(stat == 0, 'read annulus.golden.ele')

    call triangulate_compat(points, segments, holes, res, stat, &
                            min_angle=20.0_dp, quality=.true., nobisect=1)
    call check_condition(stat == 0, 'triangulate_compat succeeds')

    call check_condition(meshes_match(res, xy, tris, 'annulus q20Y'), &
                         'mesh matches golden Triangle output')

    mesh = mesh_from_domain(points, segments, hole_points=holes, &
                            min_angle=20.0_dp, mode='triangle')
    call check_condition(mesh%data%n_vertices == size(xy, 2) .and. &
                         mesh%data%n_triangles == size(tris, 2), &
                         'mesh_from_domain mode=triangle matches golden counts')

    call check_summary('Triangle compat golden')
end program test_triangle_compat_golden
