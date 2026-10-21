program test_parametric_laplace_bem_topology
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: &
        assemble_parametric_laplace_bem, bem_operators_t
    use parametric_bem_test_geometry, only: test_sphere_t
    use fortfem_sphere_surface_mesh, only: generate_sphere_surface_mesh
    implicit none
    type(test_sphere_t) :: sphere
    type(bem_operators_t) :: op
    integer, allocatable :: bad(:, :)
    integer :: status, swap
    logical :: passed
    passed = .true.
    call generate_sphere_surface_mesh(1.0_dp, 0, sphere%vertices, sphere%faces)
    bad = sphere%faces(:, 2:)
    call assemble_parametric_laplace_bem(sphere, bad, 4, op, status)
    call record(status /= 0, 'open shell is rejected')
    bad = sphere%faces
    swap = bad(1,1)
    bad(1,1) = bad(2,1)
    bad(2,1) = swap
    call assemble_parametric_laplace_bem(sphere, bad, 4, op, status)
    call record(status /= 0, 'inconsistent edge orientation is rejected')
    call assemble_parametric_laplace_bem(sphere, sphere%faces, 4, op, status)
    call record(status == 0, 'closed sphere is accepted')
    call record(op%genus == 0, 'sphere genus is zero')
    call check_summary('Parametric BEM topology')
    if (.not. passed) error stop 1
contains
    subroutine record(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed = passed .and. ok
        call check_condition(ok, label)
    end subroutine
end program
