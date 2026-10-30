program test_elasticity_planar_acoustic_baffled_interface
    !! Finite elastic block whose top face is a finite interface in an infinite
    !! hard wall over a lossy fluid. The periodic DtN (lfac = 1) makes the two
    !! end samples of the interface neighbours through the fluid; the baffled
    !! DtN (dtn_period_factor = lfac) removes that coupling, and the solution
    !! converges as lfac grows because the periodic images decay through the
    !! fluid loss. The form itself is checked against the Rayleigh integral in
    !! test_planar_acoustic_baffled_dtn; here the solver must reproduce the
    !! periodic solution for lfac = 1, converge geometrically in lfac, and the
    !! periodic artefact must sit at the end samples.
    use check, only: check_condition, check_summary
    use fortfem_boundary, only: solve_elasticity_planar_acoustic_dtn_p1
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 32, ny = 8, m = nx + 1
    integer, parameter :: nv = (nx + 1)*(ny + 1)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: h = 1.0_dp/nx, period = m*h
    real(dp), parameter :: omega = 4.0_dp*pi, rho_f = 1.0_dp
    real(dp), parameter :: young = 40.0_dp, nu = 0.3_dp, rho_s = 2.0_dp
    complex(dp), parameter :: c = (1.0_dp, -0.05_dp)
    real(dp), parameter :: factors(6) = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, &
        16.0_dp, 32.0_dp]
    real(dp) :: vertices(2, nv), normals(2, m)
    integer :: triangles(3, 2*nx*ny), inodes(m), dofs(2*(nx + 1))
    complex(dp) :: values(2*(nx + 1)), load(2*nv), pinc(m)
    complex(dp) :: periodic(2*nv), solution(2*nv, 6), top(m, 6)
    real(dp) :: change(5), end_artefact, centre_artefact
    integer :: i, status
    logical :: all_passed

    all_passed = .true.
    call build_problem()
    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, c, rho_f, period, m/2, cmplx(young, 0.0_dp, dp), nu, &
        rho_s, load, pinc, dofs, values, periodic, status)
    call record(status == 0, "periodic interface solves")
    do i = 1, size(factors)
        call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, &
            inodes, normals, omega, c, rho_f, period, m/2, &
            cmplx(young, 0.0_dp, dp), nu, rho_s, load, pinc, dofs, values, &
            solution(:, i), status, dtn_period_factor=factors(i))
        call record(status == 0, "baffled interface solves")
        top(:, i) = solution(2*inodes, i)
    end do
    call record(maxval(abs(solution(:, 1) - periodic)) == 0.0_dp, &
        "lfac = 1 is the periodic interface")
    do i = 1, 5
        change(i) = maxval(abs(top(:, i + 1) - top(:, i)))/maxval(abs(top(:, 6)))
    end do
    write (*, "(a,5es10.3)") "relative change of u_n, lfac 1->2 ... 16->32: ", &
        change
    call record(all(change(2:5) < 0.3_dp*change(1:4)), &
        "baffled solution converges geometrically in lfac")
    call record(change(5) < 1.0e-4_dp, "lfac = 32 agrees with lfac = 16")
    end_artefact = max(abs(top(1, 1) - top(1, 6)), abs(top(m, 1) - top(m, 6)))
    centre_artefact = abs(top(m/2 + 1, 1) - top(m/2 + 1, 6))
    write (*, "(a,2es10.3)") "periodic - baffled u_n at the end, centre: ", &
        end_artefact/maxval(abs(top(:, 6))), centre_artefact/maxval(abs(top(:, 6)))
    call record(end_artefact > 3.0_dp*centre_artefact, &
        "the periodic artefact is largest at the end samples")
    call record(maxval(abs(top(:, 6) - top(m:1:-1, 6))) <= &
        1.0e-10_dp*maxval(abs(top(:, 6))), &
        "baffled solution of the symmetric block is mirror symmetric")

    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, c, rho_f, period, m/2, cmplx(young, 0.0_dp, dp), nu, &
        rho_s, load, pinc, dofs, values, periodic, status, &
        dtn_period_factor=1.5_dp)
    call record(status == 1, "lfac*m must be an integer")
    call solve_elasticity_planar_acoustic_dtn_p1(vertices, triangles, inodes, &
        normals, omega, c, rho_f, period, m/2, cmplx(young, 0.0_dp, dp), nu, &
        rho_s, load, pinc, dofs, values, periodic, status, &
        dtn_period_factor=0.5_dp)
    call record(status == 1, "lfac < 1 is rejected")

    call check_summary("Elastic block under a baffled planar DtN interface")
    if (.not. all_passed) error stop 1

contains

    pure integer function node(col, row)
        integer, intent(in) :: col, row

        node = 1 + col + (nx + 1)*row
    end function node

    subroutine build_problem()
        integer :: col, row, t, ll, k

        do row = 0, ny
            do col = 0, nx
                vertices(:, node(col, row)) = h*real([col, row], dp)
            end do
        end do
        t = 0
        do row = 0, ny - 1
            do col = 0, nx - 1
                ll = node(col, row)
                t = t + 1
                triangles(:, t) = [ll, ll + 1, ll + nx + 2]
                t = t + 1
                triangles(:, t) = [ll, ll + nx + 2, ll + nx + 1]
            end do
        end do
        ! Mirror-symmetric diagonals would need a symmetric mesh; use the
        ! union-jack pattern so the block is symmetric about x = 1/2.
        t = 0
        do row = 0, ny - 1
            do col = 0, nx - 1
                ll = node(col, row)
                t = t + 2
                if (col >= nx/2) then
                    triangles(:, t - 1) = [ll, ll + 1, ll + nx + 1]
                    triangles(:, t) = [ll + 1, ll + nx + 2, ll + nx + 1]
                end if
            end do
        end do
        do k = 1, m
            inodes(k) = node(k - 1, ny)
        end do
        normals(1, :) = 0.0_dp
        normals(2, :) = 1.0_dp
        pinc = (1.0_dp, 0.0_dp)
        do k = 0, nx
            dofs(2*k + 1) = 2*node(k, 0) - 1
            dofs(2*k + 2) = 2*node(k, 0)
        end do
        values = (0.0_dp, 0.0_dp)
        load = (0.0_dp, 0.0_dp)
    end subroutine build_problem

    subroutine record(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record

end program test_elasticity_planar_acoustic_baffled_interface
