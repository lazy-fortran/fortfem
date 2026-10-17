program test_meridian_hermitian_forms
    use fortfem_kinds, only: dp
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build, &
        meridian_stiffness
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t
    use meridian_test_meshes, only: make_warped_rectangle_mesh
    use fortsparse, only: csc_t, csc_z_t, csc_from_triplet, &
        fortsparse_status_t, sparse_solver_t, sparse_factor, sparse_solve, &
        sparse_destroy, FORTSPARSE_BACKEND_NATIVE, FORTSPARSE_SYM_SPD, &
        FORTSPARSE_SYM_HERMITIAN
    implicit none

    call curved_mesh_forms()
    call synthetic_gram()
    print '(a)', 'Canonical Hermitian meridian forms and native SPD/HPD: PASS'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 1
        end if
    end subroutine require

    function dense(a) result(values)
        type(csc_z_t), intent(in) :: a
        complex(dp) :: values(a%nrow, a%ncol)
        integer :: j, p
        values = (0.0_dp, 0.0_dp)
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                values(a%row_idx(p), j) = a%val(p)
            end do
        end do
    end function dense

    function real_dense(a) result(values)
        type(csc_t), intent(in) :: a
        real(dp) :: values(a%nrow, a%ncol)
        integer :: j, p
        values = 0.0_dp
        do j = 1, a%ncol
            do p = a%col_ptr(j), a%col_ptr(j + 1) - 1
                values(a%row_idx(p), j) = a%val(p)
            end do
        end do
    end function real_dense

    subroutine synthetic_gram()
        ! Independent dense products check the values, not just mirrored storage.
        ! A full-rank synthetic differential avoids the physical gradient kernel.
        real(dp), parameter :: ar(4, 3) = reshape([ &
            1.0_dp, -0.7_dp, 0.3_dp, 0.0_dp, &
            0.5_dp, 1.0_dp, -0.6_dp, 0.9_dp, &
            0.0_dp, 0.2_dp, 1.0_dp, -0.4_dp], [4, 3])
        real(dp), parameter :: ai(4, 3) = reshape([ &
            0.2_dp, 0.0_dp, -0.3_dp, 0.4_dp, &
            0.0_dp, -0.1_dp, 0.7_dp, 0.2_dp, &
            -0.5_dp, 0.3_dp, 0.0_dp, 0.6_dp], [4, 3])
        real(dp), parameter :: m(4, 4) = reshape([ &
            2.1_dp, 0.2_dp, -0.3_dp, 0.1_dp, &
            0.2_dp, 3.7_dp, 0.4_dp, -0.2_dp, &
            -0.3_dp, 0.4_dp, 1.2_dp, 0.3_dp, &
            0.1_dp, -0.2_dp, 0.3_dp, 4.9_dp], [4, 4])
        integer, parameter :: modes(3) = [0, 1, -3]
        integer :: rows(12), columns(12), mr(16), mc(16)
        integer :: i, j, q, space, n, info
        type(meridian_feec_t) :: fs
        type(csc_z_t) :: k
        type(csc_t) :: kr
        type(fortsparse_status_t) :: status
        type(sparse_solver_t) :: solver
        complex(dp) :: z(4, 3), expected(3, 3), actual(3, 3)
        complex(dp) :: x(3), rhs(3), solution(3)
        real(dp) :: xr(3), br(3), solved(3)

        q = 0
        do j = 1, 3
            do i = 1, 4
                q = q + 1
                rows(q) = i
                columns(q) = j
            end do
        end do
        q = 0
        do j = 1, 4
            do i = 1, 4
                q = q + 1
                mr(q) = i
                mc(q) = j
            end do
        end do
        call csc_from_triplet(4, 3, rows, columns, reshape(ar, [12]), &
            fs%grad%a, status)
        call csc_from_triplet(4, 3, rows, columns, reshape(ai, [12]), &
            fs%grad%b, status)
        fs%curl = fs%grad
        fs%div = fs%grad
        do space = 1, 3
            call csc_from_triplet(4, 4, mr, mc, reshape(m, [16]), &
                fs%mass(space), status)
        end do
        x = cmplx([1.0_dp, -2.0_dp, 0.5_dp], &
            [0.2_dp, 0.7_dp, -0.3_dp], dp)
        do space = 0, 2
            do q = 1, size(modes)
                n = modes(q)
                call meridian_stiffness(fs, space, n, k, info)
                call require(info == 0, 'Gram assembly succeeds')
                actual = dense(k)
                z = cmplx(ar, real(n, dp)*ai, dp)
                expected = matmul(conjg(transpose(z)), matmul(m, z))
                call require(maxval(abs(actual - expected)) &
                    < 128*epsilon(1.0_dp)*maxval(abs(expected)), &
                    'Sparse Gram agrees with independently ordered dense product')
                call require(all(actual == conjg(transpose(actual))), &
                    'Gram is exactly Hermitian, including real diagonal')
                solver%backend_id = FORTSPARSE_BACKEND_NATIVE
                solver%options%symmetry = FORTSPARSE_SYM_HERMITIAN
                call sparse_factor(solver, k, status)
                call require(status%code == 0, 'Native HPD admits canonical Gram')
                rhs = matmul(actual, x)
                call sparse_solve(solver, rhs, solution, status)
                call require(status%code == 0, 'Native HPD solve succeeds')
                call require(maxval(abs(solution - x)) < 1.0e-12_dp, &
                    'Native HPD recovers independent manufactured vector')
                call sparse_destroy(solver)
                if (n /= 0) cycle
                kr%nrow = k%nrow
                kr%ncol = k%ncol
                kr%nnz = k%nnz
                kr%col_ptr = k%col_ptr
                kr%row_idx = k%row_idx
                kr%val = real(k%val, dp)
                solver%options%symmetry = FORTSPARSE_SYM_SPD
                call sparse_factor(solver, kr, status)
                call require(status%code == 0, 'Native real SPD admits n=0 Gram')
                xr = real(x, dp)
                br = matmul(real(actual, dp), xr)
                call sparse_solve(solver, br, solved, status)
                call require(status%code == 0, 'Native real SPD solve succeeds')
                call require(maxval(abs(solved - xr)) < 1.0e-12_dp, &
                    'Native real SPD recovers independent manufactured vector')
                call sparse_destroy(solver)
            end do
        end do
    end subroutine synthetic_gram

    subroutine curved_mesh_forms()
        ! Shared global dofs receive duplicate element contributions. Compression
        ! must not sum the two transpose entries in unrelated floating orders.
        type(meridian_mesh_t) :: mesh
        type(meridian_feec_t) :: fs
        type(csc_z_t) :: k
        real(dp), allocatable :: mass(:, :)
        complex(dp), allocatable :: stiffness(:, :)
        integer :: space, mode, info

        call make_warped_rectangle_mesh(2, 2, mesh)
        call meridian_feec_build(mesh, 2, fs, info)
        call require(info == 0, 'Curved quadratic shared-dof fixture builds')
        do space = 0, 3
            mass = real_dense(fs%mass(space))
            call require(all(mass == transpose(mass)), &
                'Compressed physical mass is exactly symmetric')
        end do
        do space = 0, 2
            do mode = 0, 2, 2
                call meridian_stiffness(fs, space, mode, k, info)
                call require(info == 0, 'Curved Gram assembly succeeds')
                stiffness = dense(k)
                call require(all(stiffness == conjg(transpose(stiffness))), &
                    'Curved shared-dof Gram is exactly Hermitian')
            end do
        end do
    end subroutine curved_mesh_forms
end program test_meridian_hermitian_forms
