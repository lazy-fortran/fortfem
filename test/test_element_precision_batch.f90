program test_element_precision_batch
    use, intrinsic :: iso_fortran_env, only: real32, real64
    use check, only: check_condition, check_summary
    use fortfem_element_batch, only: p1_element_batch
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
                                      csr_assemble_batched
!$  use omp_lib, only: omp_set_num_threads
    use fortfem_generated_p2_jet_fp64, only: generated_p2_jet_fp64
    use fortfem_generated_p2_jet_fp32, only: generated_p2_jet_fp32
    use fortfem_generated_triangle_factors_fp64, only: &
        generated_triangle_factors_fp64
    use fortfem_generated_triangle_factors_fp32, only: &
        generated_triangle_factors_fp32
    implicit none
    integer, parameter :: sp = real32, dp = real64, nel = 257
    real(dp) :: jac(nel, 4), diffusion(nel), reaction(nel), blocks(nel, 3, 3)
    real(dp) :: widened(nel, 3, 3), reference(3, 3), gradients(2, 3), det, area
    real(dp) :: kappa, bound, u, gamma, dense(4, 4), error
    real(sp) :: jac32(nel, 4), d32(nel), r32(nel), blocks32(nel, 3, 3)
    type(csr_matrix_t) :: matrix
    real(dp) :: jet(6, 6), nodal(6), x, y, f(5), expected(5), v(2)
    real(sp) :: jet32(6, 6), f32(5)
    logical :: cell_ok(3)
    integer :: status(nel), e, i, j, p, conn(3, 2), thread

    do e = 1, nel
        jac(e, :) = [1.0_dp + real(e, dp)/nel, 0.125_dp, &
                     -0.25_dp, 0.75_dp + real(e, dp)/(2*nel)]
        if (mod(e, 2) == 0) jac(e, 3:4) = -jac(e, 3:4)
        diffusion(e) = 0.5_dp + real(e, dp)/nel
        reaction(e) = 0.2_dp
    end do
    jac32 = real(jac, sp); d32 = real(diffusion, sp); r32 = real(reaction, sp)
    do thread = 1, 4, 3
!$      call omp_set_num_threads(thread)
        call p1_element_batch(jac, diffusion, reaction, blocks, status)
        call check_condition(all(status == 0), 'fp64 regular geometry')
        call p1_element_batch(jac32, d32, r32, blocks32, status)
        call check_condition(all(status == 0), 'fp32 regular geometry')
        call p1_element_batch(jac32, d32, r32, widened, status)
        call check_condition(all(status == 0), 'fp32 storage / fp64 compute')
        cell_ok = .true.
        do e = 1, nel
            det = jac(e, 1)*jac(e, 4) - jac(e, 2)*jac(e, 3)
            area = abs(det)/2
            ! Independent physical barycentric-gradient oracle.
            gradients(:, 2) = [jac(e, 4), -jac(e, 3)]/det
            gradients(:, 3) = [-jac(e, 2), jac(e, 1)]/det
            gradients(:, 1) = -sum(gradients(:, 2:3), dim=2)
            do j = 1, 3
                do i = 1, 3
                    reference(i, j) = diffusion(e)*area* &
                                      dot_product(gradients(:, i), gradients(:, j)) + &
                                      reaction(e)*area/12
                    if (i == j) reference(i, j) = reference(i, j) + &
                                                  reaction(e)*area/12
                end do
            end do
            kappa = sum(abs(jac(e, :)))**2/abs(det)
            u = epsilon(1.0_dp)/2
            gamma = 96*u/(1 - 96*u)
            bound = gamma*kappa**2*maxval(abs(reference))
            cell_ok(1) = cell_ok(1) .and. &
                         maxval(abs(blocks(e, :, :) - reference)) <= bound
            u = real(epsilon(1.0_sp), dp)/2
            gamma = 96*u/(1 - 96*u)
            bound = gamma*kappa**2*maxval(abs(reference))
            cell_ok(2) = cell_ok(2) .and. &
                         maxval(abs(real(blocks32(e, :, :), dp) - reference)) <= bound
            cell_ok(3) = cell_ok(3) .and. &
                         maxval(abs(widened(e, :, :) - reference)) <= bound
        end do
        call check_condition(cell_ok(1), 'fp64 independent cell oracle')
        call check_condition(cell_ok(2), 'fp32 cell roundoff model')
        call check_condition(cell_ok(3), 'storage conversion in cell model')
        conn(:, 1) = [1, 2, 3]; conn(:, 2) = [2, 4, 3]
        call csr_build_pattern(4, conn, .false., matrix)
        call csr_assemble_batched(matrix, conn, blocks32(:2, :, :))
        dense = 0
        do e = 1, 2
            do j = 1, 3
                do i = 1, 3
                    dense(conn(i, e), conn(j, e)) = &
                        dense(conn(i, e), conn(j, e)) + real(blocks32(e, i, j), dp)
                end do
            end do
        end do
        error = 0
        do i = 1, 4
            do p = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                error = max(error, abs(matrix%val(p) - dense(i, matrix%colidx(p))))
            end do
        end do
        call check_condition(error == 0, 'shared-node fp64 CSR accumulation')
        call csr_assemble_batched(matrix, conn, blocks(:2, :, :))
        dense = 0
        do e = 1, 2
            do j = 1, 3
                do i = 1, 3
                    dense(conn(i, e), conn(j, e)) = &
                        dense(conn(i, e), conn(j, e)) + blocks(e, i, j)
                end do
            end do
        end do
        error = 0
        do i = 1, 4
            do p = matrix%rowptr(i), matrix%rowptr(i + 1) - 1
                error = max(error, abs(matrix%val(p) - dense(i, matrix%colidx(p))))
            end do
        end do
        call check_condition(error == 0, 'fp64 batch CSR dense oracle')
    end do
    jac(1, :) = 0
    call p1_element_batch(jac, diffusion, reaction, blocks, status)
    call check_condition(status(1) /= 0, 'singular cell rejected before divide')
    call check_condition(all(blocks(1, :, :) == 0), 'invalid cell zeroed')
    call p1_element_batch(jac(:0, :), diffusion(:0), reaction(:0), &
                          blocks(:0, :, :), status(:0))
    nodal = [0.0_dp, 1.0_dp, 2.0_dp, 0.25_dp, 1.0_dp, 0.5_dp]
    do i = 1, 13
        x = real(i, dp)/31; y = real(13 - i, dp)/29
        call generated_p2_jet_fp64(x, y, jet)
        call generated_p2_jet_fp32(real(x, sp), real(y, sp), jet32)
        expected = [x*x + x*y + 2*y*y, 2*x + y, x + 4*y, 2.0_dp, 1.0_dp]
        call check_condition(maxval(abs(matmul(jet(:5, :), nodal) - expected)) &
                          < 32*epsilon(1.0_dp), 'P2 polynomial jets independent oracle')
        call check_condition(abs(dot_product(jet(6, :), nodal) - 4) < &
                             32*epsilon(1.0_dp), 'P2 eta Hessian')
        call check_condition(maxval(abs(real(jet32, dp) - jet)) < &
                             32*real(epsilon(1.0_sp), dp), 'fp32 basis derivatives')
    end do
    v = [0.3_dp, -0.7_dp]
    call generated_triangle_factors_fp64(2.0_dp, 0.25_dp, -0.5_dp, 1.0_dp, &
                                         v(1), v(2), f(1), f(2), f(3), f(4), f(5))
    expected(1) = 2.125_dp
    expected(2:3) = [(v(1) - 0.25_dp*v(2)), (0.5_dp*v(1) + 2*v(2))]/expected(1)
    expected(4:5) = [(2*v(1) - 0.5_dp*v(2)), (0.25_dp*v(1) + v(2))]/expected(1)
    call check_condition(maxval(abs(f - expected)) < 8*epsilon(1.0_dp), &
                         'covariant and contravariant Piola independent oracle')
    call generated_triangle_factors_fp32(2.0_sp, 0.25_sp, -0.5_sp, 1.0_sp, &
                 real(v(1), sp), real(v(2), sp), f32(1), f32(2), f32(3), f32(4), f32(5))
    call check_condition(maxval(abs(real(f32, dp) - expected)) < &
                         8*real(epsilon(1.0_sp), dp), 'fp32 Piola factors')
    call check_summary('element precision batch')
end program test_element_precision_batch
