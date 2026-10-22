program test_pcg_scale_invariance
    use fortfem_kinds, only: dp
    use fortfem_advanced_solvers, only: solve_sparse, solver_options_t, solver_stats_t
    use fortfem_sparse_matrix, only: sparse_matrix_t, sparse_from_dense
    use check, only: check_condition, check_summary
    implicit none
    type(sparse_matrix_t) :: matrix
    type(solver_options_t) :: options
    type(solver_stats_t) :: stats
    real(dp) :: a(3,3), exact(3), rhs(3), solution(3), scale
    integer :: k
    a = reshape([4._dp,1._dp,0._dp,1._dp,3._dp,1._dp,0._dp,1._dp,2._dp], [3,3])
    call sparse_from_dense(a,matrix)
    options%method = 'pcg'
    options%preconditioner = 'ichol_row'
    options%fill_level = 8
    options%tolerance = 1.e-12_dp
    do k = 0,1
        scale = 10._dp**(-8*k)
        exact = scale*[1._dp,-2._dp,3._dp]
        rhs = matmul(a,exact)
        solution = 0
        call solve_sparse(matrix,rhs,solution,options,stats)
        call check_condition(stats%converged,'PCG converges for both RHS scales')
        call check_condition(sqrt(sum((solution-exact)**2))/scale < 1.e-11_dp, &
            'PCG preserves solution accuracy under RHS scaling')
    end do
    call check_summary('PCG scale invariance')
end program
