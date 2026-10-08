program test_affine_map_invariance
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_triangle_affine_map, only: invert_triangle_affine_map, &
        invert_triangle_affine_map_jvp, invert_triangle_affine_map_vjp
    use fortfem_tetra_affine_map, only: invert_tetra_affine_map, &
        invert_tetra_affine_map_jvp, invert_tetra_affine_map_vjp
    implicit none
    real(dp), allocatable :: v(:, :), vd(:, :), vb(:, :), point(:), pd(:), pb(:)
    real(dp), allocatable :: reference(:), expected(:), rd(:), rb(:), weights(:)
    integer :: dimension, sample, row, column, status, shift

    do dimension = 2, 3
        allocate(v(dimension, dimension + 1), vd(dimension, dimension + 1))
        allocate(vb(dimension, dimension + 1), weights(dimension + 1))
        allocate(point(dimension), pd(dimension), pb(dimension))
        allocate(reference(dimension), expected(dimension), rd(dimension), rb(dimension))
        do row = 1, dimension
            expected(row) = 0.05_dp*row
            rb(row) = (-1.0_dp)**row*(row + 0.5_dp)
        end do
        weights(1) = 1.0_dp - sum(expected)
        weights(2:) = expected
        do sample = 1, 9
            do row = 1, dimension
                v(row, 1) = 0.1_dp*row - 0.07_dp*sample
                do column = 1, dimension
                    v(row, column + 1) = v(row, 1) + 0.03_dp*(row - column)
                    if (row == column) v(row, column + 1) = &
                        v(row, column + 1) + 1.0_dp + 0.05_dp*sample
                end do
                do column = 1, dimension + 1
                    vd(row, column) = 0.02_dp*row - 0.015_dp*column
                end do
            end do
            point = matmul(v, weights)
            pd = matmul(vd, weights)
            if (dimension == 2) then
                call invert_triangle_affine_map(v, point, reference, status)
            else
                call invert_tetra_affine_map(v, point, reference, status)
            end if
            call check_condition(status == 0, "Affine inverse accepts oriented simplex")
            call check_condition(maxval(abs(reference - expected)) < 2.0e-14_dp, &
                "Affine inverse recovers prescribed barycentric coordinates")
            if (dimension == 2) then
                call invert_triangle_affine_map_jvp(v, point, vd, pd, rd, status)
            else
                call invert_tetra_affine_map_jvp(v, point, vd, pd, rd, status)
            end if
            call check_condition(status == 0, "Affine inverse tangent succeeds")
            call check_condition(maxval(abs(rd)) < 2.0e-14_dp, &
                "Joint vertex/point motion preserves barycentric coordinates")
            if (dimension == 2) then
                call invert_triangle_affine_map_vjp(v, point, rb, vb, pb, status)
            else
                call invert_tetra_affine_map_vjp(v, point, rb, vb, pb, status)
            end if
            call check_condition(status == 0, "Affine inverse reverse product succeeds")
            call check_condition(maxval(abs(sum(vb, dim=2) + pb)) < 2.0e-14_dp, &
                "Affine inverse reverse product preserves translation invariance")
            call check_condition(abs(sum(vb*vd) + dot_product(pb, pd)) < 2.0e-14_dp, &
                "Reverse product annihilates joint coordinate-preserving motion")
        end do
        ! Dyadic simplices and coordinates retain exact differences even far
        ! from the origin. Expanding their determinant in global coordinates
        ! destroys that property and is an inadmissible generated replacement.
        do shift = 20, 40, 10
            do row = 1, dimension
                v(row, 1) = (-1.0_dp)**row*2.0_dp**shift
            end do
            do column = 1, dimension
                v(:, column + 1) = v(:, 1)
                v(column, column + 1) = v(column, column + 1) + 2.0_dp
            end do
            point = v(:, 1) + 0.25_dp
            if (dimension == 2) then
                call invert_triangle_affine_map(v, point, reference, status)
            else
                call invert_tetra_affine_map(v, point, reference, status)
            end if
            call check_condition(status == 0, "Translated dyadic simplex is valid")
            call check_condition(maxval(abs(reference - 0.125_dp)) == 0.0_dp, &
                "Large translation preserves exactly representable local coordinates")
        end do
        deallocate(v, vd, vb, point, pd, pb, reference, expected, rd, rb, weights)
    end do
    call check_summary("Affine inverse independent invariance")
end program test_affine_map_invariance
