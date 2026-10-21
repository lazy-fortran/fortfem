module fortfem_harmonic_field_energy
    !! Physical harmonic representatives are supplied by the caller. Their
    !! number is the volume Betti number g, not the surface Betti number 2g.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    use fortfem_generated_field_inner_product, only: field_inner_product
    implicit none
    private
    public :: assemble_field_gram, couple_harmonic_sector
contains
    subroutine assemble_field_gram(weights, fields, gram, status)
        !! Batch storage: Cartesian component, quadrature point, field.
        real(dp), intent(in) :: weights(:), fields(:, :, :)
        real(dp), intent(out) :: gram(:, :)
        integer, intent(out) :: status
        real(dp) :: value, total
        integer :: n, i, j, q

        status = 1
        gram = 0.0_dp
        n = size(fields, 3)
        if (n < 1 .or. size(weights) < 1) return
        if (size(fields, 1) /= 3 .or. size(fields, 2) /= size(weights)) return
        if (size(gram, 1) /= n .or. size(gram, 2) /= n) return
        if (.not. all(ieee_is_finite(weights))) return
        if (.not. all(ieee_is_finite(fields))) return
        if (any(weights < 0.0_dp) .or. .not. any(weights > 0.0_dp)) return
        do j = 1, n
            do i = 1, j
                total = 0.0_dp
                !$omp parallel do default(none) shared(fields, weights, i, j) &
                !$omp private(value) reduction(+:total) schedule(static)
                do q = 1, size(weights)
                    call field_inner_product(weights(q), fields(1, q, i), &
                        fields(2, q, i), fields(3, q, i), fields(1, q, j), &
                        fields(2, q, j), fields(3, q, j), value)
                    total = total+value
                end do
                !$omp end parallel do
                gram(i, j) = total
                gram(j, i) = total
            end do
        end do
        if (.not. all(ieee_is_finite(gram))) then
            gram = 0.0_dp
            return
        end if
        status = 0
    end subroutine

    subroutine couple_harmonic_sector(scalar, cross, harmonic, coordinates, &
            prescribed, u, coefficients, load, energy, status)
        !! E = u^T S u + 2 u^T C c + c^T G c; coordinates*c = prescribed.
        !! load is half dE/du at fixed prescribed physical coordinates.
        !! S,C,G must describe the same domain and field representatives.
        real(dp), intent(in) :: scalar(:, :), cross(:, :), harmonic(:, :)
        real(dp), intent(in) :: coordinates(:, :), prescribed(:), u(:)
        real(dp), intent(out) :: coefficients(:), load(:), energy
        integer, intent(out) :: status
        real(dp), allocatable :: c(:), su(:), cc(:), gc(:)
        real(dp) :: scalar_energy, cross_energy, sector_energy, total, tolerance
        integer :: n, g, info, i

        status = 1
        coefficients = 0.0_dp
        load = 0.0_dp
        energy = 0.0_dp
        n = size(u)
        g = size(prescribed)
        if (n < 1 .or. g < 1) return
        if (size(scalar, 1) /= n .or. size(scalar, 2) /= n) return
        if (size(cross, 1) /= n .or. size(cross, 2) /= g) return
        if (size(harmonic, 1) /= g .or. size(harmonic, 2) /= g) return
        if (size(coordinates, 1) /= g .or. size(coordinates, 2) /= g) return
        if (size(coefficients) /= g .or. size(load) /= n) return
        if (.not. all(ieee_is_finite(scalar))) return
        if (.not. all(ieee_is_finite(cross))) return
        if (.not. all(ieee_is_finite(harmonic))) return
        if (.not. all(ieee_is_finite(coordinates))) return
        if (.not. all(ieee_is_finite(prescribed))) return
        if (.not. all(ieee_is_finite(u))) return
        tolerance = 64*epsilon(1.0_dp)*real(n+g, dp)
        if (maxval(abs(scalar-transpose(scalar))) > &
            tolerance*maxval(abs(scalar))) return
        if (maxval(abs(harmonic-transpose(harmonic))) > &
            tolerance*maxval(abs(harmonic))) return
        do i = 1, n
            if (scalar(i, i) < 0.0_dp) return
        end do
        do i = 1, g
            if (harmonic(i, i) <= 0.0_dp) return
        end do
        allocate(c(g), su(n), cc(n), gc(g))
        call dense_solve(coordinates, prescribed, c, info)
        status = 2
        if (info /= 0) return
        if (.not. all(ieee_is_finite(c))) return
        su = matmul(scalar, u)
        cc = matmul(cross, c)
        gc = matmul(harmonic, c)
        scalar_energy = dot_product(u, su)
        cross_energy = dot_product(u, cc)
        sector_energy = dot_product(c, gc)
        total = scalar_energy+2*cross_energy+sector_energy
        if (.not. ieee_is_finite(total)) return
        tolerance = tolerance*(abs(scalar_energy)+2*abs(cross_energy)+ &
            abs(sector_energy))
        if (total < -tolerance) return
        if (.not. all(ieee_is_finite(su+cc))) return
        coefficients = c
        load = su+cc
        energy = max(0.0_dp, total)
        status = 0
    end subroutine
end module
