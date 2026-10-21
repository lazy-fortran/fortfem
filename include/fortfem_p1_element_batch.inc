        integer :: e, nel
        real(ck) :: a, b, c, d, det, scale

        nel = size(jac,1)
        if (size(jac,2) /= 4 .or. size(diffusion) /= nel .or. &
            size(reaction) /= nel .or. size(blocks,1) /= nel .or. &
            size(blocks,2) /= 3 .or. size(blocks,3) /= 3 .or. &
            size(status) /= nel) error stop 'p1_element_batch: shape mismatch'
        !$omp parallel do simd default(none) schedule(static) &
        !$omp shared(jac, diffusion, reaction, blocks, status, nel) &
        !$omp private(e, a, b, c, d, det, scale)
        do e = 1, nel
            status(e) = 1
            blocks(e,:,:) = 0
            if (.not. all(ieee_is_finite(jac(e,:)))) cycle
            if (.not. ieee_is_finite(diffusion(e))) cycle
            if (.not. ieee_is_finite(reaction(e))) cycle
            a = real(jac(e,1),ck); b = real(jac(e,2),ck)
            c = real(jac(e,3),ck); d = real(jac(e,4),ck)
            call measure(a,b,c,d,det,scale)
            status(e) = 2
            if (.not. ieee_is_finite(scale)) cycle
            if (abs(det) <= epsilon(1.0_ck)*scale) cycle
            call cell(a,b,c,d,real(diffusion(e),ck),real(reaction(e),ck), &
                blocks(e,1,1),blocks(e,2,1),blocks(e,3,1), &
                blocks(e,1,2),blocks(e,2,2),blocks(e,3,2), &
                blocks(e,1,3),blocks(e,2,3),blocks(e,3,3))
            status(e) = 0
            if (.not. all(ieee_is_finite(blocks(e,:,:)))) then
                status(e) = 3
                blocks(e,:,:) = 0
            end if
        end do
        !$omp end parallel do simd
