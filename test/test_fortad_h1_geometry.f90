program test_fortad_h1_geometry
    !! fortad's H1 geometry products against the fortsym-generated ones.
    !!
    !! fortfem's production kernel for this element-matrix entry is generated
    !! symbolically. This checks that the fortad path is an alternative and not
    !! a variant: for the same primal, at the same point, both products must
    !! agree to rounding.
    !!
    !! Agreement with fortsym alone would pass if the primal under
    !! tools/fortad/kernels stated the wrong computation - both fortad products
    !! would then be consistently wrong together. So the value is also checked
    !! against central differences of the primal, and the two products against
    !! each other through the adjoint identity.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_generated_bspline_h1_geometry_jvp, only: &
        generated_bspline_h1_geometry_jvp
    use fortfem_generated_bspline_h1_geometry_vjp, only: &
        generated_bspline_h1_geometry_vjp
    use fortfem_fortad_h1_geometry_jvp, only: fortfem_h1_geometry_jvp_fortad
    use fortfem_fortad_h1_geometry_vjp, only: fortfem_h1_geometry_vjp_fortad
    implicit none

    ! A well-conditioned cell: the Jacobian is far from singular, so the
    ! division by the determinant is not what the comparison measures.
    real(dp), parameter :: J11 = 1.3_dp, J21 = -0.4_dp, J12 = 0.2_dp, J22 = 1.1_dp
    real(dp), parameter :: GR1 = 0.7_dp, GR2 = -1.2_dp
    real(dp), parameter :: GC1 = -0.3_dp, GC2 = 0.9_dp
    real(dp), parameter :: BR = 0.45_dp, BC = 0.62_dp
    real(dp), parameter :: KS = 1.7_dp, KM = 0.8_dp, QW = 0.25_dp
    real(dp), parameter :: V(4) = [0.31_dp, -0.77_dp, 1.05_dp, 0.42_dp]

    integer :: failures

    failures = 0
    call check_products(failures)
    call check_differences(failures)

    if (failures == 0) then
        print *, "test_fortad_h1_geometry: all cases passed"
    else
        print *, "test_fortad_h1_geometry: ", failures, " case(s) FAILED"
        error stop 1
    end if

contains

    subroutine check_products(failures)
        integer, intent(inout) :: failures
        real(dp) :: tangent_ref, tangent_new
        real(dp) :: g_ref(4), g_new(4)
        real(dp), parameter :: u = 1.4_dp

        call generated_bspline_h1_geometry_jvp(J11, J21, J12, J22, GR1, GR2, &
                                               GC1, GC2, BR, BC, KS, KM, QW, &
                                               V(1), V(2), V(3), V(4), &
                                               tangent_ref)
        call fortfem_h1_geometry_jvp_fortad(J11, V(1), J21, V(2), J12, V(3), &
                                            J22, V(4), GR1, GR2, GC1, GC2, &
                                            BR, BC, KS, KM, QW, tangent_new)
        call same("jvp", tangent_ref, tangent_new, failures)

        call generated_bspline_h1_geometry_vjp(J11, J21, J12, J22, GR1, GR2, &
                                               GC1, GC2, BR, BC, KS, KM, QW, &
                                               u, g_ref(1), g_ref(2), &
                                               g_ref(3), g_ref(4))
        call fortfem_h1_geometry_vjp_fortad(J11, J21, J12, J22, GR1, GR2, &
                                            GC1, GC2, BR, BC, KS, KM, QW, &
                                            u, g_new(1), g_new(2), &
                                            g_new(3), g_new(4))
        call same("vjp jacobian_11", g_ref(1), g_new(1), failures)
        call same("vjp jacobian_21", g_ref(2), g_new(2), failures)
        call same("vjp jacobian_12", g_ref(3), g_new(3), failures)
        call same("vjp jacobian_22", g_ref(4), g_new(4), failures)

        call same("adjoint identity", dot_product(g_new, V), tangent_new*u, &
                  failures)
    end subroutine check_products

    subroutine check_differences(failures)
        !! The tangent against central differences of the primal.
        !!
        !! This is what makes the comparison mean something. Agreement with
        !! fortsym would hold even if the fortad primal stated a different
        !! computation, as long as it stated it consistently; differences of
        !! the primal in tools/fortad/kernels catch that, because they use no
        !! derivative machinery at all.
        integer, intent(inout) :: failures
        real(dp) :: h, plus, minus, tangent, fd

        call fortfem_h1_geometry_jvp_fortad(J11, V(1), J21, V(2), J12, V(3), &
                                            J22, V(4), GR1, GR2, GC1, GC2, &
                                            BR, BC, KS, KM, QW, tangent)
        h = 1.0e-6_dp
        call primal(J11 + h*V(1), J21 + h*V(2), J12 + h*V(3), J22 + h*V(4), plus)
        call primal(J11 - h*V(1), J21 - h*V(2), J12 - h*V(3), J22 - h*V(4), minus)
        fd = (plus - minus)/(2.0_dp*h)
        if (abs(tangent - fd) > 1.0e-6_dp*max(1.0_dp, abs(fd))) then
            print *, "FAIL central differences: ", fd, " vs ", tangent
            failures = failures + 1
        else
            print *, "pass central differences"
        end if
    end subroutine check_differences

    subroutine primal(j11, j21, j12, j22, contribution)
        !! The same computation the generated kernels differentiate, written
        !! out here so the difference check owes nothing to either generator.
        real(dp), intent(in) :: j11, j21, j12, j22
        real(dp), intent(out) :: contribution
        real(dp) :: determinant, row_1, row_2, column_1, column_2

        determinant = j11*j22 - j12*j21
        row_1 = (j22*GR1 - j21*GR2)/determinant
        row_2 = (-j12*GR1 + j11*GR2)/determinant
        column_1 = (j22*GC1 - j21*GC2)/determinant
        column_2 = (-j12*GC1 + j11*GC2)/determinant
        contribution = QW*determinant*(KS*(row_1*column_1 + row_2*column_2) &
                                       + KM*BR*BC)
    end subroutine primal

    subroutine same(label, reference, produced, failures)
        !! Agreement to rounding. The two kernels evaluate the same expression
        !! in a different order, so they are not bit-identical.
        character(*), intent(in) :: label
        real(dp), intent(in) :: reference, produced
        integer, intent(inout) :: failures

        if (abs(produced - reference) > 1.0e-12_dp*max(1.0_dp, abs(reference))) then
            print *, "FAIL ", label, ": ", reference, " vs ", produced
            failures = failures + 1
        else
            print *, "pass ", label
        end if
    end subroutine same

end program test_fortad_h1_geometry
