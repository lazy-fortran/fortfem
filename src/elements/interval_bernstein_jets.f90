module fortfem_interval_bernstein_jets
    !! Degree dispatch for the generated Bernstein jets on [0, 1]:
    !! value(i), first(i), second(i) of B_(i-1)^p(t), i = 1..p + 1.
    use fortfem_kinds, only: dp
    use fortfem_generated_interval_bernstein_jet_0, only: &
        generated_interval_bernstein_jet_0
    use fortfem_generated_interval_bernstein_jet_1, only: &
        generated_interval_bernstein_jet_1
    use fortfem_generated_interval_bernstein_jet_2, only: &
        generated_interval_bernstein_jet_2
    use fortfem_generated_interval_bernstein_jet_3, only: &
        generated_interval_bernstein_jet_3
    use fortfem_generated_interval_bernstein_jet_4, only: &
        generated_interval_bernstein_jet_4
    use fortfem_generated_interval_bernstein_jet_5, only: &
        generated_interval_bernstein_jet_5
    use fortfem_generated_interval_bernstein_jet_6, only: &
        generated_interval_bernstein_jet_6
    use fortfem_generated_interval_bernstein_jet_7, only: &
        generated_interval_bernstein_jet_7
    use fortfem_generated_interval_bernstein_jet_8, only: &
        generated_interval_bernstein_jet_8
    implicit none
    private

    integer, parameter, public :: INTERVAL_BERNSTEIN_MAX_DEGREE = 8

    public :: interval_bernstein_jet

contains

    pure subroutine interval_bernstein_jet(degree, t, value, first, second, status)
        integer, intent(in) :: degree
        real(dp), intent(in) :: t
        real(dp), intent(out) :: value(:), first(:), second(:)
        integer, intent(out) :: status

        status = 1
        if (degree < 0 .or. degree > INTERVAL_BERNSTEIN_MAX_DEGREE) return
        if (size(value) /= degree + 1 .or. size(first) /= degree + 1 .or. &
            size(second) /= degree + 1) return
        select case (degree)
        case (0)
            call generated_interval_bernstein_jet_0(t, value, first, second)
        case (1)
            call generated_interval_bernstein_jet_1(t, value, first, second)
        case (2)
            call generated_interval_bernstein_jet_2(t, value, first, second)
        case (3)
            call generated_interval_bernstein_jet_3(t, value, first, second)
        case (4)
            call generated_interval_bernstein_jet_4(t, value, first, second)
        case (5)
            call generated_interval_bernstein_jet_5(t, value, first, second)
        case (6)
            call generated_interval_bernstein_jet_6(t, value, first, second)
        case (7)
            call generated_interval_bernstein_jet_7(t, value, first, second)
        case (8)
            call generated_interval_bernstein_jet_8(t, value, first, second)
        end select
        status = 0
    end subroutine interval_bernstein_jet

end module fortfem_interval_bernstein_jets
