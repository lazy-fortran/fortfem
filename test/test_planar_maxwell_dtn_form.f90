program test_planar_maxwell_dtn_form
    use check, only: check_condition, check_summary
    use fortfem_api, only: apply_planar_maxwell_dtn, &
        assemble_planar_maxwell_dtn_form
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: nx = 5, ny = 5
    real(dp), parameter :: length_x = 3.0_dp, length_y = 2.5_dp
    real(dp), parameter :: wave_number = 4.7_dp
    complex(dp), allocatable :: form(:, :), samples(:)
    complex(dp) :: derivative(2, nx, ny), trace(2, nx, ny)
    complex(dp) :: weak_expected(2*nx*ny)
    real(dp) :: weight
    integer :: component, i, index, j, status
    logical :: all_passed

    all_passed = .true.
    index = 0
    do j = 1, ny
        do i = 1, nx
            do component = 1, 2
                index = index + 1
                trace(component, i, j) = cmplx( &
                    sin(0.31_dp*real(index, dp)), &
                    cos(0.23_dp*real(index, dp)), dp)
            end do
        end do
    end do
    call apply_planar_maxwell_dtn( &
        trace, wave_number, length_x, length_y, derivative, status)
    if (status /= 0) error stop "strong Maxwell DtN oracle failed"
    call assemble_planar_maxwell_dtn_form( &
        nx, ny, wave_number, length_x, length_y, form, status)
    if (status /= 0) error stop "Maxwell DtN form assembly failed"
    samples = reshape(trace, [size(trace)])
    weight = length_x*length_y/real(nx*ny, dp)
    weak_expected = weight*reshape(derivative, [size(derivative)])
    call record_condition( &
        maxval(abs(matmul(form, samples) - weak_expected)) < 3.0e-12_dp, &
        "Maxwell DtN weak form matches weighted strong application")
    call record_condition(maxval(abs(form - transpose(form))) < 3.0e-12_dp, &
        "reciprocal planar Maxwell DtN form is complex symmetric")
    call check_summary("Planar Maxwell DtN weak form")
    if (.not. all_passed) error stop 1

contains

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        all_passed = all_passed .and. condition
        call check_condition(condition, description)
    end subroutine record_condition

end program test_planar_maxwell_dtn_form
