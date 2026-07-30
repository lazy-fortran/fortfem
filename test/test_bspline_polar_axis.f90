program test_bspline_polar_axis
    use check, only: check_condition, check_summary
    use fortfem_api, only: build_bspline_polar_h1_extraction
    use fortfem_kinds, only: dp
    implicit none

    integer, parameter :: azimuth_count = 8
    integer, parameter :: radial_count = 5
    real(dp), allocatable :: extraction(:, :), polar_points(:, :)
    real(dp), allocatable :: tensor_points(:, :)
    real(dp) :: angle, radius(azimuth_count)
    integer :: azimuth, status

    call build_bspline_polar_h1_extraction( &
        azimuth_count, radial_count, extraction, status)
    call check_condition(status == 0, &
        "Type-1 polar H1 extraction accepts a periodic radial patch")
    call check_condition(all(shape(extraction) == [ &
        azimuth_count*(radial_count - 2) + 3, &
        azimuth_count*radial_count]), &
        "Polar H1 dimension matches Toshniwal-Hughes equation 74")
    call check_condition(minval(extraction) >= -epsilon(1.0_dp) .and. &
        maxval(extraction) <= 1.0_dp + epsilon(1.0_dp), &
        "Polar extraction is nonnegative")
    call check_condition(maxval(abs(sum(extraction, dim=1) - 1.0_dp)) < &
        2.0e-15_dp, "Polar extraction is column stochastic")
    call check_condition(gram_is_positive_definite(extraction), &
        "Polar extraction has full row rank")

    allocate(polar_points(2, size(extraction, 1)))
    polar_points = 0.0_dp
    polar_points(:, 1) = [0.0_dp, 2.0_dp/sqrt(3.0_dp)]
    polar_points(:, 2) = [-1.0_dp, -1.0_dp/sqrt(3.0_dp)]
    polar_points(:, 3) = [1.0_dp, -1.0_dp/sqrt(3.0_dp)]
    tensor_points = matmul(polar_points, extraction)
    call check_condition(maxval(abs( &
        tensor_points(:, :azimuth_count) - &
        spread(tensor_points(:, 1), 2, azimuth_count))) < 2.0e-15_dp, &
        "Polar extraction collapses the complete periodic edge to one point")
    do azimuth = 1, azimuth_count
        radius(azimuth) = norm2( &
            tensor_points(:, azimuth_count + azimuth) - &
            tensor_points(:, azimuth))
    end do
    call check_condition(minval(radius) > 0.0_dp .and. &
        maxval(abs(radius - radius(1))) < 2.0e-15_dp, &
        "First polar control ring is a nondegenerate uniform fan")
    angle = signed_area( &
        tensor_points(:, azimuth_count + 1), &
        tensor_points(:, azimuth_count + 2), tensor_points(:, 1))
    call check_condition(abs(angle) > 1.0e-2_dp, &
        "Polar control fan has a well-defined orientation")

    call check_summary("Isogeometric magnetic-axis polar extraction")

contains

    logical function gram_is_positive_definite(matrix) result(positive)
        real(dp), intent(in) :: matrix(:, :)
        real(dp) :: factor(size(matrix, 1), size(matrix, 1))
        real(dp) :: gram(size(matrix, 1), size(matrix, 1)), pivot
        integer :: column, row

        gram = matmul(matrix, transpose(matrix))
        factor = 0.0_dp
        positive = .true.
        do column = 1, size(gram, 1)
            pivot = gram(column, column) - &
                dot_product(factor(column, :column - 1), &
                factor(column, :column - 1))
            if (pivot <= 1.0e3_dp*epsilon(1.0_dp)) then
                positive = .false.
                return
            end if
            factor(column, column) = sqrt(pivot)
            do row = column + 1, size(gram, 1)
                factor(row, column) = ( &
                    gram(row, column) - &
                    dot_product(factor(row, :column - 1), &
                    factor(column, :column - 1)))/factor(column, column)
            end do
        end do
    end function gram_is_positive_definite

    pure real(dp) function signed_area(first, second, origin) result(area)
        real(dp), intent(in) :: first(2), second(2), origin(2)

        area = (first(1) - origin(1))*(second(2) - origin(2)) - &
            (first(2) - origin(2))*(second(1) - origin(1))
    end function signed_area

end program test_bspline_polar_axis
