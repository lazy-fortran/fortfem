module fortfem_bspline_polar
    !! Polar-axis extraction for periodic tensor-product H1 splines.
    !!
    !! This is the Type-1 scalar extraction of Toshniwal and Hughes,
    !! CMAME 376 (2021) 113576, equations (64)--(70). The companion
    !! differential-form extractions are deliberately separate: a scalar
    !! edge collapse alone is not a polar FEEC complex.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: build_bspline_polar_h1_extraction
    public :: build_bspline_polar_feec_2d_operators
    public :: build_bspline_polar_feec_2d_extractions
    public :: evaluate_periodic_bspline_basis

contains

    subroutine evaluate_periodic_bspline_basis( &
            basis_count, degree, coordinate, values, derivatives, status)
        !! Uniform periodic B-spline basis on [0,1).
        integer, intent(in) :: basis_count, degree
        real(dp), intent(in) :: coordinate
        real(dp), allocatable, intent(out) :: values(:), derivatives(:)
        integer, intent(out) :: status

        real(dp) :: scaled_coordinate, value
        integer :: extended_basis, first_basis, periodic_basis

        status = 1
        if (degree < 0) return
        if (basis_count < degree + 1) return
        if (coordinate < 0.0_dp .or. coordinate > 1.0_dp) return
        scaled_coordinate = real(basis_count, dp)*coordinate
        if (coordinate == 1.0_dp) scaled_coordinate = 0.0_dp
        allocate(values(basis_count), derivatives(basis_count))
        values = 0.0_dp
        derivatives = 0.0_dp
        first_basis = floor(scaled_coordinate) - degree
        do extended_basis = first_basis, first_basis + degree
            periodic_basis = modulo(extended_basis, basis_count) + 1
            value = cardinal_bspline( &
                degree, scaled_coordinate - real(extended_basis, dp))
            values(periodic_basis) = values(periodic_basis) + value
            if (degree == 0) cycle
            derivatives(periodic_basis) = derivatives(periodic_basis) + &
                real(basis_count, dp)*( &
                cardinal_bspline( &
                degree - 1, scaled_coordinate - &
                real(extended_basis, dp)) - &
                cardinal_bspline( &
                degree - 1, scaled_coordinate - &
                real(extended_basis + 1, dp)))
        end do
        status = 0
    end subroutine evaluate_periodic_bspline_basis

    recursive pure real(dp) function cardinal_bspline(degree, coordinate) &
            result(value)
        integer, intent(in) :: degree
        real(dp), intent(in) :: coordinate

        if (degree == 0) then
            if (coordinate >= 0.0_dp .and. coordinate < 1.0_dp) then
                value = 1.0_dp
            else
                value = 0.0_dp
            end if
            return
        end if
        value = coordinate/real(degree, dp)* &
            cardinal_bspline(degree - 1, coordinate) + &
            (real(degree + 1, dp) - coordinate)/real(degree, dp)* &
            cardinal_bspline(degree - 1, coordinate - 1.0_dp)
    end function cardinal_bspline

    subroutine build_bspline_polar_feec_2d_operators( &
            azimuth_count, radial_count, gradient, curl, status)
        !! Type-1 polar 0--1--2 form incidence complex.
        !!
        !! The two pole-edge rows encode differences between the three polar
        !! H1 coefficients. Spoke rows use the same barycentric fan as the H1
        !! extraction. This is the locally exact incidence construction
        !! underlying Toshniwal--Hughes equations (74)--(101).
        integer, intent(in) :: azimuth_count, radial_count
        real(dp), allocatable, intent(out) :: gradient(:, :), curl(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: extraction(:, :)
        real(dp) :: delta(3)
        integer :: angular_edge, azimuth, face, next_azimuth, radial
        integer :: current_vertex, next_vertex, previous_vertex
        integer :: current_radial, next_radial, previous_angular
        integer :: radial_edge_count, special_edge

        call build_bspline_polar_h1_extraction( &
            azimuth_count, radial_count, extraction, status)
        if (status /= 0) return
        allocate(gradient( &
            2 + 2*azimuth_count*(radial_count - 2), &
            3 + azimuth_count*(radial_count - 2)))
        allocate(curl( &
            azimuth_count*(radial_count - 2), size(gradient, 1)))
        gradient = 0.0_dp
        curl = 0.0_dp

        radial_edge_count = azimuth_count*(radial_count - 2)
        do radial = 1, radial_count - 2
            do azimuth = 1, azimuth_count
                current_vertex = polar_vertex(radial, azimuth)
                current_radial = &
                    (radial - 1)*azimuth_count + azimuth
                if (radial == 1) then
                    gradient(current_radial, 1:3) = &
                        -extraction(1:3, azimuth_count + azimuth)
                else
                    previous_vertex = polar_vertex(radial - 1, azimuth)
                    gradient(current_radial, previous_vertex) = -1.0_dp
                end if
                gradient(current_radial, current_vertex) = 1.0_dp

                next_azimuth = modulo(azimuth, azimuth_count) + 1
                next_vertex = polar_vertex(radial, next_azimuth)
                angular_edge = radial_edge_count + 2 + &
                    (radial - 1)*azimuth_count + azimuth
                gradient(angular_edge, current_vertex) = -1.0_dp
                gradient(angular_edge, next_vertex) = 1.0_dp

                face = (radial - 1)*azimuth_count + azimuth
                next_radial = &
                    (radial - 1)*azimuth_count + next_azimuth
                curl(face, current_radial) = 1.0_dp
                curl(face, next_radial) = -1.0_dp
                curl(face, angular_edge) = 1.0_dp
                if (radial == 1) then
                    delta = &
                        extraction(1:3, azimuth_count + next_azimuth) - &
                        extraction(1:3, azimuth_count + azimuth)
                    curl(face, radial_edge_count + 1) = -delta(2)
                    curl(face, radial_edge_count + 2) = -delta(3)
                else
                    previous_angular = angular_edge - azimuth_count
                    curl(face, previous_angular) = -1.0_dp
                end if
            end do
        end do
        special_edge = radial_edge_count + 1
        gradient(special_edge, 1:2) = [-1.0_dp, 1.0_dp]
        gradient(special_edge + 1, [1, 3]) = [-1.0_dp, 1.0_dp]
        status = 0

    contains

        pure integer function polar_vertex(radial_ring, azimuth_index) &
                result(index)
            integer, intent(in) :: radial_ring, azimuth_index

            index = 3 + (radial_ring - 1)*azimuth_count + azimuth_index
        end function polar_vertex

    end subroutine build_bspline_polar_feec_2d_operators

    subroutine build_bspline_polar_feec_2d_extractions( &
            azimuth_count, radial_count, h1_extraction, hcurl_extraction, &
            l2_extraction, status)
        !! Sparse-block basis extractions from the periodic tensor space.
        !!
        !! Rows are polar basis functions and columns are tensor-product basis
        !! functions. This is the Type-1 C1 construction of
        !! Toshniwal--Hughes (2021), equations (68), (90), and (102).
        integer, intent(in) :: azimuth_count, radial_count
        real(dp), allocatable, intent(out) :: h1_extraction(:, :)
        real(dp), allocatable, intent(out) :: hcurl_extraction(:, :)
        real(dp), allocatable, intent(out) :: l2_extraction(:, :)
        integer, intent(out) :: status

        real(dp) :: angular_delta
        integer :: angular_column_offset, angular_row_offset, azimuth
        integer :: polar_row, radial, radial_edge_count, special_row
        integer :: tensor_column, tensor_radial_edge_count

        call build_bspline_polar_h1_extraction( &
            azimuth_count, radial_count, h1_extraction, status)
        if (status /= 0) return
        radial_edge_count = azimuth_count*(radial_count - 2)
        tensor_radial_edge_count = azimuth_count*(radial_count - 1)
        allocate(hcurl_extraction( &
            2 + 2*radial_edge_count, &
            azimuth_count*(2*radial_count - 1)))
        allocate(l2_extraction( &
            radial_edge_count, azimuth_count*(radial_count - 1)))
        hcurl_extraction = 0.0_dp
        l2_extraction = 0.0_dp

        do radial = 1, radial_count - 2
            do azimuth = 1, azimuth_count
                polar_row = (radial - 1)*azimuth_count + azimuth
                tensor_column = radial*azimuth_count + azimuth
                hcurl_extraction(polar_row, tensor_column) = 1.0_dp

                angular_row_offset = radial_edge_count + 2
                polar_row = angular_row_offset + &
                    (radial - 1)*azimuth_count + azimuth
                angular_column_offset = tensor_radial_edge_count
                tensor_column = angular_column_offset + &
                    (radial + 1)*azimuth_count + azimuth
                hcurl_extraction(polar_row, tensor_column) = 1.0_dp

                polar_row = (radial - 1)*azimuth_count + azimuth
                tensor_column = radial*azimuth_count + azimuth
                l2_extraction(polar_row, tensor_column) = 1.0_dp
            end do
        end do

        do polar_row = 1, 2
            special_row = radial_edge_count + polar_row
            do azimuth = 1, azimuth_count
                hcurl_extraction(special_row, azimuth) = &
                    h1_extraction(polar_row + 1, &
                    azimuth_count + azimuth) - 1.0_dp/3.0_dp
                angular_delta = h1_extraction( &
                    polar_row + 1, azimuth_count + &
                    modulo(azimuth, azimuth_count) + 1) - &
                    h1_extraction( &
                    polar_row + 1, azimuth_count + azimuth)
                tensor_column = tensor_radial_edge_count + &
                    azimuth_count + azimuth
                hcurl_extraction( &
                    special_row, tensor_column) = angular_delta
            end do
        end do
        status = 0
    end subroutine build_bspline_polar_feec_2d_extractions

    subroutine build_bspline_polar_h1_extraction( &
            azimuth_count, radial_count, extraction, status)
        integer, intent(in) :: azimuth_count, radial_count
        real(dp), allocatable, intent(out) :: extraction(:, :)
        integer, intent(out) :: status

        real(dp), parameter :: one_third = 1.0_dp/3.0_dp
        real(dp) :: angle, barycentric(3), pi
        integer :: azimuth, polar_dof, radial, tensor_dof

        status = 1
        if (azimuth_count < 3 .or. radial_count < 5) return
        allocate(extraction( &
            azimuth_count*(radial_count - 2) + 3, &
            azimuth_count*radial_count))
        extraction = 0.0_dp
        pi = acos(-1.0_dp)
        do azimuth = 1, azimuth_count
            tensor_dof = azimuth
            extraction(1:3, tensor_dof) = one_third

            angle = 2.0_dp*pi - &
                real(2*azimuth - 1, dp)*pi/real(azimuth_count, dp)
            barycentric = [ &
                one_third + cos(angle)/3.0_dp, &
                one_third - cos(angle)/6.0_dp + &
                sqrt(3.0_dp)*sin(angle)/6.0_dp, &
                one_third - cos(angle)/6.0_dp - &
                sqrt(3.0_dp)*sin(angle)/6.0_dp]
            tensor_dof = azimuth_count + azimuth
            extraction(1:3, tensor_dof) = barycentric
        end do
        polar_dof = 3
        do radial = 3, radial_count
            do azimuth = 1, azimuth_count
                polar_dof = polar_dof + 1
                tensor_dof = azimuth + (radial - 1)*azimuth_count
                extraction(polar_dof, tensor_dof) = 1.0_dp
            end do
        end do
        status = 0
    end subroutine build_bspline_polar_h1_extraction

end module fortfem_bspline_polar
