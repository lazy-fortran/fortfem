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

contains

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
        integer :: azimuth, face, next_azimuth, radial, row
        integer :: current_vertex, next_vertex
        integer :: bottom_edge, current_spoke, next_spoke, top_edge

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

        gradient(1, 1:2) = [-1.0_dp, 1.0_dp]
        gradient(2, [1, 3]) = [-1.0_dp, 1.0_dp]
        do azimuth = 1, azimuth_count
            next_azimuth = modulo(azimuth, azimuth_count) + 1
            current_vertex = polar_vertex(1, azimuth)
            next_vertex = polar_vertex(1, next_azimuth)
            current_spoke = 2 + azimuth
            bottom_edge = 2 + azimuth_count + azimuth

            gradient(current_spoke, 1:3) = &
                -extraction(1:3, azimuth_count + azimuth)
            gradient(current_spoke, current_vertex) = 1.0_dp
            gradient(bottom_edge, current_vertex) = -1.0_dp
            gradient(bottom_edge, next_vertex) = 1.0_dp

            delta = extraction(1:3, azimuth_count + next_azimuth) - &
                extraction(1:3, azimuth_count + azimuth)
            next_spoke = 2 + next_azimuth
            face = azimuth
            curl(face, current_spoke) = -1.0_dp
            curl(face, next_spoke) = 1.0_dp
            curl(face, bottom_edge) = -1.0_dp
            curl(face, 1) = delta(2)
            curl(face, 2) = delta(3)
        end do

        row = 2 + 2*azimuth_count
        face = azimuth_count
        do radial = 1, radial_count - 3
            do azimuth = 1, azimuth_count
                next_azimuth = modulo(azimuth, azimuth_count) + 1
                current_vertex = polar_vertex(radial, azimuth)
                next_vertex = polar_vertex(radial + 1, azimuth)
                row = row + 1
                gradient(row, current_vertex) = -1.0_dp
                gradient(row, next_vertex) = 1.0_dp
            end do
            do azimuth = 1, azimuth_count
                next_azimuth = modulo(azimuth, azimuth_count) + 1
                current_vertex = polar_vertex(radial + 1, azimuth)
                next_vertex = polar_vertex(radial + 1, next_azimuth)
                row = row + 1
                gradient(row, current_vertex) = -1.0_dp
                gradient(row, next_vertex) = 1.0_dp

                face = face + 1
                bottom_edge = &
                    2 + 2*azimuth_count*radial - azimuth_count + azimuth
                current_spoke = &
                    2 + 2*azimuth_count*radial + azimuth
                next_spoke = 2 + 2*azimuth_count*radial + next_azimuth
                top_edge = row
                curl(face, bottom_edge) = 1.0_dp
                curl(face, next_spoke) = 1.0_dp
                curl(face, top_edge) = -1.0_dp
                curl(face, current_spoke) = -1.0_dp
            end do
        end do
        status = 0

    contains

        pure integer function polar_vertex(radial_ring, azimuth_index) &
                result(index)
            integer, intent(in) :: radial_ring, azimuth_index

            index = 3 + (radial_ring - 1)*azimuth_count + azimuth_index
        end function polar_vertex

    end subroutine build_bspline_polar_feec_2d_operators

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
