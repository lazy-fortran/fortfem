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

contains

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
