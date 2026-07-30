module fortfem_cartesian_pml_geometry
    !! Elementwise polynomial stretch for Cartesian PML layers.
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: build_cartesian_pml_element_stretch

contains

    subroutine build_cartesian_pml_element_stretch( &
            vertices, cells, physical_min, physical_max, outer_min, &
            outer_max, wave_number, sigma_max, polynomial_degree, stretch, &
            status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: physical_min(:), physical_max(:)
        real(dp), intent(in) :: outer_min(:), outer_max(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), allocatable, intent(out) :: stretch(:, :)
        integer, intent(out) :: status

        real(dp) :: centroid, distance, layer_width, sigma
        integer :: cell, coordinate, dimension

        status = 1
        if (allocated(stretch)) deallocate(stretch)
        dimension = size(vertices, 1)
        if (dimension < 1 .or. dimension > 3) return
        if (size(vertices, 2) < 1) return
        if (size(cells, 1) < 2 .or. size(cells, 2) < 1) return
        if (size(physical_min) /= dimension) return
        if (size(physical_max) /= dimension) return
        if (size(outer_min) /= dimension) return
        if (size(outer_max) /= dimension) return
        if (any(cells < 1) .or. any(cells > size(vertices, 2))) return
        if (any(physical_min >= physical_max)) return
        if (any(outer_min >= physical_min)) return
        if (any(outer_max <= physical_max)) return
        if (wave_number <= 0.0_dp .or. sigma_max <= 0.0_dp) return
        if (polynomial_degree < 1) return
        do coordinate = 1, dimension
            if (any(vertices(coordinate, :) < outer_min(coordinate))) return
            if (any(vertices(coordinate, :) > outer_max(coordinate))) return
        end do

        allocate(stretch(dimension, size(cells, 2)))
        stretch = cmplx(1.0_dp, 0.0_dp, dp)
        do cell = 1, size(cells, 2)
            do coordinate = 1, dimension
                centroid = sum(vertices(coordinate, cells(:, cell)))/ &
                    real(size(cells, 1), dp)
                distance = 0.0_dp
                layer_width = 1.0_dp
                if (centroid < physical_min(coordinate)) then
                    distance = physical_min(coordinate) - centroid
                    layer_width = &
                        physical_min(coordinate) - outer_min(coordinate)
                else if (centroid > physical_max(coordinate)) then
                    distance = centroid - physical_max(coordinate)
                    layer_width = &
                        outer_max(coordinate) - physical_max(coordinate)
                end if
                if (distance <= 0.0_dp) cycle
                sigma = sigma_max*(distance/layer_width)**polynomial_degree
                stretch(coordinate, cell) = &
                    cmplx(1.0_dp, sigma/wave_number, dp)
            end do
        end do
        status = 0
    end subroutine build_cartesian_pml_element_stretch

end module fortfem_cartesian_pml_geometry
