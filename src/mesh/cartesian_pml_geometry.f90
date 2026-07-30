module fortfem_cartesian_pml_geometry
    !! Elementwise polynomial stretch for Cartesian PML layers.
    use fortfem_kinds, only: dp
    use fortfem_generated_cartesian_pml_attenuation_jvp, only: &
        generated_cartesian_pml_attenuation_jvp
    use fortfem_generated_cartesian_pml_attenuation_vjp, only: &
        generated_cartesian_pml_attenuation_vjp
    implicit none
    private

    public :: build_cartesian_pml_element_stretch
    public :: build_cartesian_pml_element_stretch_jvp
    public :: build_cartesian_pml_element_stretch_vjp

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
        if (.not. valid_pml_geometry_inputs( &
            vertices, cells, physical_min, physical_max, outer_min, outer_max, &
            wave_number, sigma_max, polynomial_degree)) return

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

    subroutine build_cartesian_pml_element_stretch_jvp( &
            vertices, cells, physical_min, physical_max, outer_min, &
            outer_max, wave_number, sigma_max, polynomial_degree, &
            vertices_dot, physical_min_dot, physical_max_dot, outer_min_dot, &
            outer_max_dot, wave_number_dot, sigma_max_dot, stretch_dot, status)
        real(dp), intent(in) :: vertices(:, :), vertices_dot(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: physical_min(:), physical_max(:)
        real(dp), intent(in) :: outer_min(:), outer_max(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        real(dp), intent(in) :: physical_min_dot(:), physical_max_dot(:)
        real(dp), intent(in) :: outer_min_dot(:), outer_max_dot(:)
        real(dp), intent(in) :: wave_number_dot, sigma_max_dot
        complex(dp), allocatable, intent(out) :: stretch_dot(:, :)
        integer, intent(out) :: status

        real(dp) :: attenuation_dot, centroid, centroid_dot
        real(dp) :: distance, distance_dot, layer_width, layer_width_dot
        integer :: cell, coordinate, dimension

        status = 1
        if (allocated(stretch_dot)) deallocate(stretch_dot)
        dimension = size(vertices, 1)
        if (.not. valid_pml_geometry_inputs( &
            vertices, cells, physical_min, physical_max, outer_min, outer_max, &
            wave_number, sigma_max, polynomial_degree)) return
        if (any(shape(vertices_dot) /= shape(vertices))) return
        if (size(physical_min_dot) /= dimension) return
        if (size(physical_max_dot) /= dimension) return
        if (size(outer_min_dot) /= dimension) return
        if (size(outer_max_dot) /= dimension) return

        allocate (stretch_dot(dimension, size(cells, 2)))
        stretch_dot = cmplx(0.0_dp, 0.0_dp, dp)
        do cell = 1, size(cells, 2)
            do coordinate = 1, dimension
                centroid = sum(vertices(coordinate, cells(:, cell)))/ &
                    real(size(cells, 1), dp)
                centroid_dot = sum(vertices_dot(coordinate, cells(:, cell)))/ &
                    real(size(cells, 1), dp)
                distance = 0.0_dp
                distance_dot = 0.0_dp
                layer_width = 1.0_dp
                layer_width_dot = 0.0_dp
                if (centroid < physical_min(coordinate)) then
                    distance = physical_min(coordinate) - centroid
                    distance_dot = physical_min_dot(coordinate) - centroid_dot
                    layer_width = &
                        physical_min(coordinate) - outer_min(coordinate)
                    layer_width_dot = &
                        physical_min_dot(coordinate) - outer_min_dot(coordinate)
                else if (centroid > physical_max(coordinate)) then
                    distance = centroid - physical_max(coordinate)
                    distance_dot = centroid_dot - physical_max_dot(coordinate)
                    layer_width = &
                        outer_max(coordinate) - physical_max(coordinate)
                    layer_width_dot = &
                        outer_max_dot(coordinate) - physical_max_dot(coordinate)
                end if
                if (distance <= 0.0_dp) cycle
                call generated_cartesian_pml_attenuation_jvp( &
                    distance, layer_width, wave_number, sigma_max, &
                    real(polynomial_degree, dp), distance_dot, layer_width_dot, &
                    wave_number_dot, sigma_max_dot, attenuation_dot)
                stretch_dot(coordinate, cell) = &
                    cmplx(0.0_dp, attenuation_dot, dp)
            end do
        end do
        status = 0
    end subroutine build_cartesian_pml_element_stretch_jvp

    subroutine build_cartesian_pml_element_stretch_vjp( &
            vertices, cells, physical_min, physical_max, outer_min, &
            outer_max, wave_number, sigma_max, polynomial_degree, stretch_bar, &
            vertices_bar, physical_min_bar, physical_max_bar, outer_min_bar, &
            outer_max_bar, wave_number_bar, sigma_max_bar, status)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: physical_min(:), physical_max(:)
        real(dp), intent(in) :: outer_min(:), outer_max(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: stretch_bar(:, :)
        real(dp), intent(out) :: vertices_bar(:, :)
        real(dp), intent(out) :: physical_min_bar(:), physical_max_bar(:)
        real(dp), intent(out) :: outer_min_bar(:), outer_max_bar(:)
        real(dp), intent(out) :: wave_number_bar, sigma_max_bar
        integer, intent(out) :: status

        real(dp) :: attenuation_bar, centroid, centroid_bar
        real(dp) :: distance, distance_bar, layer_width, layer_width_bar
        real(dp) :: local_sigma_max_bar, local_wave_number_bar
        integer :: cell, coordinate, dimension, vertex

        status = 1
        vertices_bar = 0.0_dp
        physical_min_bar = 0.0_dp
        physical_max_bar = 0.0_dp
        outer_min_bar = 0.0_dp
        outer_max_bar = 0.0_dp
        wave_number_bar = 0.0_dp
        sigma_max_bar = 0.0_dp
        dimension = size(vertices, 1)
        if (.not. valid_pml_geometry_inputs( &
            vertices, cells, physical_min, physical_max, outer_min, outer_max, &
            wave_number, sigma_max, polynomial_degree)) return
        if (any(shape(stretch_bar) /= [dimension, size(cells, 2)])) return
        if (any(shape(vertices_bar) /= shape(vertices))) return
        if (size(physical_min_bar) /= dimension) return
        if (size(physical_max_bar) /= dimension) return
        if (size(outer_min_bar) /= dimension) return
        if (size(outer_max_bar) /= dimension) return

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
                attenuation_bar = aimag(stretch_bar(coordinate, cell))
                call generated_cartesian_pml_attenuation_vjp( &
                    distance, layer_width, wave_number, sigma_max, &
                    real(polynomial_degree, dp), attenuation_bar, distance_bar, &
                    layer_width_bar, local_wave_number_bar, local_sigma_max_bar)
                wave_number_bar = wave_number_bar + local_wave_number_bar
                sigma_max_bar = sigma_max_bar + local_sigma_max_bar
                if (centroid < physical_min(coordinate)) then
                    centroid_bar = -distance_bar
                    physical_min_bar(coordinate) = &
                        physical_min_bar(coordinate) + distance_bar + &
                        layer_width_bar
                    outer_min_bar(coordinate) = &
                        outer_min_bar(coordinate) - layer_width_bar
                else
                    centroid_bar = distance_bar
                    physical_max_bar(coordinate) = &
                        physical_max_bar(coordinate) - distance_bar - &
                        layer_width_bar
                    outer_max_bar(coordinate) = &
                        outer_max_bar(coordinate) + layer_width_bar
                end if
                do vertex = 1, size(cells, 1)
                    vertices_bar(coordinate, cells(vertex, cell)) = &
                        vertices_bar(coordinate, cells(vertex, cell)) + &
                        centroid_bar/real(size(cells, 1), dp)
                end do
            end do
        end do
        status = 0
    end subroutine build_cartesian_pml_element_stretch_vjp

    pure logical function valid_pml_geometry_inputs( &
            vertices, cells, physical_min, physical_max, outer_min, outer_max, &
            wave_number, sigma_max, polynomial_degree) result(valid)
        real(dp), intent(in) :: vertices(:, :)
        integer, intent(in) :: cells(:, :)
        real(dp), intent(in) :: physical_min(:), physical_max(:)
        real(dp), intent(in) :: outer_min(:), outer_max(:)
        real(dp), intent(in) :: wave_number, sigma_max
        integer, intent(in) :: polynomial_degree

        integer :: coordinate, dimension

        valid = .false.
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
        valid = .true.
    end function valid_pml_geometry_inputs

end module fortfem_cartesian_pml_geometry
