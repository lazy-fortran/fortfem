module fortfem_scalar_helmholtz_pml_slab_1d
    !! P1 finite elements for a scalar Helmholtz slab terminated by a
    !! polynomial complex-coordinate perfectly matched layer.
    use fortfem_cartesian_helmholtz_pml, only: &
        cartesian_scalar_helmholtz_pml_coefficients
    use fortfem_kinds, only: dp
    use fortnum_linalg, only: dense_solve
    implicit none
    private

    public :: solve_scalar_helmholtz_pml_slab_1d

contains

    subroutine solve_scalar_helmholtz_pml_slab_1d( &
            nodes, physical_end, wave_number, sigma_max, polynomial_degree, &
            left_value, solution, status)
        real(dp), intent(in) :: nodes(:)
        real(dp), intent(in) :: physical_end, wave_number, sigma_max
        integer, intent(in) :: polynomial_degree
        complex(dp), intent(in) :: left_value
        complex(dp), allocatable, intent(out) :: solution(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:)
        complex(dp) :: gradient(3), mass, stretch(3)
        complex(dp) :: element_matrix(2, 2)
        real(dp) :: element_length, midpoint, pml_width, sigma
        integer :: element, info, local_status

        status = 1
        if (allocated(solution)) deallocate(solution)
        if (size(nodes) < 3) return
        if (any(nodes(2:) <= nodes(:size(nodes) - 1))) return
        if (physical_end <= nodes(1) .or. physical_end >= nodes(size(nodes))) &
            return
        if (wave_number <= 0.0_dp .or. sigma_max <= 0.0_dp) return
        if (polynomial_degree < 1) return

        allocate(matrix(size(nodes), size(nodes)))
        allocate(right_hand_side(size(nodes)), solution(size(nodes)))
        matrix = cmplx(0.0_dp, 0.0_dp, dp)
        right_hand_side = cmplx(0.0_dp, 0.0_dp, dp)
        pml_width = nodes(size(nodes)) - physical_end
        do element = 1, size(nodes) - 1
            element_length = nodes(element + 1) - nodes(element)
            midpoint = 0.5_dp*(nodes(element) + nodes(element + 1))
            sigma = 0.0_dp
            if (midpoint > physical_end) then
                sigma = sigma_max*((midpoint - physical_end)/pml_width)** &
                    polynomial_degree
            end if
            stretch = cmplx(1.0_dp, 0.0_dp, dp)
            stretch(1) = cmplx(1.0_dp, sigma/wave_number, dp)
            call cartesian_scalar_helmholtz_pml_coefficients( &
                stretch, gradient, mass, local_status)
            if (local_status /= 0) return
            element_matrix = gradient(1)/element_length*reshape( &
                [1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp], [2, 2]) - &
                wave_number**2*mass*element_length/6.0_dp*reshape( &
                [2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
            matrix(element:element + 1, element:element + 1) = &
                matrix(element:element + 1, element:element + 1) + &
                element_matrix
        end do

        right_hand_side = right_hand_side - matrix(:, 1)*left_value
        matrix(:, 1) = cmplx(0.0_dp, 0.0_dp, dp)
        matrix(1, :) = cmplx(0.0_dp, 0.0_dp, dp)
        matrix(1, 1) = cmplx(1.0_dp, 0.0_dp, dp)
        right_hand_side(1) = left_value
        matrix(:, size(nodes)) = cmplx(0.0_dp, 0.0_dp, dp)
        matrix(size(nodes), :) = cmplx(0.0_dp, 0.0_dp, dp)
        matrix(size(nodes), size(nodes)) = cmplx(1.0_dp, 0.0_dp, dp)
        right_hand_side(size(nodes)) = cmplx(0.0_dp, 0.0_dp, dp)
        call dense_solve(matrix, right_hand_side, solution, info)
        if (info /= 0) then
            status = 2
            return
        end if
        status = 0
    end subroutine solve_scalar_helmholtz_pml_slab_1d

end module fortfem_scalar_helmholtz_pml_slab_1d
