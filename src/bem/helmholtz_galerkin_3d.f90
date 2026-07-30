module fortfem_helmholtz_galerkin_3d
    !! P0 Galerkin single-layer BEM for the outgoing three-dimensional
    !! Helmholtz kernel. The singularity subtraction
    !! exp(i*k*r)/r = 1/r + (exp(i*k*r)-1)/r leaves the singular Laplace
    !! contribution to the analytical diagonal treatment.
    use fortfem_kinds, only: dp
    use fortfem_laplace_galerkin_3d, only: &
        assemble_laplace_single_layer_p0_3d
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    implicit none

    private

    public :: assemble_helmholtz_single_layer_p0_3d
    public :: solve_helmholtz_dirichlet_p0_3d

    interface
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            complex(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine zgesv
    end interface

contains

    subroutine solve_helmholtz_dirichlet_p0_3d( &
            vertices, triangles, boundary_value, wave_number, &
            quadrature_degree, density, status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), intent(in) :: boundary_value
        complex(dp), allocatable, intent(out) :: density(:)
        integer, intent(out) :: status

        complex(dp), allocatable :: matrix(:, :), right_hand_side(:, :)
        integer, allocatable :: pivots(:)
        integer :: element, info

        status = 1
        call assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, &
            status)
        if (status /= 0) return
        allocate( &
            density(size(triangles, 2)), &
            right_hand_side(size(triangles, 2), 1), &
            pivots(size(triangles, 2)))
        do element = 1, size(triangles, 2)
            right_hand_side(element, 1) = boundary_value*triangle_area( &
                vertices(:, triangles(:, element)))
        end do
        call zgesv( &
            size(triangles, 2), 1, matrix, size(triangles, 2), pivots, &
            right_hand_side, size(triangles, 2), info)
        if (info /= 0) then
            status = 2
            return
        end if
        density = right_hand_side(:, 1)
        status = 0
    end subroutine solve_helmholtz_dirichlet_p0_3d

    subroutine assemble_helmholtz_single_layer_p0_3d( &
            vertices, triangles, wave_number, quadrature_degree, matrix, &
            status)
        real(dp), intent(in) :: vertices(:, :), wave_number
        integer, intent(in) :: triangles(:, :), quadrature_degree
        complex(dp), allocatable, intent(out) :: matrix(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: laplace_matrix(:, :)
        real(dp), allocatable :: eta(:), weights(:), xi(:)
        real(dp) :: first_jacobian, first_point_value(3), radius
        real(dp) :: second_jacobian
        real(dp) :: first_vertices(3, 3), second_vertices(3, 3)
        integer :: first, first_point, quadrature_status, second
        integer :: second_point

        status = 1
        if (wave_number < 0.0_dp) return
        call assemble_laplace_single_layer_p0_3d( &
            vertices, triangles, quadrature_degree, laplace_matrix, status)
        if (status /= 0) return
        call triangle_duffy_quadrature( &
            quadrature_degree, xi, eta, weights, quadrature_status)
        if (quadrature_status /= 0) return
        allocate(matrix(size(laplace_matrix, 1), size(laplace_matrix, 2)))
        matrix = cmplx(laplace_matrix, 0.0_dp, dp)
        do first = 1, size(triangles, 2)
            first_vertices = vertices(:, triangles(:, first))
            first_jacobian = 2.0_dp*triangle_area(first_vertices)
            do second = 1, first
                second_vertices = vertices(:, triangles(:, second))
                second_jacobian = 2.0_dp*triangle_area(second_vertices)
                do first_point = 1, size(weights)
                    first_point_value = triangle_point( &
                        first_vertices, xi(first_point), eta(first_point))
                    do second_point = 1, size(weights)
                        radius = norm2(first_point_value - triangle_point( &
                            second_vertices, xi(second_point), &
                            eta(second_point)))
                        matrix(first, second) = matrix(first, second) + &
                            first_jacobian*second_jacobian* &
                            weights(first_point)*weights(second_point)* &
                            smooth_remainder(wave_number, radius)
                    end do
                end do
                matrix(second, first) = matrix(first, second)
            end do
        end do
        status = 0
    end subroutine assemble_helmholtz_single_layer_p0_3d

    pure function smooth_remainder(wave_number, radius) result(value)
        real(dp), intent(in) :: wave_number, radius
        complex(dp) :: value

        if (radius <= 64.0_dp*epsilon(1.0_dp)) then
            value = cmplx(0.0_dp, wave_number, dp)/(4.0_dp*acos(-1.0_dp))
        else
            value = (exp(cmplx(0.0_dp, wave_number*radius, dp)) - 1.0_dp)/ &
                (4.0_dp*acos(-1.0_dp)*radius)
        end if
    end function smooth_remainder

    pure function triangle_point(vertices, xi, eta) result(point)
        real(dp), intent(in) :: vertices(3, 3), xi, eta
        real(dp) :: point(3)

        point = vertices(:, 1) + xi*(vertices(:, 2) - vertices(:, 1)) + &
            eta*(vertices(:, 3) - vertices(:, 1))
    end function triangle_point

    pure function triangle_area(vertices) result(area)
        real(dp), intent(in) :: vertices(3, 3)
        real(dp) :: area

        area = 0.5_dp*norm2(cross_product( &
            vertices(:, 2) - vertices(:, 1), &
            vertices(:, 3) - vertices(:, 1)))
    end function triangle_area

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

end module fortfem_helmholtz_galerkin_3d
