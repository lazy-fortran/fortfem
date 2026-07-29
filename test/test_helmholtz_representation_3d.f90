program test_helmholtz_representation_3d
    use check, only: check_condition, check_summary
    use fortfem_api, only: evaluate_helmholtz_representation_triangles_3d, &
        generate_torus_surface_mesh
    use fortfem_kinds, only: dp
    implicit none

    real(dp), parameter :: major_radius = 2.0_dp
    real(dp), parameter :: minor_radius = 0.6_dp
    real(dp), parameter :: wave_number = 0.2_dp
    real(dp), allocatable :: vertices(:, :)
    integer, allocatable :: triangles(:, :)
    complex(dp), allocatable :: trace(:), flux(:)
    complex(dp) :: exact, numerical
    real(dp) :: centroid(3), displacement(3), first_edge(3), normal(3)
    real(dp) :: second_edge(3), source(3), target(3), radius
    integer :: element, point, status
    logical :: all_passed

    all_passed = .true.
    source = [major_radius, 0.0_dp, 0.0_dp]
    target = [0.0_dp, 0.0_dp, 0.0_dp]
    call generate_torus_surface_mesh( &
        major_radius, minor_radius, 96, 72, vertices, triangles)
    allocate(trace(size(vertices, 2)), flux(size(triangles, 2)))
    do point = 1, size(vertices, 2)
        displacement = vertices(:, point) - source
        radius = norm2(displacement)
        trace(point) = exp(cmplx(0.0_dp, wave_number*radius, dp))/radius
    end do
    do element = 1, size(triangles, 2)
        centroid = sum(vertices(:, triangles(:, element)), dim=2)/3.0_dp
        first_edge = vertices(:, triangles(2, element)) - &
            vertices(:, triangles(1, element))
        second_edge = vertices(:, triangles(3, element)) - &
            vertices(:, triangles(1, element))
        normal = cross_product(first_edge, second_edge)
        normal = normal/norm2(normal)
        displacement = centroid - source
        radius = norm2(displacement)
        flux(element) = exp(cmplx(0.0_dp, wave_number*radius, dp))* &
            cmplx(-1.0_dp/radius**2, wave_number/radius, dp)* &
            dot_product(displacement/radius, normal)
    end do

    call evaluate_helmholtz_representation_triangles_3d( &
        vertices, triangles, trace, flux, target, wave_number, 8, &
        numerical, status)
    exact = exp(cmplx(0.0_dp, wave_number*major_radius, dp))/major_radius
    call record_condition(status == 0, &
        "Three-dimensional Helmholtz representation evaluates")
    call record_condition(abs(numerical - exact) < 2.0e-2_dp*abs(exact), &
        "Toroidal Helmholtz BEM reproduces an outgoing point source")

    call check_summary("Three-dimensional Helmholtz representation")
    if (.not. all_passed) error stop 1

contains

    pure function cross_product(first, second) result(product)
        real(dp), intent(in) :: first(3), second(3)
        real(dp) :: product(3)

        product = [ &
            first(2)*second(3) - first(3)*second(2), &
            first(3)*second(1) - first(1)*second(3), &
            first(1)*second(2) - first(2)*second(1)]
    end function cross_product

    subroutine record_condition(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        call check_condition(condition, description)
        all_passed = all_passed .and. condition
    end subroutine record_condition
end program test_helmholtz_representation_3d
