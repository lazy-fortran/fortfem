module fortfem_gauss_quadrature_2d
    use fortfem_kinds, only: dp
    implicit none
    private

    public :: gauss_quadrature_triangle_t
    public :: get_gauss_quadrature_triangle

    type :: gauss_quadrature_triangle_t
        integer :: n_points
        real(dp), allocatable :: xi(:)
        real(dp), allocatable :: eta(:)
        real(dp), allocatable :: weights(:)
    contains
        procedure :: init
        procedure :: destroy
    end type

contains

    subroutine init(this, order)
        class(gauss_quadrature_triangle_t), intent(out) :: this
        integer, intent(in) :: order

        select case(order)
        case(1)
            call init_order_1(this)
        case(2)
            call init_order_2(this)
        case(3)
            call init_order_3(this)
        case(4)
            call init_order_4(this)
        case(5)
            call init_order_5(this)
        case(6)
            call init_order_6(this)
        case(7)
            call init_order_7(this)
        case default
            print *, "Error: Gauss quadrature order", order, "not implemented"
            print *, "Available orders: 1-7"
            stop 1
        end select
    end subroutine

    subroutine destroy(this)
        class(gauss_quadrature_triangle_t), intent(inout) :: this

        if (allocated(this%xi)) deallocate(this%xi)
        if (allocated(this%eta)) deallocate(this%eta)
        if (allocated(this%weights)) deallocate(this%weights)
        this%n_points = 0
    end subroutine

    ! Get quadrature rule based on polynomial degree to integrate exactly
    function get_gauss_quadrature_triangle(degree) result(quad)
        integer, intent(in) :: degree
        type(gauss_quadrature_triangle_t) :: quad
        integer :: order

        ! Map polynomial degree to quadrature order
        ! order p integrates polynomials of degree 2p-1 exactly
        if (degree <= 1) then
            order = 1
        else if (degree <= 2) then
            order = 2
        else if (degree <= 3) then
            order = 3
        else if (degree <= 4) then
            order = 4
        else if (degree <= 5) then
            order = 5
        else if (degree <= 6) then
            order = 6
        else
            order = 7
        end if

        call quad%init(order)
    end function

    ! Order 1: degree 1 exact, 1 point (centroid)
    subroutine init_order_1(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this

        this%n_points = 1
        allocate(this%xi(1), this%eta(1), this%weights(1))

        this%xi(1) = 1.0_dp/3.0_dp
        this%eta(1) = 1.0_dp/3.0_dp
        this%weights(1) = 0.5_dp ! Area of reference triangle
    end subroutine

    ! Order 2: degree 2 exact, 3 points (edge midpoints)
    subroutine init_order_2(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this

        this%n_points = 3
        allocate(this%xi(3), this%eta(3), this%weights(3))

        this%xi = [0.5_dp, 0.0_dp, 0.5_dp]
        this%eta = [0.0_dp, 0.5_dp, 0.5_dp]
        this%weights = [1.0_dp/6.0_dp, 1.0_dp/6.0_dp, 1.0_dp/6.0_dp]
    end subroutine

    ! Order 3: degree 3 exact, 4 points
    subroutine init_order_3(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this

        this%n_points = 4
        allocate(this%xi(4), this%eta(4), this%weights(4))

        ! Centroid point
        this%xi(1) = 1.0_dp/3.0_dp
        this%eta(1) = 1.0_dp/3.0_dp
        this%weights(1) = -27.0_dp/96.0_dp

        ! Three points on edges
        this%xi(2) = 0.6_dp
        this%eta(2) = 0.2_dp
        this%weights(2) = 25.0_dp/96.0_dp

        this%xi(3) = 0.2_dp
        this%eta(3) = 0.6_dp
        this%weights(3) = 25.0_dp/96.0_dp

        this%xi(4) = 0.2_dp
        this%eta(4) = 0.2_dp
        this%weights(4) = 25.0_dp/96.0_dp
    end subroutine

    ! Order 4: degree 4 exact, 6 points
    subroutine init_order_4(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this
        real(dp), parameter :: a = 0.445948490915965_dp
        real(dp), parameter :: b = 0.091576213509771_dp
        real(dp), parameter :: w1 = 0.111690794839005_dp
        real(dp), parameter :: w2 = 0.054975871827661_dp

        this%n_points = 6
        allocate(this%xi(6), this%eta(6), this%weights(6))

        ! Three points near vertices
        this%xi(1) = a
        this%eta(1) = 1.0_dp - 2.0_dp*a
        this%weights(1) = w1

        this%xi(2) = a
        this%eta(2) = a
        this%weights(2) = w1

        this%xi(3) = 1.0_dp - 2.0_dp*a
        this%eta(3) = a
        this%weights(3) = w1

        ! Three points near edge midpoints
        this%xi(4) = b
        this%eta(4) = 1.0_dp - 2.0_dp*b
        this%weights(4) = w2

        this%xi(5) = b
        this%eta(5) = b
        this%weights(5) = w2

        this%xi(6) = 1.0_dp - 2.0_dp*b
        this%eta(6) = b
        this%weights(6) = w2
    end subroutine

    ! Order 5: degree 5 exact, 7 points (Dunavant)
    subroutine init_order_5(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this
        real(dp), parameter :: a1 = 0.101286507323456_dp
        real(dp), parameter :: a2 = 0.797426985353087_dp
        real(dp), parameter :: b1 = 0.470142064105115_dp
        real(dp), parameter :: b2 = 0.059715871789770_dp
        real(dp), parameter :: w0 = 0.225000000000000_dp
        real(dp), parameter :: w1 = 0.125939180544827_dp
        real(dp), parameter :: w2 = 0.132394152788506_dp

        this%n_points = 7
        allocate(this%xi(7), this%eta(7), this%weights(7))

        ! Centroid
        this%xi(1) = 1.0_dp/3.0_dp
        this%eta(1) = 1.0_dp/3.0_dp
        this%weights(1) = w0 * 0.5_dp

        ! First group: three points
        this%xi(2) = a2
        this%eta(2) = a1
        this%weights(2) = w1 * 0.5_dp

        this%xi(3) = a1
        this%eta(3) = a2
        this%weights(3) = w1 * 0.5_dp

        this%xi(4) = a1
        this%eta(4) = a1
        this%weights(4) = w1 * 0.5_dp

        ! Second group: three points
        this%xi(5) = b2
        this%eta(5) = b1
        this%weights(5) = w2 * 0.5_dp

        this%xi(6) = b1
        this%eta(6) = b2
        this%weights(6) = w2 * 0.5_dp

        this%xi(7) = b1
        this%eta(7) = b1
        this%weights(7) = w2 * 0.5_dp
    end subroutine

    ! Order 6: degree 6 exact, 12 points
    subroutine init_order_6(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this
        real(dp), parameter :: a = 0.063089014491502_dp
        real(dp), parameter :: b = 0.249286745170910_dp
        real(dp), parameter :: c = 0.310352451033785_dp
        real(dp), parameter :: d = 0.053145049844816_dp
        real(dp), parameter :: w1 = 0.050844906370207_dp
        real(dp), parameter :: w2 = 0.116786275726379_dp
        real(dp), parameter :: w3 = 0.082851075618374_dp

        this%n_points = 12
        allocate(this%xi(12), this%eta(12), this%weights(12))

        ! First group: 3 points
        this%xi(1) = a
        this%eta(1) = a
        this%weights(1) = w1 * 0.5_dp

        this%xi(2) = 1.0_dp - 2.0_dp*a
        this%eta(2) = a
        this%weights(2) = w1 * 0.5_dp

        this%xi(3) = a
        this%eta(3) = 1.0_dp - 2.0_dp*a
        this%weights(3) = w1 * 0.5_dp

        ! Second group: 3 points
        this%xi(4) = b
        this%eta(4) = b
        this%weights(4) = w2 * 0.5_dp

        this%xi(5) = 1.0_dp - 2.0_dp*b
        this%eta(5) = b
        this%weights(5) = w2 * 0.5_dp

        this%xi(6) = b
        this%eta(6) = 1.0_dp - 2.0_dp*b
        this%weights(6) = w2 * 0.5_dp

        ! Third group: 6 points
        this%xi(7) = c
        this%eta(7) = d
        this%weights(7) = w3 * 0.5_dp

        this%xi(8) = d
        this%eta(8) = c
        this%weights(8) = w3 * 0.5_dp

        this%xi(9) = 1.0_dp - c - d
        this%eta(9) = c
        this%weights(9) = w3 * 0.5_dp

        this%xi(10) = 1.0_dp - c - d
        this%eta(10) = d
        this%weights(10) = w3 * 0.5_dp

        this%xi(11) = c
        this%eta(11) = 1.0_dp - c - d
        this%weights(11) = w3 * 0.5_dp

        this%xi(12) = d
        this%eta(12) = 1.0_dp - c - d
        this%weights(12) = w3 * 0.5_dp
    end subroutine

    ! Order 7: degree 7 exact, 13 points
    subroutine init_order_7(this)
        type(gauss_quadrature_triangle_t), intent(out) :: this
        real(dp), parameter :: a = 0.0651301029022_dp
        real(dp), parameter :: b = 0.2603459660790_dp
        real(dp), parameter :: c = 0.0839477740996_dp
        real(dp), parameter :: d = 0.3128654960049_dp
        real(dp), parameter :: e = 0.0486903154253_dp
        real(dp), parameter :: w0 = -0.1495700444677_dp
        real(dp), parameter :: w1 = 0.0533472356089_dp
        real(dp), parameter :: w2 = 0.0771137608903_dp
        real(dp), parameter :: w3 = 0.1756152574332_dp

        this%n_points = 13
        allocate(this%xi(13), this%eta(13), this%weights(13))

        ! Centroid
        this%xi(1) = 1.0_dp/3.0_dp
        this%eta(1) = 1.0_dp/3.0_dp
        this%weights(1) = w0 * 0.5_dp

        ! First group: 3 points
        this%xi(2) = a
        this%eta(2) = a
        this%weights(2) = w1 * 0.5_dp

        this%xi(3) = 1.0_dp - 2.0_dp*a
        this%eta(3) = a
        this%weights(3) = w1 * 0.5_dp

        this%xi(4) = a
        this%eta(4) = 1.0_dp - 2.0_dp*a
        this%weights(4) = w1 * 0.5_dp

        ! Second group: 3 points
        this%xi(5) = b
        this%eta(5) = b
        this%weights(5) = w2 * 0.5_dp

        this%xi(6) = 1.0_dp - 2.0_dp*b
        this%eta(6) = b
        this%weights(6) = w2 * 0.5_dp

        this%xi(7) = b
        this%eta(7) = 1.0_dp - 2.0_dp*b
        this%weights(7) = w2 * 0.5_dp

        ! Third group: 6 points
        this%xi(8) = c
        this%eta(8) = d
        this%weights(8) = w3 * 0.5_dp

        this%xi(9) = d
        this%eta(9) = c
        this%weights(9) = w3 * 0.5_dp

        this%xi(10) = 1.0_dp - c - d
        this%eta(10) = c
        this%weights(10) = w3 * 0.5_dp

        this%xi(11) = 1.0_dp - c - d
        this%eta(11) = d
        this%weights(11) = w3 * 0.5_dp

        this%xi(12) = c
        this%eta(12) = 1.0_dp - c - d
        this%weights(12) = w3 * 0.5_dp

        this%xi(13) = d
        this%eta(13) = 1.0_dp - c - d
        this%weights(13) = w3 * 0.5_dp
    end subroutine

end module fortfem_gauss_quadrature_2d
