module meridian_fourier_apply_fields
    !! Phi-dependent symmetric positive definite tensor (Fourier degree <= 2).
    use fortfem_kinds, only: dp
    use fortfem_meridian_fourier_scalar, only: meridian_fourier_tensor_t
    implicit none
    private
    public :: apply_tensor_t

    type, extends(meridian_fourier_tensor_t) :: apply_tensor_t
    contains
        procedure :: tensor => apply_eval
    end type

contains

    subroutine apply_eval(self, t, x, k, status)
        class(apply_tensor_t), intent(in) :: self
        integer, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: k(:, :, :)
        integer, intent(out) :: status
        integer :: p
        do p = 1, size(x, 2)
            k(:, :, p) = 0.0_dp
            k(1, 1, p) = 2.0_dp + 0.2_dp*x(1, p) + 0.5_dp*cos(x(3, p))
            k(2, 2, p) = 1.5_dp + 0.1_dp*x(2, p) + 0.4_dp*cos(2*x(3, p))
            k(3, 3, p) = 1.0_dp + 0.1_dp*x(1, p) + 0.3_dp*sin(x(3, p))
            k(1, 2, p) = 0.3_dp*sin(x(3, p)); k(2, 1, p) = k(1, 2, p)
            k(1, 3, p) = 0.1_dp + 0.05_dp*x(1, p); k(3, 1, p) = k(1, 3, p)
            k(2, 3, p) = 0.2_dp*cos(x(3, p)); k(3, 2, p) = k(2, 3, p)
        end do
        status = 0
    end subroutine apply_eval

end module meridian_fourier_apply_fields

program test_meridian_fourier_apply
    !! The matrix-free FFT apply equals the assembled mode-coupled operator
    !! (mass on V1 and V2, curl-curl) for a random vector, to roundoff, for
    !! even and odd phi grids.
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use fortfem_mesh_2d, only: mesh_2d_t
    use fortfem_meridian_curved_mesh, only: meridian_mesh_t, &
        meridian_mesh_from_map
    use fortfem_meridian_feec_n, only: meridian_feec_t, meridian_feec_build
    use fortfem_meridian_fourier_apply, only: meridian_fourier_op_t, &
        meridian_fourier_op_init, meridian_fourier_mass_apply, &
        meridian_fourier_curl_curl_apply
    use fortfem_meridian_fourier_curl, only: meridian_fourier_mass, &
        meridian_fourier_curl_curl, meridian_fourier_curl
    use fortsparse, only: csc_t, csc_matvec, csc_transpose, fortsparse_status_t
    use meridian_fourier_apply_fields, only: apply_tensor_t
    implicit none

    type(mesh_2d_t) :: pmesh, topo
    type(meridian_mesh_t) :: mm
    type(meridian_feec_t) :: fs
    type(apply_tensor_t) :: coef
    integer :: st, t, case, space
    integer, parameter :: modes_c(3) = [3, 3, 5], nphi_c(3) = [16, 7, 24]
    real(dp), allocatable :: pv(:, :, :)

    call pmesh%create_rectangular(5, 4, 2.0_dp, 3.0_dp, -0.5_dp, 0.5_dp)
    topo = pmesh
    allocate (pv(2, 3, topo%n_triangles))
    do t = 1, topo%n_triangles
        pv(:, :, t) = pmesh%vertices(:, pmesh%triangles(:, t))
    end do
    call meridian_mesh_from_map(topo, pv, 1, mm, st, ident_map)
    call check_condition(st == 0, "build rectangle mesh")
    call meridian_feec_build(mm, 2, fs, st)
    call check_condition(st == 0, "build P2 complex")

    do case = 1, 3
        do space = 1, 2
            call compare_mass(space, modes_c(case), nphi_c(case))
        end do
        call compare_curl_curl(modes_c(case), nphi_c(case))
    end do
    call check_summary("meridian_fourier_apply")

contains

    subroutine ident_map(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x(2)
        x = p
    end subroutine ident_map

    subroutine compare_mass(space, modes, nphi)
        integer, intent(in) :: space, modes, nphi
        type(csc_t) :: m
        type(meridian_fourier_op_t) :: op
        real(dp), allocatable :: u(:), ya(:), ym(:)
        integer :: st
        character(len=64) :: label

        write (label, "(a,i0,a,i0,a,i0)") "mass space ", space, " N ", modes, &
            " nphi ", nphi
        call meridian_fourier_mass(fs, space, modes, nphi, coef, m, st)
        call check_condition(st == 0, "assembled "//trim(label))
        call meridian_fourier_op_init(fs, space, modes, nphi, coef, op, st)
        call check_condition(st == 0, "apply setup "//trim(label))
        allocate (u(m%ncol), ya(m%ncol), ym(m%ncol))
        call random_number(u)
        u = u - 0.5_dp
        ym = csc_matvec(m, u)
        call meridian_fourier_mass_apply(op, u, ya)
        call check_condition(maxval(abs(ya - ym)) <= 1.0e-12_dp*maxval(abs(ym)), &
            "apply equals assembled "//trim(label))
    end subroutine compare_mass

    subroutine compare_curl_curl(modes, nphi)
        integer, intent(in) :: modes, nphi
        type(csc_t) :: s, c, ct
        type(fortsparse_status_t) :: fst
        type(meridian_fourier_op_t) :: op
        real(dp), allocatable :: u(:), ya(:), ym(:)
        integer :: st
        character(len=64) :: label

        write (label, "(a,i0,a,i0)") "curl-curl N ", modes, " nphi ", nphi
        call meridian_fourier_curl_curl(fs, modes, nphi, coef, s, st)
        call check_condition(st == 0, "assembled "//trim(label))
        call meridian_fourier_curl(fs, modes, c, st)
        call csc_transpose(c, ct, fst)
        call meridian_fourier_op_init(fs, 2, modes, nphi, coef, op, st)
        allocate (u(s%ncol), ya(s%ncol), ym(s%ncol))
        call random_number(u)
        u = u - 0.5_dp
        ym = csc_matvec(s, u)
        call meridian_fourier_curl_curl_apply(op, c, ct, u, ya)
        call check_condition(maxval(abs(ya - ym)) <= 1.0e-11_dp*maxval(abs(ym)), &
            "apply equals assembled "//trim(label))
    end subroutine compare_curl_curl

end program test_meridian_fourier_apply
