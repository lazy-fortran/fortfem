module fortfem_direct_fourier_transform
    !! Small, deterministic direct Fourier transform contract.
    !!
    !! The plan stores a caller-selected set of integer modes and sample
    !! angles.  It performs the normalized dense transform directly, with no
    !! FFT assumptions and no plasma-specific coordinate convention:
    !!
    !!   c_m = N**(-1/2) sum_j exp(-i*m*theta_j) f_j.
    !!
    !! The adjoint uses the conjugate transpose of the same matrix.  A square
    !! uniform plan with modes 0,...,N-1 therefore has an exact deterministic
    !! round trip up to floating-point arithmetic.  Chunk bounds are metadata
    !! only; they make bounded execution explicit for future large plans.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortsparse, only: fortsparse_status_t, status_set, &
        FORTSPARSE_INVALID_MATRIX, FORTSPARSE_OK
    implicit none
    private

    type, public :: direct_fourier_plan_t
        private
        integer :: sample_count = 0
        integer :: mode_count = 0
        integer :: chunk_size = 0
        integer :: chunk_count = 0
        real(dp), allocatable :: angles(:)
        real(dp), allocatable :: sample_weights(:)
        integer, allocatable :: mode_indices(:)
        integer, allocatable :: chunk_start(:)
        integer, allocatable :: chunk_end(:)
        complex(dp), allocatable :: kernel(:, :)
    contains
        procedure, private :: assign_direct_fourier_plan
        generic :: assignment(=) => assign_direct_fourier_plan
    end type direct_fourier_plan_t

    public :: initialize_direct_fourier_plan
    public :: validate_direct_fourier_plan
    public :: direct_fourier_plan_metadata
    public :: direct_fourier_plan_chunk_bounds
    public :: direct_fourier_forward
    public :: direct_fourier_adjoint
    public :: direct_fourier_forward_jvp
    public :: direct_fourier_forward_vjp
    public :: direct_fourier_plan_sample_count
    public :: direct_fourier_plan_mode_count

contains

    subroutine assign_direct_fourier_plan(lhs, rhs)
        class(direct_fourier_plan_t), intent(out) :: lhs
        type(direct_fourier_plan_t), intent(in) :: rhs

        lhs%sample_count = rhs%sample_count
        lhs%mode_count = rhs%mode_count
        lhs%chunk_size = rhs%chunk_size
        lhs%chunk_count = rhs%chunk_count
        if (allocated(rhs%angles)) allocate(lhs%angles, source=rhs%angles)
        if (allocated(rhs%sample_weights)) then
            allocate(lhs%sample_weights, source=rhs%sample_weights)
        end if
        if (allocated(rhs%mode_indices)) allocate(lhs%mode_indices, source=rhs%mode_indices)
        if (allocated(rhs%chunk_start)) allocate(lhs%chunk_start, source=rhs%chunk_start)
        if (allocated(rhs%chunk_end)) allocate(lhs%chunk_end, source=rhs%chunk_end)
        if (allocated(rhs%kernel)) allocate(lhs%kernel, source=rhs%kernel)
    end subroutine assign_direct_fourier_plan

    subroutine initialize_direct_fourier_plan( &
            plan, angles, mode_indices, chunk_size, sample_weights, status)
        type(direct_fourier_plan_t), intent(out) :: plan
        real(dp), intent(in) :: angles(:)
        integer, intent(in) :: mode_indices(:)
        integer, intent(in) :: chunk_size
        real(dp), intent(in), optional :: sample_weights(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: sample, mode, selected_size
        real(dp) :: phase, normalization

        call clear_direct_fourier_plan(plan)
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "direct Fourier plan received incompatible metadata")
        if (.not. valid_raw_plan(angles, mode_indices, chunk_size, sample_weights)) return
        plan%sample_count = size(angles)
        plan%mode_count = size(mode_indices)
        selected_size = chunk_size
        plan%chunk_size = selected_size
        plan%chunk_count = (plan%sample_count + selected_size - 1) / selected_size
        allocate(plan%angles(plan%sample_count), plan%sample_weights(plan%sample_count), &
            plan%mode_indices(plan%mode_count), plan%chunk_start(plan%chunk_count), &
            plan%chunk_end(plan%chunk_count), plan%kernel(plan%mode_count, plan%sample_count))
        plan%angles = angles
        plan%mode_indices = mode_indices
        if (present(sample_weights)) then
            plan%sample_weights = sample_weights
        else
            plan%sample_weights = 1.0_dp
        end if
        normalization = 1.0_dp/sqrt(real(plan%sample_count, dp))
        do mode = 1, plan%mode_count
            do sample = 1, plan%sample_count
                phase = real(plan%mode_indices(mode), dp)*plan%angles(sample)
                plan%kernel(mode, sample) = normalization*cmplx( &
                    cos(phase), -sin(phase), dp)
            end do
        end do
        do sample = 1, plan%chunk_count
            plan%chunk_start(sample) = (sample - 1)*selected_size + 1
            plan%chunk_end(sample) = min(sample*selected_size, plan%sample_count)
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine initialize_direct_fourier_plan

    logical function validate_direct_fourier_plan(plan, status) result(valid)
        type(direct_fourier_plan_t), intent(in) :: plan
        type(fortsparse_status_t), intent(out) :: status

        valid = .false.
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "direct Fourier plan has inconsistent metadata")
        if (.not. allocated(plan%angles) .or. .not. allocated(plan%mode_indices) .or. &
            .not. allocated(plan%sample_weights) .or. .not. allocated(plan%kernel) .or. &
            .not. allocated(plan%chunk_start) .or. .not. allocated(plan%chunk_end)) return
        if (.not. valid_raw_plan( &
            plan%angles, plan%mode_indices, plan%chunk_size, plan%sample_weights)) return
        if (plan%sample_count /= size(plan%angles) .or. &
            plan%mode_count /= size(plan%mode_indices) .or. &
            size(plan%kernel, 1) /= plan%mode_count .or. &
            size(plan%kernel, 2) /= plan%sample_count) return
        if (plan%chunk_count /= (plan%sample_count + plan%chunk_size - 1) / &
            plan%chunk_size .or. size(plan%chunk_start) /= plan%chunk_count .or. &
            size(plan%chunk_end) /= plan%chunk_count) return
        if (any(.not. ieee_is_finite(real(plan%kernel, dp))) .or. &
            any(.not. ieee_is_finite(aimag(plan%kernel)))) return
        if (plan%chunk_start(1) /= 1 .or. &
            plan%chunk_end(plan%chunk_count) /= plan%sample_count) return
        if (any(plan%chunk_start < 1) .or. any(plan%chunk_end > plan%sample_count) .or. &
            any(plan%chunk_start > plan%chunk_end) .or. &
            any(plan%chunk_end - plan%chunk_start + 1 > plan%chunk_size)) return
        if (plan%chunk_count > 1) then
            if (any(plan%chunk_start(2:) /= plan%chunk_end(:plan%chunk_count-1) + 1)) return
        end if
        valid = .true.
        call status_set(status, FORTSPARSE_OK, "")
    end function validate_direct_fourier_plan

    subroutine direct_fourier_plan_metadata( &
            plan, angles, mode_indices, sample_weights, chunk_size, chunk_start, &
            chunk_end, status)
        type(direct_fourier_plan_t), intent(in) :: plan
        real(dp), allocatable, intent(out) :: angles(:), sample_weights(:)
        integer, allocatable, intent(out) :: mode_indices(:), chunk_start(:), chunk_end(:)
        integer, intent(out) :: chunk_size
        type(fortsparse_status_t), intent(out) :: status

        chunk_size = 0
        if (allocated(angles)) deallocate(angles)
        if (allocated(sample_weights)) deallocate(sample_weights)
        if (allocated(mode_indices)) deallocate(mode_indices)
        if (allocated(chunk_start)) deallocate(chunk_start)
        if (allocated(chunk_end)) deallocate(chunk_end)
        if (.not. validate_direct_fourier_plan(plan, status)) return
        chunk_size = plan%chunk_size
        allocate(angles(size(plan%angles)), sample_weights(size(plan%sample_weights)), &
            mode_indices(size(plan%mode_indices)), chunk_start(size(plan%chunk_start)), &
            chunk_end(size(plan%chunk_end)))
        angles = plan%angles
        sample_weights = plan%sample_weights
        mode_indices = plan%mode_indices
        chunk_start = plan%chunk_start
        chunk_end = plan%chunk_end
    end subroutine direct_fourier_plan_metadata

    subroutine direct_fourier_plan_chunk_bounds( &
            plan, chunk, first_sample, last_sample, status)
        type(direct_fourier_plan_t), intent(in) :: plan
        integer, intent(in) :: chunk
        integer, intent(out) :: first_sample, last_sample
        type(fortsparse_status_t), intent(out) :: status

        first_sample = 0
        last_sample = 0
        call status_set(status, FORTSPARSE_INVALID_MATRIX, &
            "direct Fourier chunk is out of bounds")
        if (.not. validate_direct_fourier_plan(plan, status)) return
        if (chunk < 1 .or. chunk > plan%chunk_count) return
        first_sample = plan%chunk_start(chunk)
        last_sample = plan%chunk_end(chunk)
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine direct_fourier_plan_chunk_bounds

    subroutine direct_fourier_forward(plan, samples, coefficients, status)
        type(direct_fourier_plan_t), intent(in) :: plan
        complex(dp), intent(in) :: samples(:)
        complex(dp), intent(out) :: coefficients(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, sample

        coefficients = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_direct_fourier_plan(plan, status)) return
        if (size(samples) /= plan%sample_count .or. size(coefficients) /= plan%mode_count .or. &
            .not. finite_complex(samples)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "direct Fourier forward received incompatible vectors")
            return
        end if
        do mode = 1, plan%mode_count
            do sample = 1, plan%sample_count
                coefficients(mode) = coefficients(mode) + plan%kernel(mode, sample)* &
                    samples(sample)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine direct_fourier_forward

    subroutine direct_fourier_adjoint(plan, coefficients, samples, status)
        type(direct_fourier_plan_t), intent(in) :: plan
        complex(dp), intent(in) :: coefficients(:)
        complex(dp), intent(out) :: samples(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, sample

        samples = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_direct_fourier_plan(plan, status)) return
        if (size(coefficients) /= plan%mode_count .or. size(samples) /= plan%sample_count .or. &
            .not. finite_complex(coefficients)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "direct Fourier adjoint received incompatible vectors")
            return
        end if
        do sample = 1, plan%sample_count
            do mode = 1, plan%mode_count
                samples(sample) = samples(sample) + conjg(plan%kernel(mode, sample))* &
                    coefficients(mode)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine direct_fourier_adjoint

    subroutine direct_fourier_forward_jvp( &
            plan, samples, sample_tangent, angle_tangent, coefficient_tangent, status)
        !! Differentiate the fixed-mode direct transform with respect to its
        !! complex samples and real sample angles.
        type(direct_fourier_plan_t), intent(in) :: plan
        complex(dp), intent(in) :: samples(:), sample_tangent(:)
        real(dp), intent(in) :: angle_tangent(:)
        complex(dp), intent(out) :: coefficient_tangent(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, sample
        complex(dp) :: kernel_tangent

        coefficient_tangent = cmplx(0.0_dp, 0.0_dp, dp)
        if (.not. validate_direct_fourier_plan(plan, status)) return
        if (size(samples) /= plan%sample_count .or. &
            size(sample_tangent) /= plan%sample_count .or. &
            size(angle_tangent) /= plan%sample_count .or. &
            size(coefficient_tangent) /= plan%mode_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "direct Fourier JVP received incompatible vectors")
            return
        end if
        if (.not. finite_complex(samples) .or. &
            .not. finite_complex(sample_tangent) .or. &
            .not. all(ieee_is_finite(angle_tangent))) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "direct Fourier JVP received non-finite vectors")
            return
        end if
        do mode = 1, plan%mode_count
            do sample = 1, plan%sample_count
                kernel_tangent = plan%kernel(mode, sample)*cmplx( &
                    0.0_dp, -real(plan%mode_indices(mode), dp)* &
                    angle_tangent(sample), dp)
                coefficient_tangent(mode) = coefficient_tangent(mode) + &
                    plan%kernel(mode, sample)*sample_tangent(sample) + &
                    kernel_tangent*samples(sample)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine direct_fourier_forward_jvp

    subroutine direct_fourier_forward_vjp( &
            plan, samples, coefficient_cotangent, sample_cotangent, &
            angle_cotangent, status)
        !! Apply the real-part complex adjoint of the fixed-mode transform.
        type(direct_fourier_plan_t), intent(in) :: plan
        complex(dp), intent(in) :: samples(:), coefficient_cotangent(:)
        complex(dp), intent(out) :: sample_cotangent(:)
        real(dp), intent(out) :: angle_cotangent(:)
        type(fortsparse_status_t), intent(out) :: status

        integer :: mode, sample
        complex(dp) :: kernel_tangent

        sample_cotangent = cmplx(0.0_dp, 0.0_dp, dp)
        angle_cotangent = 0.0_dp
        if (.not. validate_direct_fourier_plan(plan, status)) return
        if (size(samples) /= plan%sample_count .or. &
            size(coefficient_cotangent) /= plan%mode_count .or. &
            size(sample_cotangent) /= plan%sample_count .or. &
            size(angle_cotangent) /= plan%sample_count) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "direct Fourier VJP received incompatible vectors")
            return
        end if
        if (.not. finite_complex(samples) .or. &
            .not. finite_complex(coefficient_cotangent)) then
            call status_set(status, FORTSPARSE_INVALID_MATRIX, &
                "direct Fourier VJP received non-finite vectors")
            return
        end if
        do sample = 1, plan%sample_count
            do mode = 1, plan%mode_count
                sample_cotangent(sample) = sample_cotangent(sample) + &
                    conjg(plan%kernel(mode, sample))*coefficient_cotangent(mode)
                kernel_tangent = plan%kernel(mode, sample)*cmplx( &
                    0.0_dp, -real(plan%mode_indices(mode), dp), dp)*samples(sample)
                angle_cotangent(sample) = angle_cotangent(sample) + real( &
                    conjg(coefficient_cotangent(mode))*kernel_tangent, dp)
            end do
        end do
        call status_set(status, FORTSPARSE_OK, "")
    end subroutine direct_fourier_forward_vjp

    integer function direct_fourier_plan_sample_count(plan)
        type(direct_fourier_plan_t), intent(in) :: plan

        direct_fourier_plan_sample_count = plan%sample_count
    end function direct_fourier_plan_sample_count

    integer function direct_fourier_plan_mode_count(plan)
        type(direct_fourier_plan_t), intent(in) :: plan

        direct_fourier_plan_mode_count = plan%mode_count
    end function direct_fourier_plan_mode_count

    logical function valid_raw_plan( &
            angles, mode_indices, chunk_size, sample_weights) result(valid)
        real(dp), intent(in) :: angles(:)
        integer, intent(in) :: mode_indices(:), chunk_size
        real(dp), intent(in), optional :: sample_weights(:)
        integer :: mode, other, sample, other_sample
        real(dp), parameter :: periodic_tolerance = 128.0_dp*epsilon(1.0_dp)

        valid = .false.
        if (size(angles) < 1 .or. size(mode_indices) < 1 .or. &
            chunk_size < 1 .or. chunk_size > size(angles) .or. &
            .not. all(ieee_is_finite(angles))) return
        if (present(sample_weights)) then
            if (size(sample_weights) /= size(angles) .or. &
                .not. all(ieee_is_finite(sample_weights)) .or. &
                any(sample_weights <= 0.0_dp)) return
        end if
        do sample = 1, size(angles)
            do other_sample = sample + 1, size(angles)
                if (abs(sin(0.5_dp*(angles(sample) - angles(other_sample)))) <= &
                    periodic_tolerance) return
            end do
        end do
        do mode = 1, size(mode_indices)
            do other = mode + 1, size(mode_indices)
                if (mode_indices(mode) == mode_indices(other)) return
            end do
        end do
        valid = .true.
    end function valid_raw_plan

    logical function finite_complex(values) result(valid)
        complex(dp), intent(in) :: values(:)

        valid = all(ieee_is_finite(real(values, dp))) .and. &
            all(ieee_is_finite(aimag(values)))
    end function finite_complex

    subroutine clear_direct_fourier_plan(plan)
        type(direct_fourier_plan_t), intent(out) :: plan

        plan%sample_count = 0
        plan%mode_count = 0
        plan%chunk_size = 0
        plan%chunk_count = 0
        if (allocated(plan%angles)) deallocate(plan%angles)
        if (allocated(plan%sample_weights)) deallocate(plan%sample_weights)
        if (allocated(plan%mode_indices)) deallocate(plan%mode_indices)
        if (allocated(plan%chunk_start)) deallocate(plan%chunk_start)
        if (allocated(plan%chunk_end)) deallocate(plan%chunk_end)
        if (allocated(plan%kernel)) deallocate(plan%kernel)
    end subroutine clear_direct_fourier_plan

end module fortfem_direct_fourier_transform
