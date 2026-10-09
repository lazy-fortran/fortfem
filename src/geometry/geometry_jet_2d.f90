module fortfem_geometry_jet_2d
    !! Transform a scalar jet using caller-supplied sampled geometry.
    !! J(a,i) = dT_a/dxi_i and map_hessian(i,j,a) = d2T_a/dxi_i dxi_j.
    !! Success checks finite inputs/results and a positive computed determinant.
    !! It does not certify global invertibility, conformity, conditioning or roundoff.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: geometry_jet_success = 0
    integer, parameter, public :: geometry_jet_invalid_input = 1
    integer, parameter, public :: geometry_jet_nonpositive_jacobian = 2
    integer, parameter, public :: geometry_jet_nonfinite_result = 3
    public :: transform_scalar_jet_2d

contains

    pure subroutine transform_scalar_jet_2d(point, jacobian, map_hessian, &
            reference_gradient, reference_hessian, determinant, &
            inverse_jacobian, gradient, hessian, status)
        real(dp), intent(in) :: point(2), jacobian(2, 2), map_hessian(2, 2, 2)
        real(dp), intent(in) :: reference_gradient(2), reference_hessian(2, 2)
        real(dp), intent(out) :: determinant, inverse_jacobian(2, 2)
        real(dp), intent(out) :: gradient(2), hessian(2, 2)
        integer, intent(out) :: status
        real(dp) :: det, inverse(2, 2), physical_gradient(2)
        real(dp) :: corrected_hessian(2, 2), physical_hessian(2, 2)

        determinant = 0.0_dp
        inverse_jacobian = 0.0_dp
        gradient = 0.0_dp
        hessian = 0.0_dp
        status = geometry_jet_invalid_input
        if (.not. all(ieee_is_finite(point))) return
        if (.not. all(ieee_is_finite(jacobian))) return
        if (.not. all(ieee_is_finite(map_hessian))) return
        if (.not. all(ieee_is_finite(reference_gradient))) return
        if (.not. all(ieee_is_finite(reference_hessian))) return

        status = geometry_jet_nonfinite_result
        det = jacobian(1, 1)*jacobian(2, 2) &
            - jacobian(1, 2)*jacobian(2, 1)
        if (.not. ieee_is_finite(det)) return
        if (det <= 0.0_dp) then
            status = geometry_jet_nonpositive_jacobian
            return
        end if
        inverse(1, 1) = jacobian(2, 2)/det
        inverse(1, 2) = -jacobian(1, 2)/det
        inverse(2, 1) = -jacobian(2, 1)/det
        inverse(2, 2) = jacobian(1, 1)/det
        if (.not. all(ieee_is_finite(inverse))) return

        physical_gradient = matmul(transpose(inverse), reference_gradient)
        if (.not. all(ieee_is_finite(physical_gradient))) return
        ! Differentiate grad_x u: the map curvature is essential even for P2.
        corrected_hessian = reference_hessian &
            - physical_gradient(1)*map_hessian(:, :, 1) &
            - physical_gradient(2)*map_hessian(:, :, 2)
        if (.not. all(ieee_is_finite(corrected_hessian))) return
        physical_hessian = matmul(transpose(inverse), &
            matmul(corrected_hessian, inverse))
        if (.not. all(ieee_is_finite(physical_hessian))) return

        determinant = det
        inverse_jacobian = inverse
        gradient = physical_gradient
        hessian = physical_hessian
        status = geometry_jet_success
    end subroutine transform_scalar_jet_2d

end module fortfem_geometry_jet_2d
