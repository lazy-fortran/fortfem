! Extracted from src/generated/fortfem_tensor_power_split_products.f90 by
! tools/fortad/extract.py. It is the primal fortsym generated
! its products from, with the offload directives removed:
! they are the consumer's choice, not part of what the
! function means.

pure subroutine generated_tensor_power_split(tensor_11, tensor_12, tensor_13, tensor_21, tensor_22, &
        tensor_23, tensor_31, tensor_32, tensor_33, vector_1, vector_2, vector_3, symmetric_power, &
        skew_power, total_power)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: tensor_11, tensor_12, tensor_13, tensor_21, tensor_22, tensor_23, &
        tensor_31, tensor_32, tensor_33, vector_1, vector_2, vector_3
    real(dp), intent(out) :: symmetric_power, skew_power, total_power
    real(dp) :: t1, t2, t3

    t1 = vector_1**2
    t2 = vector_2**2
    t3 = vector_3**2
    symmetric_power = tensor_11*t1*2*5.0000000000000000E-001_dp + tensor_22*t2*2* &
        5.0000000000000000E-001_dp + tensor_33*t3*2*5.0000000000000000E-001_dp + vector_1*vector_2* &
        (tensor_12 + tensor_21)*2*5.0000000000000000E-001_dp + vector_1*vector_3*(tensor_13 + tensor_31)*2* &
        5.0000000000000000E-001_dp + vector_2*vector_3*(tensor_23 + tensor_32)*2* &
        5.0000000000000000E-001_dp + 0.0000000000000000E+000_dp
    skew_power = vector_1*vector_2*(tensor_12 - tensor_21)*5.0000000000000000E-001_dp + vector_1* &
        vector_2*(tensor_21 - tensor_12)*5.0000000000000000E-001_dp + vector_1*vector_3*(tensor_13 - &
        tensor_31)*5.0000000000000000E-001_dp + vector_1*vector_3*(tensor_31 - tensor_13)* &
        5.0000000000000000E-001_dp + vector_2*vector_3*(tensor_23 - tensor_32)*5.0000000000000000E-001_dp + &
        vector_2*vector_3*(tensor_32 - tensor_23)*5.0000000000000000E-001_dp + 0.0000000000000000E+000_dp
    total_power = tensor_11*t1 + tensor_12*vector_1*vector_2 + tensor_13*vector_1*vector_3 + &
        tensor_21*vector_1*vector_2 + tensor_22*t2 + tensor_23*vector_2*vector_3 + tensor_31*vector_1* &
        vector_3 + tensor_32*vector_2*vector_3 + tensor_33*t3 + 0.0000000000000000E+000_dp

end subroutine generated_tensor_power_split
