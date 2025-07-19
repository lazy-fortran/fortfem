program basis_functions
    ! Visualize finite element basis functions
    use fortfem
    use fortplot
    implicit none
    
    real(dp) :: x(100), phi(100)
    integer :: i
    
    ! Plot P1 basis on [0,1]
    do i = 1, 100
        x(i) = real(i-1, dp) / 99.0_dp
        phi(i) = p1_basis_1d(x(i), node=1)  ! Hat function at x=0
    end do
    
    call figure()
    call plot(x, phi)
    call xlabel('x')
    call ylabel('phi(x)')
    call title('P1 Basis Function')
    call savefig('p1_basis.png')
    
end program basis_functions