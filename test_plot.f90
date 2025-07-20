program test_plot
    use fortplot
    implicit none
    
    real(8) :: x(4), y(4)
    
    x = [0.0d0, 1.0d0, 1.0d0, 0.0d0]
    y = [0.0d0, 0.0d0, 1.0d0, 1.0d0]
    
    call figure()
    call plot(x, y)
    call savefig("test.png")
    
    print *, "Test completed"
end program test_plot