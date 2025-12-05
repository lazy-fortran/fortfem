program test_vector_plotting
   use fortfem_kinds, only: dp
   use fortfem_api, only: mesh_t, vector_function_space_t, &
                          vector_trial_function_t, vector_test_function_t, vector_function_t, &
                          vector_bc_t, form_expr_t, unit_square_mesh, vector_function_space, &
                          vector_function, vector_trial_function, vector_test_function, &
                          vector_bc, inner, curl, dx, solve, plot, operator(*), operator(==), &
                          operator(+)
   use check, only: check_condition, check_summary
   implicit none

   call test_vector_streamplot_output()

   call check_summary("Vector plotting")

contains

   subroutine test_vector_streamplot_output()
      type(mesh_t) :: mesh
      type(vector_function_space_t) :: Vh
      type(vector_trial_function_t) :: E
      type(vector_test_function_t) :: F
      type(vector_function_t) :: J, Eh
      type(vector_bc_t) :: bc
      type(form_expr_t) :: a, L
      logical :: exists
      character(len=*), parameter :: filename = &
                                     "build/test_vector_plot.png"

      mesh = unit_square_mesh(4)
      Vh = vector_function_space(mesh, "Nedelec", 1)

      E = vector_trial_function(Vh)
      F = vector_test_function(Vh)
      J = vector_function(Vh)
      J%values(:, 1) = 1.0_dp
      J%values(:, 2) = 0.0_dp

      a = inner(curl(E), curl(F))*dx + inner(E, F)*dx
      L = inner(J, F)*dx

      bc = vector_bc(Vh, [0.0_dp, 0.0_dp], "tangential")
      Eh = vector_function(Vh)

      call solve(a == L, Eh, bc)

      call plot(Eh, filename=filename, title="Vector FEM Solution", &
                plot_type="streamplot")

      inquire (file=filename, exist=exists)
      call check_condition(exists, "Vector plot creates output file")
   end subroutine test_vector_streamplot_output

end program test_vector_plotting
