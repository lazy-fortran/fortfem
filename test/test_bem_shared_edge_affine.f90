module shared_edge_affine_test_geometry
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    implicit none
    type, extends(bem_surface_t) :: rectangle_pair_t
        real(dp) :: vertices(3,3,2)
    contains
        procedure :: evaluate
    end type
contains
    pure subroutine evaluate(self,panel,uv,x,a,b,status)
        class(rectangle_pair_t), intent(in) :: self
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        real(dp), intent(out) :: x(3),a(3),b(3)
        integer, intent(out) :: status
        status=1
        x=0; a=0; b=0
        if (panel<1 .or. panel>2) return
        a=self%vertices(:,2,panel)-self%vertices(:,1,panel)
        b=self%vertices(:,3,panel)-self%vertices(:,1,panel)
        x=self%vertices(:,1,panel)+uv(1)*a+uv(2)*b
        status=0
    end subroutine
end module


program test_bem_shared_edge_affine
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use shared_edge_affine_test_geometry, only: rectangle_pair_t
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
        bem_pair_report_t, bem_shared_edge_chart_t, &
        integrate_parametric_bem_p2_pair
    implicit none
    ! Independent physical clipped-polygon quadratic overlap moments,
    ! integrated over displacement-polar coordinates at angular order256.
    ! Target [(0,0),(40,0),(0,1)], source [(40,0),(40,1),(0,1)].
    ! G=1/(4*pi*r); P1 flux and P2 trace, original native vertex order.
    real(dp), parameter :: reference_v(3,3)=reshape([ &
        3.51797794581342638e-01_dp, 7.63937288026291350e-01_dp , &
        3.78370791512707028e-01_dp, 3.34011981834287597e-01_dp , &
        6.75654411141379807e-01_dp, 3.51797794581346523e-01_dp , &
        6.75654411141376809e-01_dp, 8.85988866843433032e-01_dp , &
        7.63937288026289019e-01_dp], [3,3])
    real(dp), parameter :: reference_w(6,6)=reshape([ &
        -6.12364905719449126e-01_dp, 0.00000000000000000e+00_dp , &
        7.17763701070968985e-01_dp, 2.92067434437624485e+00_dp , &
        -2.92067434437624485e+00_dp, -1.05398795351511865e-01_dp , &
        3.66600145396174304e-01_dp, 7.77151369798792513e-04_dp , &
        -6.12364905719426922e-01_dp, -1.50915444756459238e+00_dp , &
        1.50793102747816921e+00_dp, 2.46211029039886675e-01_dp , &
        7.77151369798792621e-04_dp, -4.70204980684104344e-04_dp , &
        0.00000000000000000e+00_dp, -3.06946389114688169e-04_dp , &
        -3.90410819496203425e-03_dp, 3.90410819496203425e-03_dp , &
        2.46211029039899137e-01_dp, 3.90410819496205853e-03_dp , &
        -1.05398795351542063e-01_dp, -1.41664742509303054e+00_dp , &
        1.41630973807370819e+00_dp, -1.44378654864014466e-01_dp , &
        -1.50915444756457284e+00_dp, -3.06946389114687301e-04_dp , &
        2.92067434437622886e+00_dp, -1.41745087897588782e+01_dp , &
        1.41799432644293777e+01_dp, -1.41664742509304475e+00_dp , &
        1.50793102747814967e+00_dp, -3.90410819496205853e-03_dp , &
        -2.92067434437622886e+00_dp, 1.41799432644293777e+01_dp , &
        -1.41796055774100562e+01_dp, 1.41630973807372240e+00_dp], [6,6])
    real(dp), parameter :: reference_sum=5.18115062768845380e+00_dp
    type(rectangle_pair_t) :: pair
    type(bem_shared_edge_chart_t) :: chart,invalid_chart
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    real(dp) :: v(3,3),k(3,6),w(6,6),correction(3,6),errors(3)
    integer :: status
    pair%vertices(:,:,1)=reshape([0.0_dp,0.0_dp,0.0_dp, &
        40.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,3])
    pair%vertices(:,:,2)=reshape([40.0_dp,0.0_dp,0.0_dp, &
        40.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,3])
    options%relative_tolerance=1.e-8_dp
    options%maximum_degree=64
    ! Native atol1e-9, degree8:4:64 and work100M remain unchanged.
    chart%vertices(:,1)=[2,3,1]
    chart%vertices(:,2)=[1,3,2]
    chart%aspect=40
    invalid_chart=chart
    invalid_chart%vertices(:,1)=[2,2,1]
    call integrate_parametric_bem_p2_pair(pair,1,2,options,v,k,w,correction, &
        report,status,edge_chart=invalid_chart)
    call check_condition(status==1,'reject duplicate chart vertices')
    invalid_chart=chart
    invalid_chart%aspect=-1
    call integrate_parametric_bem_p2_pair(pair,1,2,options,v,k,w,correction, &
        report,status,edge_chart=invalid_chart)
    call check_condition(status==1,'reject nonpositive chart aspect')
    call integrate_parametric_bem_p2_pair(pair,1,1,options,v,k,w,correction, &
        report,status,edge_chart=chart)
    call check_condition(status==1,'reject self pair edge chart')
    call integrate_parametric_bem_p2_pair(pair,1,2,options,v,k,w,correction, &
        report,status,edge_chart=chart)
    call check_condition(status==0.and.report%accepted_pairs==1.and. &
        report%error_ratio<=1,'original native stopping gate')
    call check_condition(report%maximum_degree_used<=64.and. &
        report%kernel_evaluations+1848000<=options%maximum_evaluations, &
        'original degree and total work including moment controls')
    errors(1)=maxval(abs(v-reference_v))/maxval(abs(reference_v))
    errors(2)=maxval(abs(w-reference_w))/maxval(abs(reference_w))
    errors(3)=abs(sum(v)/reference_sum-1)
    call check_condition(all(ieee_is_finite(errors)).and.all(errors<1.e-7_dp), &
        'independent full V W and scalar overlap oracle')
    call check_condition(max(maxval(abs(k)),maxval(abs(correction)))<1.e-13_dp, &
        'coplanar K C oracle')
    call check_summary('explicit affine shared-edge pair')
end program
