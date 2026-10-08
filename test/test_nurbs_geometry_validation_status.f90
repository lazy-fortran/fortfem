program test_nurbs_geometry_validation_status
    use fortfem_kinds, only: dp
    use fortfem_bspline_feec, only: evaluate_nurbs_surface_geometry, &
        evaluate_nurbs_surface_geometry_jvp,evaluate_nurbs_surface_geometry_vjp, &
        evaluate_nurbs_volume_geometry,evaluate_nurbs_volume_geometry_jvp, &
        evaluate_nurbs_volume_geometry_vjp
    use check, only: check_condition,check_summary
    implicit none
    real(dp) :: surface(3,2,2),volume(3,2,2,2)
    real(dp), allocatable :: ws(:,:),wv(:,:,:)
    integer :: x,y,z,scenario
    do y=1,2
        do x=1,2
            surface(:,x,y)=[real(x-1,dp),real(y-1,dp),real(x+2*y-3,dp)]
            do z=1,2
                volume(:,x,y,z)=[real(x-1,dp),real(y-1,dp),real(z-1,dp)]
            end do
        end do
    end do
    do scenario=0,5
        if(scenario==4)then
            allocate(ws(1,2),wv(1,2,2))
        else
            allocate(ws(2,2),wv(2,2,2))
        end if
        ws=1;wv=1
        select case(scenario)
        case(1)
            ws=-1;wv=-1
        case(2)
            ws=0;wv=0
        case(3)
            ws=tiny(1.0_dp)/2;wv=tiny(1.0_dp)/2
        end select
        if(scenario==5)then
            call check_methods(surface(:2,:,:),volume(:2,:,:,:),ws,wv,.false.)
        else
            call check_methods(surface,volume,ws,wv,scenario==0)
        end if
        deallocate(ws,wv)
    end do
    call check_summary('NURBS geometry validation status and valid affine controls')
contains
    subroutine check_methods(s,v,weights_s,weights_v,valid)
        real(dp), intent(in) :: s(:,:,:),v(:,:,:,:),weights_s(:,:),weights_v(:,:,:)
        logical, intent(in) :: valid
        real(dp), parameter :: knots(4)=[0._dp,0._dp,1._dp,1._dp]
        real(dp), parameter :: translation(3)=[.07_dp,-.03_dp,.11_dp]
        real(dp), parameter :: covector(3)=[.2_dp,-.3_dp,.4_dp]
        real(dp) :: point(3),js(3,2),jv(3,3),expected_js(3,2),identity(3,3)
        real(dp) :: s_dot(size(s,1),size(s,2),size(s,3))
        real(dp) :: v_dot(size(v,1),size(v,2),size(v,3),size(v,4))
        real(dp) :: s_bar(size(s,1),size(s,2),size(s,3))
        real(dp) :: v_bar(size(v,1),size(v,2),size(v,3),size(v,4))
        real(dp) :: ws_dot(size(weights_s,1),size(weights_s,2)),ws_bar(size(weights_s,1),size(weights_s,2))
        real(dp) :: wv_dot(size(weights_v,1),size(weights_v,2),size(weights_v,3))
        real(dp) :: wv_bar(size(weights_v,1),size(weights_v,2),size(weights_v,3))
        real(dp) :: sum_bar(3)
        integer :: i,status
        expected_js=reshape([1._dp,0._dp,1._dp,0._dp,1._dp,2._dp],[3,2])
        identity=0
        do i=1,3
            identity(i,i)=1
        end do
        ws_dot=0;wv_dot=0
        do i=1,size(s,1)
            s_dot(i,:,:)=translation(i);v_dot(i,:,:,:)=translation(i)
        end do
        call evaluate_nurbs_surface_geometry(knots,knots,1,1,s,weights_s,.2_dp,.3_dp,point,js,status)
        call record(status,valid,'Surface primal')
        if(valid)then
            call check_condition(maxval(abs(point-[.2_dp,.3_dp,.8_dp]))<2e-14_dp.and. &
                maxval(abs(js-expected_js))<2e-14_dp,'Valid surface follows independent affine map')
        else
            call check_condition(all(point==0).and.all(js==0),'Invalid surface primal leaves zero outputs')
        end if
        call evaluate_nurbs_surface_geometry_jvp(knots,knots,1,1,s,weights_s,s_dot,ws_dot, &
            .2_dp,.3_dp,point,js,status)
        call record(status,valid,'Surface tangent')
        if(valid)then
            call check_condition(maxval(abs(point-translation))<2e-14_dp.and.maxval(abs(js))<2e-14_dp, &
                'Valid surface tangent reproduces rigid translation')
        else
            call check_condition(all(point==0).and.all(js==0),'Invalid surface tangent leaves zero outputs')
        end if
        js=0
        call evaluate_nurbs_surface_geometry_vjp(knots,knots,1,1,s,weights_s,.2_dp,.3_dp, &
            covector,js,s_bar,ws_bar,status)
        call record(status,valid,'Surface reverse')
        if(valid)then
            sum_bar=sum(sum(s_bar,dim=3),dim=2)
            call check_condition(maxval(abs(sum_bar-covector))<2e-14_dp, &
                'Valid surface reverse gives the rigid-translation covector')
        else
            call check_condition(all(s_bar==0).and.all(ws_bar==0),'Invalid surface reverse leaves zero outputs')
        end if
        call evaluate_nurbs_volume_geometry(knots,knots,knots,1,1,1,v,weights_v, &
            .2_dp,.3_dp,.4_dp,point,jv,status)
        call record(status,valid,'Volume primal')
        if(valid)then
            call check_condition(maxval(abs(point-[.2_dp,.3_dp,.4_dp]))<2e-14_dp.and. &
                maxval(abs(jv-identity))<2e-14_dp,'Valid volume follows independent affine map')
        else
            call check_condition(all(point==0).and.all(jv==0),'Invalid volume primal leaves zero outputs')
        end if
        call evaluate_nurbs_volume_geometry_jvp(knots,knots,knots,1,1,1,v,weights_v,v_dot,wv_dot, &
            .2_dp,.3_dp,.4_dp,point,jv,status)
        call record(status,valid,'Volume tangent')
        if(valid)then
            call check_condition(maxval(abs(point-translation))<2e-14_dp.and.maxval(abs(jv))<2e-14_dp, &
                'Valid volume tangent reproduces rigid translation')
        else
            call check_condition(all(point==0).and.all(jv==0),'Invalid volume tangent leaves zero outputs')
        end if
        jv=0
        call evaluate_nurbs_volume_geometry_vjp(knots,knots,knots,1,1,1,v,weights_v, &
            .2_dp,.3_dp,.4_dp,covector,jv,v_bar,wv_bar,status)
        call record(status,valid,'Volume reverse')
        if(valid)then
            sum_bar=sum(sum(sum(v_bar,dim=4),dim=3),dim=2)
            call check_condition(maxval(abs(sum_bar-covector))<2e-14_dp, &
                'Valid volume reverse gives the rigid-translation covector')
        else
            call check_condition(all(v_bar==0).and.all(wv_bar==0),'Invalid volume reverse leaves zero outputs')
        end if
    end subroutine check_methods
    subroutine record(status,valid,label)
        integer, intent(in) :: status
        logical, intent(in) :: valid
        character(*), intent(in) :: label
        if(valid)then
            call check_condition(status==0,label//' accepts valid geometry')
        else
            call check_condition(status/=0,label//' reports invalid shape, weight or denominator')
        end if
    end subroutine record
end program test_nurbs_geometry_validation_status
