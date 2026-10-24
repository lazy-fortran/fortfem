module bem_pair_plane_geometry
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t
    implicit none
    type, extends(bem_surface_t) :: plane_t
        real(dp) :: vertices(3,3,2)
    contains
        procedure :: evaluate
    end type
contains
    pure subroutine evaluate(self,panel,uv,x,a,b,status)
        class(plane_t), intent(in) :: self
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
program test_parametric_bem_pair_accuracy
    use check, only: check_condition, check_summary
    use fortfem_kinds, only: dp
    use bem_pair_plane_geometry, only: plane_t
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
        bem_pair_report_t, integrate_parametric_bem_p2_pair
    implicit none
    type(plane_t) :: plane
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    real(dp) :: v(3,3),k(3,6),w(6,6),c(3,6),reference,refined,error,gap,aspect
    integer :: status,case_id
    logical :: passed
    passed=.true.
    plane%vertices(:,:,1)=reshape([0.0_dp,0.0_dp,0.0_dp, &
        1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,3])
    ! Exact right-triangle self energy: integrate the analytic edge potential.
    ! Integral_T Integral_T 1/|x-y| = (2+sqrt(2))*asinh(1)/3.
    ! This checks the raw form, before symmetry averaging or jump restoration.
    options%relative_tolerance=1.e-7_dp
    options%absolute_tolerance=1.e-12_dp
    options%maximum_degree=32
    reference=(2+sqrt(2.0_dp))*asinh(1.0_dp)/(12*acos(-1.0_dp))
    call integrate_parametric_bem_p2_pair(plane,1,1,options,v,k,w,c, &
        report,status)
    write(*,'(a,2i8,3es16.7)') 'BEM3D self status degree value exact ratio ', &
        status,report%maximum_degree_used,sum(v),reference,report%error_ratio
    call record(status==0,'coincident rule converges without raising degree cap')
    if (status==0) then
        call record(abs(sum(v)-reference)<1.e-8_dp*reference, &
            'raw self V agrees with exact triangle integral')
        call record(maxval(abs(v-transpose(v)))<1.e-11_dp, &
            'raw coincident V is symmetric')
        call record(maxval(abs(w-transpose(w)))<1.e-10_dp, &
            'raw coincident W is symmetric')
        call record(maxval(abs(k))+maxval(abs(c))<1.e-13_dp, &
            'coplanar coincident double layer vanishes')
    end if
    ! Skinny curved-surface panels also need the difference-angle peak resolved.
    aspect=40.0_dp
    plane%vertices(1,:,1)=aspect*plane%vertices(1,:,1)
    reference=aspect**2/(12*acos(-1.0_dp))* &
        ((asinh(aspect)+asinh(1/aspect))/sqrt(1+aspect**2)+ &
        asinh(aspect)/aspect+asinh(1/aspect))
    options%maximum_degree=48
    options%relative_tolerance=1.e-6_dp
    call integrate_parametric_bem_p2_pair(plane,1,1,options,v,k,w,c,report,status)
    write(*,'(a,2i8,3es16.7)') 'BEM3D skinny self status degree value exact ratio ', &
        status,report%maximum_degree_used,sum(v),reference,report%error_ratio
    call record(status==0,'skinny coincident panel converges at unchanged caps')
    if (status==0) call record(abs(sum(v)-reference)<1.e-7_dp*reference, &
        'skinny self V agrees with exact triangle integral')
    plane%vertices(1,:,1)=plane%vertices(1,:,1)/aspect
    options%initial_degree=8
    options%maximum_degree=64
    options%relative_tolerance=1.e-3_dp
    options%absolute_tolerance=1.e-9_dp
    do case_id=1,3
        gap=0
        if (case_id==1) gap=1.e-3_dp
        plane%vertices(:,:,2)=reshape([1.0_dp,-gap,0.0_dp, &
            0.0_dp,-gap,0.0_dp,1.0_dp,-1.0_dp-gap,0.0_dp],[3,3])
        if (case_id==3) plane%vertices(:,:,2)=reshape( &
            [0.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp,0.0_dp, &
            0.0_dp,-1.0_dp,0.0_dp],[3,3])
        ! Independent signed-edge analytic source potential, followed by
        ! target-only quadrature at two independently refined orders.
        call edge_reference(64,reference)
        call edge_reference(96,refined)
        call record(abs(reference-refined)<1.e-4_dp*abs(refined), &
            'independent analytic-edge reference is refined')
        call integrate_parametric_bem_p2_pair(plane,1,2,options,v,k,w,c, &
            report,status)
        write(*,'(a,4i8,4es16.7)') 'BEM3D pair case status degree evaluations ', &
            case_id,status,report%maximum_degree_used,report%kernel_evaluations, &
            sum(v),refined,report%error_ratio,abs(sum(v)-refined)/abs(refined)
        call record(status==0,'near/touching pair reaches requested tolerance')
        if (status/=0) cycle
        error=abs(sum(v)-refined)/abs(refined)
        call record(error<2.e-3_dp,'raw DG1 sum agrees with analytic-edge integral')
        call record(maxval(abs(k))<1.e-13_dp,'coplanar double layer vanishes')
    end do
    ! For parallel overlapping projections, D_n G develops a narrow peak.
    plane%vertices(:,:,2)=plane%vertices(:,:,1)
    plane%vertices(3,:,1)=0.02_dp
    options%relative_tolerance=3.e-3_dp
    options%maximum_degree=48
    call solid_angle_reference(64,reference)
    call solid_angle_reference(96,refined)
    call record(abs(reference-refined)<1.e-4_dp*abs(refined), &
        'independent near-normal solid-angle reference is refined')
    call integrate_parametric_bem_p2_pair(plane,1,2,options,v,k,w,c,report,status)
    write(*,'(a,2i8,3es16.7)') 'BEM3D normal pair status degree K reference ratio ', &
        status,report%maximum_degree_used,sum(k),refined,report%error_ratio
    call record(status==0,'resolved near-normal peak reaches requested tolerance')
    if (status==0) call record(abs(sum(k)-refined)<3.e-3_dp*abs(refined), &
        'raw double layer agrees with independent solid-angle integral')
    plane%vertices(3,:,1)=1.e-12_dp
    options%relative_tolerance=0.03_dp
    options%absolute_tolerance=1.e-7_dp
    options%maximum_degree=24
    call integrate_parametric_bem_p2_pair(plane,1,2,options,v,k,w,c,report,status)
    ! h*r/(h*h+r*r)**(3/2) integrates to one on [0,infinity).
    ! The angular factor gives 1/2 per interior target; target area is 1/2.
    ! At h=1e-12 the O(h*log(h)) boundary correction is negligible here.
    write(*,'(a,2i8,2es16.7)') 'BEM3D unresolved normal status degree K ratio ', &
        status,report%maximum_degree_used,sum(k),report%error_ratio
    call record(status==3 .or. &
        (status==0 .and. abs(sum(k)-0.25_dp)<3.e-3_dp), &
        'an unresolved normal peak cannot falsely claim quadrature accuracy')
    call check_summary('Independent near and touching pair integrals')
    if (.not. passed) error stop 1
contains
    subroutine solid_angle_reference(degree,value)
        integer, intent(in) :: degree
        real(dp), intent(out) :: value
        real(dp), allocatable :: xi(:),eta(:),weights(:)
        real(dp) :: x(3),a(3),b(3),r1(3),r2(3),r3(3),cross(3)
        real(dp) :: numerator,denominator,n1,n2,n3
        integer :: t,s
        call triangle_duffy_quadrature(degree,xi,eta,weights,s)
        if (s/=0) error stop 'solid-angle reference quadrature'
        value=0
        do t=1,size(weights)
            call plane%evaluate(1,[xi(t),eta(t)],x,a,b,s)
            r1=plane%vertices(:,1,2)-x
            r2=plane%vertices(:,2,2)-x
            r3=plane%vertices(:,3,2)-x
            cross=[r2(2)*r3(3)-r2(3)*r3(2),r2(3)*r3(1)-r2(1)*r3(3), &
                r2(1)*r3(2)-r2(2)*r3(1)]
            n1=norm2(r1); n2=norm2(r2); n3=norm2(r3)
            numerator=dot_product(r1,cross)
            denominator=n1*n2*n3+dot_product(r1,r2)*n3+ &
                dot_product(r2,r3)*n1+dot_product(r3,r1)*n2
            value=value-weights(t)*atan2(numerator,denominator)/ &
                (2*acos(-1.0_dp))
        end do
    end subroutine

    subroutine edge_reference(degree,value)
        integer, intent(in) :: degree
        real(dp), intent(out) :: value
        real(dp), allocatable :: xi(:),eta(:),weights(:)
        real(dp) :: x(3),a(3),b(3),first(2),second(2),length,height,potential
        real(dp) :: numerator,denominator
        integer :: t,e,next,s
        call triangle_duffy_quadrature(degree,xi,eta,weights,s)
        if (s/=0) error stop 'reference quadrature'
        value=0
        do t=1,size(weights)
            call plane%evaluate(1,[xi(t),eta(t)],x,a,b,s)
            potential=0
            do e=1,3
                next=mod(e,3)+1
                first=plane%vertices(1:2,e,2)-x(1:2)
                second=plane%vertices(1:2,next,2)-x(1:2)
                length=norm2(second-first)
                height=(first(1)*second(2)-first(2)*second(1))/length
                numerator=norm2(first)+norm2(second)+length
                denominator=norm2(first)+norm2(second)-length
                potential=potential+height*log(numerator/denominator)
            end do
            value=value+weights(t)*potential/(4*acos(-1.0_dp))
        end do
    end subroutine
    subroutine record(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        passed=passed.and.ok
        call check_condition(ok,label)
    end subroutine
end program
