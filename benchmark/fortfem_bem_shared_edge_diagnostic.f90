module shared_edge_diagnostic_geometry
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

program fortfem_bem_shared_edge_diagnostic
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use shared_edge_diagnostic_geometry, only: rectangle_pair_t
    use fortfem_parametric_bem_high_order, only: bem_pair_options_t, &
        bem_pair_report_t, bem_pair_degree_trace_t, integrate_parametric_bem_p2_pair
    implicit none
    type(rectangle_pair_t) :: pair
    type(bem_pair_options_t) :: options
    type(bem_pair_report_t) :: report
    type(bem_pair_degree_trace_t) :: trace
    real(dp) :: v(3,3),k(3,6),w(6,6),correction(3,6)
    real(dp) :: reference_v(3,3),reference_w(6,6),reference_sum,aspect
    real(dp) :: errors(3),ratios(4),kc_error
    character(len=4096) :: reference_path,trace_path
    character(len=64) :: schema
    integer :: status,unit,i,j,l,last,ios
    logical :: qualified

    call get_environment_variable('BEM_SHARED_EDGE_REFERENCE',reference_path, &
        status=ios)
    if (ios/=0) error stop 'missing public reference path'
    call get_environment_variable('BEM_SHARED_EDGE_TRACE',trace_path,status=ios)
    if (ios/=0) error stop 'missing trace output path'
    open(newunit=unit,file=trim(reference_path),status='old',action='read')
    read(unit,'(a)') schema
    if (trim(schema)/='B1_SHARED_EDGE_REFERENCE_V1') error stop 'reference schema'
    read(unit,*) aspect
    if (aspect/=40.0_dp) error stop 'reference aspect'
    read(unit,*) reference_v
    read(unit,*) reference_w
    read(unit,*) reference_sum
    close(unit)
    if (.not. all(ieee_is_finite(reference_v))) error stop 'reference V finite'
    if (.not. all(ieee_is_finite(reference_w))) error stop 'reference W finite'
    if (.not. ieee_is_finite(reference_sum)) error stop 'reference sum finite'
    if (minval(reference_v)<=0 .or. reference_sum<=0) error stop 'reference scale'
    pair%vertices(:,:,1)=reshape([0.0_dp,0.0_dp,0.0_dp, &
        aspect,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,3])
    pair%vertices(:,:,2)=reshape([aspect,0.0_dp,0.0_dp, &
        aspect,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp],[3,3])
    options%relative_tolerance=1.e-8_dp
    options%maximum_degree=64
    ! Other options, including atol1e-9 and work cap100M, remain native defaults.
    call integrate_parametric_bem_p2_pair(pair,1,2,options,v,k,w,correction, &
        report,status,trace)
    last=trace%count
    if (last<1 .or. last>size(trace%degree)) error stop 'trace count'
    open(newunit=unit,file=trim(trace_path),status='new',action='write')
    write(unit,'(a)') 'degree block row col value; metadata includes ratios V K W C'
    do l=1,last
        ratios=trace%ratios(:,l)
        write(unit,*) 'LEVEL',trace%degree(l),trace%native_status(l), &
            trace%resolved(l),trace%evaluations(l),ratios
        do j=1,3
            do i=1,3
                write(unit,*) trace%degree(l),' V ',i,j,trace%v(i,j,l)
            end do
        end do
        do j=1,6
            do i=1,3
                write(unit,*) trace%degree(l),' K ',i,j,trace%k(i,j,l)
                write(unit,*) trace%degree(l),' C ',i,j,trace%correction(i,j,l)
            end do
            do i=1,6
                write(unit,*) trace%degree(l),' W ',i,j,trace%w(i,j,l)
            end do
        end do
    end do
    close(unit)
    errors(1)=maxval(abs(trace%v(:,:,last)-reference_v))/maxval(abs(reference_v))
    errors(2)=maxval(abs(trace%w(:,:,last)-reference_w))/maxval(abs(reference_w))
    errors(3)=abs(sum(trace%v(:,:,last))/reference_sum-1)
    kc_error=max(maxval(abs(trace%k(:,:,last))), &
        maxval(abs(trace%correction(:,:,last))))
    qualified=status==0 .and. all(errors<1.e-7_dp) .and. kc_error<1.e-13_dp
    if (.not. all(ieee_is_finite(errors))) qualified=.false.
    write(*,*) 'NATIVE_STATUS DEGREE COUNT QUALIFIED',status, &
        report%maximum_degree_used,last,qualified
    write(*,*) 'REFERENCE_ERRORS V W SUMV KC',errors,kc_error
    write(*,*) 'LAST_BLOCK_RATIOS V K W C',trace%ratios(:,last)
    if (status/=0) error stop 3
    if (.not. qualified) error stop 'independent shared-edge reference failed'
end program
