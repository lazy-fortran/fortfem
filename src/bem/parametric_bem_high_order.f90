module fortfem_parametric_bem_high_order
    !! Continuous P2 potentials / discontinuous P1 fluxes on curved triangles.
    !! Successive pair rules estimate quadrature error; they are not certified
    !! bounds. Status 3 means the degree or per-pair kernel-work cap was reached.
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortfem_kinds, only: dp
    use fortfem_parametric_laplace_bem_3d, only: bem_surface_t, bem_operators_t
    use fortfem_equilibrated_flux_topology, only: topology_t, build_topology
    use fortfem_generated_p1_basis, only: generated_p1_jet
    use fortfem_generated_p2_basis, only: generated_p2_jet
    use fortfem_generated_bem_frame, only: bem_frame
    use fortfem_generated_bem_kernel, only: bem_kernel
    use fortfem_generated_bem_curl_density, only: bem_curl_density
    use fortfem_generated_bem_trace_product, only: bem_trace_product
    use fortfem_generated_bem_hyper_product, only: bem_hyper_product
    use fortfem_generated_bem_fan_map, only: bem_fan_map
    use fortfem_triangle_duffy_quadrature, only: triangle_duffy_quadrature
    use fortfem_bem_circulant, only: bem_circulant_t, circulant_matrix_t
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none
    private
    public :: bem_pair_options_t, bem_pair_report_t, build_p2_trace_connectivity
    public :: assemble_parametric_bem_p2, integrate_parametric_bem_p2_pair
    public :: assemble_parametric_bem_p2_circulant
    public :: assignment(=)
    type :: bem_pair_options_t
        real(dp) :: relative_tolerance = 1.e-3_dp
        real(dp) :: absolute_tolerance = 1.e-9_dp
        integer :: initial_degree = 8, maximum_degree = 48, degree_increment = 4
        integer(int64) :: maximum_evaluations = 100000000_int64
    end type
    type :: bem_pair_report_t
        integer :: maximum_degree_used = 0, accepted_pairs = 0, failed_pairs = 0
        integer(int64) :: kernel_evaluations = 0_int64
        real(dp) :: error_ratio = 0.0_dp, raw_v_symmetry = 0.0_dp
        real(dp) :: raw_w_symmetry = 0.0_dp
    end type
    type :: point_t
        real(dp) :: x(3), normal(3), jacobian, phi(6), psi(3), curl(3,6)
    end type
    type :: block_t
        real(dp) :: v(3,3), k(3,6), w(6,6), c(3,6)
    end type
    type :: rule_bank_t
        real(dp), allocatable :: uv(:,:,:), weights(:,:), nodes(:,:), line_w(:,:)
        integer, allocatable :: count(:), line_count(:), degree(:)
    end type
    interface assignment(=)
        module procedure assign_bank
    end interface
    real(dp), parameter :: scale = 1.0_dp/(4*acos(-1.0_dp))
    real(dp), parameter :: corners(2,3) = reshape( &
        [0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,3])
contains
    subroutine assign_bank(lhs,rhs)
        type(rule_bank_t), intent(out) :: lhs
        type(rule_bank_t), intent(in) :: rhs
        if (allocated(rhs%uv)) lhs%uv=rhs%uv
        if (allocated(rhs%weights)) lhs%weights=rhs%weights
        if (allocated(rhs%nodes)) lhs%nodes=rhs%nodes
        if (allocated(rhs%line_w)) lhs%line_w=rhs%line_w
        if (allocated(rhs%count)) lhs%count=rhs%count
        if (allocated(rhs%line_count)) lhs%line_count=rhs%line_count
        if (allocated(rhs%degree)) lhs%degree=rhs%degree
    end subroutine

    subroutine build_p2_trace_connectivity(faces,trace_dofs,flux_dofs,status)
        integer, intent(in) :: faces(:,:)
        integer, allocatable, intent(out) :: trace_dofs(:,:),flux_dofs(:,:)
        integer, intent(out) :: status
        type(topology_t) :: topo
        integer, allocatable :: counts(:),signs(:),labels(:)
        integer :: nv,nf,ne,i,a,e,lo,hi
        status=1
        if (size(faces,1)/=3 .or. size(faces,2)<4) return
        if (any(faces<1)) return
        nf=size(faces,2); nv=maxval(faces)
        do i=1,nf
            if (faces(1,i)==faces(2,i) .or. faces(2,i)==faces(3,i) .or. &
                faces(3,i)==faces(1,i)) return
        end do
        call build_topology(nv,faces,topo,status)
        if (status/=0) return
        status=1
        ne=maxval(topo%cell_edge)
        allocate(counts(ne),signs(ne),labels(nv))
        counts=0; signs=0
        labels=[(i,i=1,nv)]
        do i=1,nf
            do a=1,3
                e=topo%cell_edge(a,i)
                counts(e)=counts(e)+1
                signs(e)=signs(e)+topo%cell_orientation(a,i)
                lo=min(labels(faces(a,i)),labels(faces(mod(a,3)+1,i)))
                hi=max(labels(faces(a,i)),labels(faces(mod(a,3)+1,i)))
                where(labels==hi) labels=lo
            end do
        end do
        if (any(counts/=2) .or. any(signs/=0)) return
        if (any(labels/=labels(1))) return
        if (mod(2-nv+ne-nf,2)/=0 .or. 2-nv+ne-nf<0) return
        allocate(trace_dofs(6,nf),flux_dofs(3,nf))
        trace_dofs(1:3,:)=faces
        trace_dofs(4:6,:)=nv+topo%cell_edge
        do i=1,nf
            flux_dofs(:,i)=3*(i-1)+[1,2,3]
        end do
        status=0
    end subroutine

    subroutine make_bank(options,bank,status)
        type(bem_pair_options_t), intent(in) :: options
        type(rule_bank_t), intent(out) :: bank
        integer, intent(out) :: status
        real(dp), allocatable :: xi(:),eta(:),weights(:)
        integer, allocatable :: degrees(:)
        integer :: levels,maxq,maxn,l,n,d,used,previous_nodes
        status=1
        if (.not. ieee_is_finite(options%relative_tolerance)) return
        if (.not. ieee_is_finite(options%absolute_tolerance)) return
        if (options%relative_tolerance<0 .or. options%absolute_tolerance<=0) return
        if (options%initial_degree<2 .or. options%maximum_degree>128) return
        if (options%maximum_degree<options%initial_degree) return
        if (options%degree_increment<1 .or. options%degree_increment>128) return
        if (options%maximum_evaluations<1) return
        levels=1+(options%maximum_degree-options%initial_degree+ &
            options%degree_increment-1)/options%degree_increment
        ! Polynomial-degree requests can select identical Gauss rules.
        ! Count actual rule changes, including a capped final degree.
        allocate(degrees(levels))
        used=0; previous_nodes=0
        do l=1,levels
            d=min(options%maximum_degree, &
                options%initial_degree+(l-1)*options%degree_increment)
            n=(d+3)/2
            if (n==previous_nodes) cycle
            used=used+1
            degrees(used)=d
            previous_nodes=n
        end do
        levels=used
        maxn=(options%maximum_degree+3)/2
        maxq=maxn**2
        allocate(bank%uv(2,maxq,levels),bank%weights(maxq,levels), &
            bank%nodes(maxn,levels),bank%line_w(maxn,levels), &
            bank%count(levels),bank%line_count(levels),bank%degree(levels))
        do l=1,levels
            d=degrees(l)
            call triangle_duffy_quadrature(d,xi,eta,weights,status)
            if (status/=0) return
            n=size(weights)
            bank%degree(l)=d; bank%count(l)=n
            bank%uv(1,:n,l)=xi; bank%uv(2,:n,l)=eta
            bank%weights(:n,l)=weights
            n=(d+3)/2
            bank%line_count(l)=n
            call gauss_legendre_ab(n,0.0_dp,1.0_dp, &
                bank%nodes(:n,l),bank%line_w(:n,l))
        end do
        status=0
    end subroutine

    pure subroutine sample(surface,panel,uv,point,status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: uv(2)
        type(point_t), intent(out) :: point
        integer, intent(out) :: status
        real(dp) :: a(3),b(3),g(2),h(2,2),na,nb
        integer :: i
        call surface%evaluate(panel,uv,point%x,a,b,status)
        if (status/=0) return
        status=2
        if (.not. all(ieee_is_finite(point%x))) return
        if (.not. all(ieee_is_finite(a))) return
        if (.not. all(ieee_is_finite(b))) return
        na=norm2(a); nb=norm2(b)
        if (min(na,nb)<=tiny(1.0_dp)) return
        if (abs(dot_product(a,b))>=(1-64*epsilon(1.0_dp))*na*nb) return
        call bem_frame(a(1),a(2),a(3),b(1),b(2),b(3), &
            point%jacobian,point%normal(1),point%normal(2),point%normal(3))
        do i=1,3
            call generated_p1_jet(i,uv(1),uv(2),point%psi(i),g,h)
        end do
        do i=1,6
            call generated_p2_jet(i,uv(1),uv(2),point%phi(i),g,h)
            call bem_curl_density(a(1),a(2),a(3),b(1),b(2),b(3),g(1),g(2), &
                point%curl(1,i),point%curl(2,i),point%curl(3,i))
        end do
        status=0
    end subroutine

    pure subroutine closest(surface,panel,x,uv,status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: panel
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: uv(2)
        integer, intent(out) :: status
        real(dp) :: p(3),a(3),b(3),d(3),aa,ab,bb,det,step(2),t
        integer :: iteration
        uv=[1.0_dp/3,1.0_dp/3]
        do iteration=1,12
            call surface%evaluate(panel,uv,p,a,b,status)
            if (status/=0) return
            status=2
            if (.not. all(ieee_is_finite(p))) return
            if (.not. all(ieee_is_finite(a))) return
            if (.not. all(ieee_is_finite(b))) return
            d=p-x; aa=dot_product(a,a); ab=dot_product(a,b); bb=dot_product(b,b)
            det=aa*bb-ab*ab
            if (det<=64*epsilon(1.0_dp)*aa*bb) return
            step=[(bb*dot_product(a,d)-ab*dot_product(b,d))/det, &
                (aa*dot_product(b,d)-ab*dot_product(a,d))/det]
            uv=max(uv-step,0.0_dp)
            if (sum(uv)>1) then
                t=max(0.0_dp,min(1.0_dp,(uv(1)-uv(2)+1)/2))
                uv=[t,1-t]
            end if
        end do
        status=0
    end subroutine

    pure subroutine clear(block)
        type(block_t), intent(out) :: block
        block%v=0; block%k=0; block%w=0; block%c=0
    end subroutine

    pure subroutine add_integrand(target,source,factor,block,status)
        type(point_t), intent(in) :: target,source
        real(dp), intent(in) :: factor
        type(block_t), intent(inout) :: block
        integer, intent(out) :: status
        real(dp) :: d(3),green,double,value
        integer :: a,b
        status=2
        d=target%x-source%x
        if (norm2(d)<=tiny(1.0_dp)) return
        call bem_kernel(d(1),d(2),d(3),source%normal(1),source%normal(2), &
            source%normal(3),green,double)
        do b=1,3
            do a=1,3
                call bem_trace_product(scale*green,target%jacobian, &
                    source%jacobian,target%psi(a),source%psi(b),value)
                block%v(a,b)=block%v(a,b)+factor*value
            end do
        end do
        do b=1,6
            do a=1,3
                call bem_trace_product(scale*double,target%jacobian, &
                    source%jacobian,target%psi(a),source%phi(b),value)
                block%k(a,b)=block%k(a,b)+factor*value
                call bem_trace_product(scale*double,target%jacobian, &
                    source%jacobian,target%psi(a),target%phi(b),value)
                block%c(a,b)=block%c(a,b)+factor*value
            end do
            do a=1,6
                call bem_hyper_product(scale*green,target%curl(1,a), &
                    target%curl(2,a),target%curl(3,a),source%curl(1,b), &
                    source%curl(2,b),source%curl(3,b),value)
                block%w(a,b)=block%w(a,b)+factor*value
            end do
        end do
        status=0
    end subroutine

    pure subroutine integrate_rule(surface,target_panel,source_panel,bank,l, &
            target_points,block,evaluations,resolved,status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: target_panel,source_panel,l
        type(rule_bank_t), intent(in) :: bank
        type(point_t), intent(in) :: target_points(:)
        type(block_t), intent(out) :: block
        integer(int64), intent(out) :: evaluations
        logical, intent(out) :: resolved
        integer, intent(out) :: status
        type(point_t) :: source
        real(dp) :: apex(2),e1(2),e2(2),det,uv(2),factor,fan_weight
        real(dp) :: probe_x(3),probe_a(3),probe_b(3),distance,length_scale
        integer :: t,edge,i,j,n
        call clear(block)
        evaluations=0
        resolved=.true.
        n=bank%line_count(l)
        do t=1,bank%count(l)
            apex=bank%uv(:,t,l)
            if (target_panel/=source_panel) then
                call closest(surface,source_panel,target_points(t)%x,apex,status)
                if (status/=0) return
                if (minval(apex)>64*epsilon(1.0_dp) .and. &
                    sum(apex)<1-64*epsilon(1.0_dp)) then
                    ! An interior normal peak has physical width distance.
                    ! The squared-radial fan must sample within that scale;
                    ! otherwise two rules can agree while both miss the peak.
                    call surface%evaluate(source_panel,apex,probe_x, &
                        probe_a,probe_b,status)
                    if (status/=0) return
                    status=2
                    if (.not. all(ieee_is_finite(probe_x))) return
                    if (.not. all(ieee_is_finite(probe_a))) return
                    if (.not. all(ieee_is_finite(probe_b))) return
                    length_scale=max(norm2(probe_a),norm2(probe_b))
                    distance=norm2(target_points(t)%x-probe_x)
                    if (.not. ieee_is_finite(length_scale)) return
                    if (.not. ieee_is_finite(distance)) return
                    if (length_scale<=tiny(1.0_dp)) return
                    if (distance>0 .and. &
                        distance<2*length_scale*bank%nodes(1,l)**2) then
                        resolved=.false.
                        status=0
                        return
                    end if
                    status=0
                end if
            end if
            do edge=1,3
                e1=corners(:,edge)-apex
                e2=corners(:,mod(edge,3)+1)-apex
                det=abs(e1(1)*e2(2)-e1(2)*e2(1))
                if (det<64*epsilon(1.0_dp)) cycle
                do i=1,n
                    do j=1,n
                        call bem_fan_map(apex(1),apex(2),e1(1),e1(2), &
                            e2(1),e2(2),bank%nodes(i,l),bank%nodes(j,l),det, &
                            uv(1),uv(2),fan_weight)
                        call sample(surface,source_panel,uv,source,status)
                        if (status/=0) return
                        factor=bank%weights(t,l)*bank%line_w(i,l)* &
                            bank%line_w(j,l)*fan_weight
                        call add_integrand(target_points(t),source,factor, &
                            block,status)
                        evaluations=evaluations+1
                        if (status/=0) return
                    end do
                end do
            end do
        end do
        status=0
    end subroutine

    pure real(dp) function difference_ratio(current,previous,options) result(r)
        type(block_t), intent(in) :: current,previous
        type(bem_pair_options_t), intent(in) :: options
        real(dp) :: a,b
        a=options%absolute_tolerance; b=options%relative_tolerance
        r=max(maxval(abs(current%v-previous%v))/(a+b*maxval(abs(current%v))), &
            maxval(abs(current%k-previous%k))/(a+b*maxval(abs(current%k))), &
            maxval(abs(current%w-previous%w))/(a+b*maxval(abs(current%w))), &
            maxval(abs(current%c-previous%c))/(a+b*maxval(abs(current%c))))
    end function

    pure subroutine converge_pair(surface,i,j,options,bank,points,block, &
            report,status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i,j
        type(bem_pair_options_t), intent(in) :: options
        type(rule_bank_t), intent(in) :: bank
        type(point_t), intent(in) :: points(:,:)
        type(block_t), intent(out) :: block
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(block_t) :: previous
        integer(int64) :: bound,evaluations
        integer :: l
        logical :: resolved,have_previous
        call clear(block)
        report=bem_pair_report_t()
        report%error_ratio=huge(1.0_dp)
        have_previous=.false.
        do l=1,size(bank%degree)
            bound=3_int64*bank%count(l)*int(bank%line_count(l),int64)**2
            if (bound>options%maximum_evaluations-report%kernel_evaluations) exit
            call integrate_rule(surface,i,j,bank,l,points(:,l),block, &
                evaluations,resolved,status)
            report%kernel_evaluations=report%kernel_evaluations+evaluations
            report%maximum_degree_used=bank%degree(l)
            if (status/=0) then
                report%failed_pairs=1
                return
            end if
            if (.not. resolved) then
                have_previous=.false.
                cycle
            end if
            if (have_previous) then
                report%error_ratio=difference_ratio(block,previous,options)
                if (.not. ieee_is_finite(report%error_ratio)) then
                    status=2
                    report%failed_pairs=1
                    return
                end if
                if (report%error_ratio<=1) then
                    status=0
                    report%accepted_pairs=1
                    return
                end if
            end if
            previous=block
            have_previous=.true.
        end do
        status=3
        report%failed_pairs=1
    end subroutine

    subroutine integrate_parametric_bem_p2_pair(surface,i,j,options, &
            v,k,w,correction,report,status)
        !! Raw local forms; no symmetry averaging and no closed-surface jump.
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: i,j
        type(bem_pair_options_t), intent(in) :: options
        real(dp), intent(out) :: v(3,3),k(3,6),w(6,6),correction(3,6)
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(rule_bank_t) :: bank
        type(point_t), allocatable :: points(:,:)
        type(block_t) :: block
        integer :: t,l
        v=0; k=0; w=0; correction=0
        report=bem_pair_report_t()
        status=1
        if (min(i,j)<1) return
        call make_bank(options,bank,status)
        if (status/=0) return
        allocate(points(size(bank%weights,1),size(bank%degree)))
        do l=1,size(bank%degree)
            do t=1,bank%count(l)
                call sample(surface,i,bank%uv(:,t,l),points(t,l),status)
                if (status/=0) return
            end do
        end do
        call converge_pair(surface,i,j,options,bank,points,block,report,status)
        if (status/=0) return
        v=block%v; k=block%k; w=block%w; correction=block%c
    end subroutine

    subroutine assemble_parametric_bem_p2(surface,faces,options,op, &
            trace_dofs,flux_dofs,report,status)
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: faces(:,:)
        type(bem_pair_options_t), intent(in) :: options
        type(bem_operators_t), intent(out) :: op
        integer, allocatable, intent(out) :: trace_dofs(:,:),flux_dofs(:,:)
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(rule_bank_t) :: bank
        type(point_t), allocatable :: points(:,:,:)
        type(block_t) :: block
        type(bem_pair_report_t) :: pair_report
        integer :: nf,nt,nv,ne,i,j,l,t,a,b,local_status,errors,dmax,geometry_errors
        integer(int64) :: evaluations
        real(dp) :: value,ratio
        report=bem_pair_report_t()
        call build_p2_trace_connectivity(faces,trace_dofs,flux_dofs,status)
        if (status/=0) return
        call make_bank(options,bank,status)
        if (status/=0) return
        nf=size(faces,2); nt=maxval(trace_dofs); nv=maxval(faces); ne=nt-nv
        op%genus=(2-nv+ne-nf)/2
        allocate(points(size(bank%weights,1),nf,size(bank%degree)))
        errors=0
        !$omp parallel do default(shared) private(i,l,t,local_status) &
        !$omp reduction(+:errors) schedule(static)
        do i=1,nf
            do l=1,size(bank%degree)
                do t=1,bank%count(l)
                    call sample(surface,i,bank%uv(:,t,l),points(t,i,l),local_status)
                    if (local_status/=0) errors=errors+1
                end do
            end do
        end do
        !$omp end parallel do
        status=2
        if (errors/=0) return
        allocate(op%v(3*nf,3*nf),op%k(3*nf,nt),op%w(nt,nt), &
            op%mass(3*nf,nt),op%area(nf))
        op%v=0; op%k=0; op%w=0; op%mass=0; op%area=0
        ! Mass is independent of the pair tolerance. Use the finest bank rule.
        l=size(bank%degree)
        do i=1,nf
            do t=1,bank%count(l)
                op%area(i)=op%area(i)+bank%weights(t,l)*points(t,i,l)%jacobian
                do b=1,6
                    do a=1,3
                        call bem_trace_product(bank%weights(t,l), &
                            points(t,i,l)%jacobian,1.0_dp,points(t,i,l)%psi(a), &
                            points(t,i,l)%phi(b),value)
                        op%mass(flux_dofs(a,i),trace_dofs(b,i))= &
                            op%mass(flux_dofs(a,i),trace_dofs(b,i))+value
                    end do
                end do
            end do
        end do
        evaluations=0; dmax=0; ratio=0; geometry_errors=0
        !$omp parallel do default(shared) schedule(dynamic,1) &
        !$omp private(i,j,a,b,block,pair_report,local_status) &
        !$omp reduction(+:errors,evaluations,geometry_errors) reduction(max:dmax,ratio)
        do i=1,nf
            do j=1,nf
                call converge_pair(surface,i,j,options,bank,points(:,i,:), &
                    block,pair_report,local_status)
                evaluations=evaluations+pair_report%kernel_evaluations
                dmax=max(dmax,pair_report%maximum_degree_used)
                ratio=max(ratio,pair_report%error_ratio)
                if (local_status/=0) then
                    errors=errors+1
                    if (local_status==2) geometry_errors=geometry_errors+1
                    cycle
                end if
                op%v(flux_dofs(:,i),flux_dofs(:,j))=block%v
                do b=1,6
                    op%k(flux_dofs(:,i),trace_dofs(b,j))= &
                        op%k(flux_dofs(:,i),trace_dofs(b,j))+block%k(:,b)
                    op%k(flux_dofs(:,i),trace_dofs(b,i))= &
                        op%k(flux_dofs(:,i),trace_dofs(b,i))-block%c(:,b)
                    do a=1,6
                        !$omp atomic update
                        op%w(trace_dofs(a,i),trace_dofs(b,j))= &
                            op%w(trace_dofs(a,i),trace_dofs(b,j))+block%w(a,b)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        report%failed_pairs=errors
        report%accepted_pairs=nf**2-errors
        report%maximum_degree_used=dmax
        report%kernel_evaluations=evaluations
        report%error_ratio=ratio
        status=3
        if (geometry_errors/=0) status=2
        if (errors/=0) return
        if (.not. all(ieee_is_finite(op%v))) then
            status=2
            return
        end if
        if (.not. all(ieee_is_finite(op%k))) then
            status=2
            return
        end if
        if (.not. all(ieee_is_finite(op%w))) then
            status=2
            return
        end if
        report%raw_v_symmetry=maxval(abs(op%v-transpose(op%v)))/ &
            max(tiny(1.0_dp),maxval(abs(op%v)))
        report%raw_w_symmetry=maxval(abs(op%w-transpose(op%w)))/ &
            max(tiny(1.0_dp),maxval(abs(op%w)))
        op%k=op%k-0.5_dp*op%mass
        op%v=0.5_dp*(op%v+transpose(op%v))
        op%w=0.5_dp*(op%w+transpose(op%w))
        op%adjoint=transpose(op%k)
        status=0
    end subroutine

    subroutine sector_trace_map(td, nphi, map, local, sector, status)
        !! Infer the trace rotation from corresponding local panel DOFs.
        !! A free cyclic action is required (no fixed-axis trace DOFs).
        integer, intent(in) :: td(:, :), nphi
        integer, allocatable, intent(out) :: map(:, :), local(:), sector(:)
        integer, intent(out) :: status
        integer, allocatable :: rotation(:)
        integer :: nt, np, nf, i, j, a, q, r, count
        status = 1
        nf = size(td, 2); nt = maxval(td)
        if (nphi < 1) return
        if (mod(nf, nphi) /= 0 .or. mod(nt, nphi) /= 0) return
        np = nf/nphi
        allocate (rotation(nt), map(nt/nphi, nphi), local(nt), sector(nt))
        rotation = 0; local = 0; sector = 0
        do i = 1, nf
            j = modulo(i - 1 + np, nf) + 1
            do a = 1, 6
                q = td(a, i); r = td(a, j)
                if (rotation(q) /= 0) then
                    if (rotation(q) /= r) return
                end if
                rotation(q) = r
            end do
        end do
        if (any(rotation < 1)) return
        count = 0
        do i = 1, nt
            if (local(i) /= 0) cycle
            count = count + 1
            if (count > nt/nphi) return
            q = i
            do r = 1, nphi
                if (local(q) /= 0) return
                map(count, r) = q; local(q) = count; sector(q) = r
                q = rotation(q)
            end do
            if (q /= i) return
        end do
        if (count /= nt/nphi) return
        status = 0
    end subroutine

    subroutine check_sector_geometry(surface, np, nphi, status)
        !! Validate the declared z rotation at vertices, edge midpoints and
        !! centroid, including tangents. Exact symmetry between these samples
        !! remains the caller's geometry contract.
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: np, nphi
        integer, intent(out) :: status
        real(dp), parameter :: uv(2, 7) = reshape( &
                               [0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 1._dp, &
                    .5_dp, 0._dp, .5_dp, .5_dp, 0._dp, .5_dp, 1._dp/3, 1._dp/3], [2, 7])
        real(dp) :: reference(3, 3), actual(3, 3), rotated(3, 3), angle, c, s
        integer :: i, j, t, k, local_status
        status = 1
        do j = 1, nphi - 1
            angle = 2*acos(-1._dp)*real(j, dp)/real(nphi, dp)
            c = cos(angle); s = sin(angle)
            do i = 1, np
                do t = 1, 7
                    call surface%evaluate(i, uv(:, t), reference(:, 1), &
                                         reference(:, 2), reference(:, 3), local_status)
                    if (local_status /= 0) then
                        status = 2
                        return
                    end if
                    call surface%evaluate(i + j*np, uv(:, t), actual(:, 1), &
                                          actual(:, 2), actual(:, 3), local_status)
                    if (local_status /= 0) then
                        status = 2
                        return
                    end if
                    if (.not. all(ieee_is_finite(reference))) then
                        status = 2
                        return
                    end if
                    if (.not. all(ieee_is_finite(actual))) then
                        status = 2
                        return
                    end if
                    rotated = reference
                    rotated(1, :) = c*reference(1, :) - s*reference(2, :)
                    rotated(2, :) = s*reference(1, :) + c*reference(2, :)
                    do k = 1, 3
                        if (norm2(actual(:, k) - rotated(:, k)) > &
                            1.e-10_dp*max(1._dp, norm2(reference(:, k)))) return
                    end do
                end do
            end do
        end do
        status = 0
    end subroutine

    subroutine initialize_sector_matrix(matrix, rows, columns)
        type(circulant_matrix_t), intent(out) :: matrix
        integer, intent(in) :: rows(:, :), columns(:, :)
        matrix%row_dofs = rows; matrix%column_dofs = columns
        allocate (matrix%blocks(size(rows, 1), size(columns, 1), size(rows, 2)))
        matrix%blocks = 0._dp
    end subroutine

    subroutine assemble_parametric_bem_p2_circulant(surface, faces, nphi, options, &
                                              op, trace_dofs, flux_dofs, report, status)
        !! Consecutive identical sectors, with corresponding local panel
        !! vertices in the same order. Only target-sector-zero pairs are
        !! integrated. Reports count actual evaluated pairs, not their copies.
        class(bem_surface_t), intent(in) :: surface
        integer, intent(in) :: faces(:, :), nphi
        type(bem_pair_options_t), intent(in) :: options
        type(bem_circulant_t), intent(out) :: op
        integer, allocatable, intent(out) :: trace_dofs(:, :), flux_dofs(:, :)
        type(bem_pair_report_t), intent(out) :: report
        integer, intent(out) :: status
        type(rule_bank_t) :: bank
        type(point_t), allocatable :: points(:, :, :)
        type(block_t) :: block
        type(bem_pair_report_t) :: pair_report
        integer, allocatable :: trace_map(:, :), local(:), sector(:), flux_map(:, :)
        integer :: nf, np, nt, nv, ne, i, j, l, t, a, b, ia, ib, sa, sb, js, jl
        integer :: local_status, errors, dmax, geometry_errors
        integer(int64) :: evaluations
        real(dp) :: value, ratio
        report = bem_pair_report_t()
        status = 1
        if (nphi < 1) return
        if (mod(size(faces, 2), nphi) /= 0) return
        call build_p2_trace_connectivity(faces, trace_dofs, flux_dofs, status)
        if (status /= 0) return
        call sector_trace_map(trace_dofs, nphi, trace_map, local, sector, status)
        if (status /= 0) return
        nf = size(faces, 2); np = nf/nphi
        call check_sector_geometry(surface, np, nphi, status)
        if (status /= 0) return
        call make_bank(options, bank, status)
        if (status /= 0) return
        nt = maxval(trace_dofs); nv = maxval(faces); ne = nt - nv
        op%genus = (2 - nv + ne - nf)/2
        allocate (flux_map(3*np, nphi))
        do j = 1, nphi
            do i = 1, 3*np
                flux_map(i, j) = i + 3*np*(j - 1)
            end do
        end do
        call initialize_sector_matrix(op%v, flux_map, flux_map)
        call initialize_sector_matrix(op%k, flux_map, trace_map)
        call initialize_sector_matrix(op%w, trace_map, trace_map)
        call initialize_sector_matrix(op%mass, flux_map, trace_map)
        allocate (op%area(nf), &
                  points(size(bank%weights, 1), np, size(bank%degree)))
        op%area = 0._dp
        errors = 0
        !$omp parallel do default(shared) private(i,l,t,local_status) &
        !$omp reduction(+:errors) schedule(static)
        do i = 1, np
            do l = 1, size(bank%degree)
                do t = 1, bank%count(l)
                call sample(surface, i, bank%uv(:, t, l), points(t, i, l), local_status)
                    if (local_status /= 0) errors = errors + 1
                end do
            end do
        end do
        !$omp end parallel do
        status = 2
        if (errors /= 0) return
        l = size(bank%degree)
        do i = 1, np
            do t = 1, bank%count(l)
                op%area(i) = op%area(i) + bank%weights(t, l)*points(t, i, l)%jacobian
                do b = 1, 6
                    ib = local(trace_dofs(b, i)); sb = sector(trace_dofs(b, i))
                    do a = 1, 3
                        call bem_trace_product(bank%weights(t, l), &
                              points(t, i, l)%jacobian, 1._dp, points(t, i, l)%psi(a), &
                                               points(t, i, l)%phi(b), value)
                        op%mass%blocks(flux_dofs(a, i), ib, sb) = &
                            op%mass%blocks(flux_dofs(a, i), ib, sb) + value
                    end do
                end do
            end do
        end do
        do j = 1, nphi - 1
            do i = 1, np
                op%area(j*np + i) = op%area(i)
            end do
        end do
        evaluations = 0; dmax = 0; ratio = 0; geometry_errors = 0
        !$omp parallel do collapse(2) default(shared) schedule(dynamic,1) &
        !$omp private(i,j,a,b,ia,ib,sa,sb,js,jl,block,pair_report,local_status) &
        !$omp reduction(+:errors,evaluations,geometry_errors) reduction(max:dmax,ratio)
        do i = 1, np
            do j = 1, nf
                call converge_pair(surface, i, j, options, bank, points(:, i, :), &
                                   block, pair_report, local_status)
                evaluations = evaluations + pair_report%kernel_evaluations
                dmax = max(dmax, pair_report%maximum_degree_used)
                ratio = max(ratio, pair_report%error_ratio)
                if (local_status /= 0) then
                    errors = errors + 1
                    if (local_status == 2) geometry_errors = geometry_errors + 1
                    cycle
                end if
                js = (j - 1)/np + 1; jl = modulo(j - 1, np)
                do b = 1, 3
                    do a = 1, 3
                        op%v%blocks(flux_dofs(a, i), 3*jl + b, js) = block%v(a, b)
                    end do
                end do
                do b = 1, 6
                    ib = local(trace_dofs(b, j)); sb = sector(trace_dofs(b, j))
                    do a = 1, 3
                        !$omp atomic update
                        op%k%blocks(flux_dofs(a, i), ib, sb) = &
                            op%k%blocks(flux_dofs(a, i), ib, sb) + block%k(a, b)
                    end do
                    do a = 1, 6
                        ia = local(trace_dofs(a, i)); sa = sector(trace_dofs(a, i))
                        !$omp atomic update
                        op%w%blocks(ia, ib, modulo(sb - sa, nphi) + 1) = &
                          op%w%blocks(ia, ib, modulo(sb - sa, nphi) + 1) + block%w(a, b)
                    end do
                    ib = local(trace_dofs(b, i)); sb = sector(trace_dofs(b, i))
                    do a = 1, 3
                        !$omp atomic update
                        op%k%blocks(flux_dofs(a, i), ib, sb) = &
                            op%k%blocks(flux_dofs(a, i), ib, sb) - block%c(a, b)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        report%failed_pairs = errors
        report%accepted_pairs = np*nf - errors
        report%maximum_degree_used = dmax
        report%kernel_evaluations = evaluations
        report%error_ratio = ratio
        status = 3
        if (geometry_errors /= 0) status = 2
        if (errors /= 0) return
        status = 2
        if (.not. op%v%valid()) return
        if (.not. op%k%valid()) return
        if (.not. op%w%valid()) return
        op%k%blocks = op%k%blocks - .5_dp*op%mass%blocks
        call op%v%symmetrize(report%raw_v_symmetry)
        call op%w%symmetrize(report%raw_w_symmetry)
        status = 0
    end subroutine
end module fortfem_parametric_bem_high_order
