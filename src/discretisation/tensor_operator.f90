module fortfem_tensor_operator
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortfem_tensor_mesh, only: tensor_mesh_t
    use fortfem_tensor_line_tables, only: tensor_line_tables_t
    use fortfem_gll_line, only: build_gll_line, gll_rule, lagrange_table
    use fortfem_tensor_contractions, only: tensor_contract
    use fortfem_generated_tensor_metric, only: generated_tensor_metric
    use fortfem_element_to_csr, only: csr_matrix_t, csr_build_pattern, &
        csr_find_entry, build_node_to_element, csr_element_colouring
    implicit none
    private
    public :: tensor_operator_t, build_tensor_operator, apply_tensor_operator
    public :: assemble_tensor_operator, tensor_diagonal, tensor_load, tensor_l2_error
    public :: diffusion_coefficient, scalar_function, metric_at
    type :: tensor_operator_t
        type(tensor_mesh_t) :: mesh
        type(tensor_line_tables_t) :: line
        integer :: ni(3)=1, no(3)=1, ncolours=0
        logical :: precompute=.true.
        real(dp), allocatable :: jac(:,:,:,:),coef(:,:,:,:),reaction(:,:)
        real(dp), allocatable :: weight(:),points(:,:,:),volume(:,:)
        real(dp), allocatable :: metric(:,:,:,:),mass(:,:)
        integer, allocatable :: colour_ptr(:),colour_cells(:)
    end type
    abstract interface
        subroutine diffusion_coefficient(x,k,c)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: k(:,:),c
        end subroutine
        function scalar_function(x) result(value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp) :: value
        end function
    end interface
contains
    subroutine metric_at(j,k,w,c,g,m,det)
        real(dp), intent(in) :: j(3,3),k(3,3),w,c
        real(dp), intent(out) :: g(3,3),m,det
        call generated_tensor_metric(j(1,1),j(2,1),j(3,1), &
            j(1,2),j(2,2),j(3,2),j(1,3),j(2,3),j(3,3), &
            k(1,1),k(2,1),k(3,1),k(1,2),k(2,2),k(3,2), &
            k(1,3),k(2,3),k(3,3),w,c,g(1,1),g(2,1),g(3,1), &
            g(1,2),g(2,2),g(3,2),g(1,3),g(2,3),g(3,3),m,det)
    end subroutine

    subroutine build_tensor_operator(mesh,coefficient,op,status,precompute,nq)
        type(tensor_mesh_t), intent(in) :: mesh
        procedure(diffusion_coefficient) :: coefficient
        type(tensor_operator_t), intent(out) :: op
        integer, intent(out) :: status
        logical, intent(in), optional :: precompute
        integer, intent(in), optional :: nq
        type(csr_matrix_t) :: adjacency
        real(dp), allocatable :: nodes(:),w(:),bg(:,:),dg(:,:), &
            a(:,:,:),t1(:),t2(:),values(:),identity(:,:)
        integer, allocatable :: colour(:),fill(:)
        real(dp) :: g(3,3),m,det,kmat(3,3),c
        integer :: ng(3),mq,n,e,comp,dir,i,j,l,point
        status=1
        if (mesh%dim<2 .or. mesh%dim>3) return
        op%mesh=mesh
        if (present(precompute)) op%precompute=precompute
        n=mesh%degree+1
        mq=n
        if (present(nq)) mq=nq
        call build_gll_line(mesh%degree,1,op%line,status,nq=mq, &
            x0=-1.0_dp,x1=1.0_dp)
        if (status/=0) return
        op%ni(1:mesh%dim)=n
        op%no(1:mesh%dim)=mq
        ng=1
        ng(1:mesh%dim)=mesh%geometry_degree+1
        point=product(op%no)
        allocate(op%jac(3,3,point,mesh%nelem), &
            op%coef(3,3,point,mesh%nelem),op%reaction(point,mesh%nelem), &
            op%points(mesh%dim,point,mesh%nelem),op%weight(point), &
            op%volume(point,mesh%nelem))
        if (op%precompute) allocate(op%metric(3,3,point,mesh%nelem), &
            op%mass(point,mesh%nelem))
        call gll_rule(mesh%geometry_degree,nodes,w,status)
        if (status/=0) return
        call lagrange_table(nodes,op%line%xq(:,1),bg,dg)
        allocate(a(mq,ng(1),3),t1(maxval([ng,op%no])**3), &
            t2(maxval([ng,op%no])**3),values(point),identity(1,1))
        identity=1
        do l=1,op%no(3)
            do j=1,mq
                do i=1,mq
                    point=i+mq*(j-1)+mq**2*(l-1)
                    op%weight(point)=op%line%wq(i,1)*op%line%wq(j,1)
                    if (mesh%dim==3) op%weight(point)= &
                        op%weight(point)*op%line%wq(l,1)
                end do
            end do
        end do
        op%jac=0
        if (mesh%dim==2) op%jac(3,3,:,:)=1
        do e=1,mesh%nelem
            do comp=1,mesh%dim
                call contract(bg,bg,bg,mesh%geometry(comp,:,e),values)
                op%points(comp,:,e)=values
                do dir=1,mesh%dim
                    a(:,:,1)=bg; a(:,:,2)=bg; a(:,:,3)=bg
                    a(:,:,dir)=dg
                    call contract(a(:,:,1),a(:,:,2),a(:,:,3), &
                        mesh%geometry(comp,:,e),values)
                    op%jac(comp,dir,:,e)=values
                end do
            end do
            do point=1,product(op%no)
                kmat=0
                kmat(3,3)=1
                call coefficient(op%points(:,point,e), &
                    kmat(1:mesh%dim,1:mesh%dim),c)
                if (maxval(abs(kmat-transpose(kmat)))>1.e-13_dp) return
                op%coef(:,:,point,e)=kmat
                op%reaction(point,e)=c
                call metric_at(op%jac(:,:,point,e),kmat,op%weight(point), &
                    c,g,m,det)
                if (det<=0) then
                    status=2
                    return
                end if
                op%volume(point,e)=op%weight(point)*det
                if (op%precompute) then
                    op%metric(:,:,point,e)=g
                    op%mass(point,e)=m
                end if
            end do
        end do
        ! Colour from node-to-cell adjacency only, without a dense CSR pattern.
        call build_node_to_element(mesh%ndof,mesh%dof, &
            adjacency%elem_ptr,adjacency%elem_idx)
        call csr_element_colouring(adjacency,mesh%dof,colour,op%ncolours)
        allocate(op%colour_ptr(op%ncolours+1), &
            op%colour_cells(mesh%nelem),fill(op%ncolours))
        op%colour_ptr=0
        do e=1,mesh%nelem
            op%colour_ptr(colour(e)+1)=op%colour_ptr(colour(e)+1)+1
        end do
        op%colour_ptr(1)=1
        do i=1,op%ncolours
            op%colour_ptr(i+1)=op%colour_ptr(i+1)+op%colour_ptr(i)
        end do
        fill=op%colour_ptr(1:op%ncolours)
        do e=1,mesh%nelem
            i=colour(e)
            op%colour_cells(fill(i))=e
            fill(i)=fill(i)+1
        end do
        status=0
    contains
        subroutine contract(a1,a2,a3,u,v)
            real(dp), intent(in) :: a1(:,:),a2(:,:),a3(:,:),u(:)
            real(dp), intent(out) :: v(:)
            if (mesh%dim==2) then
                call tensor_contract(ng,op%no,a1,a2,identity,u,v,t1,t2)
            else
                call tensor_contract(ng,op%no,a1,a2,a3,u,v,t1,t2)
            end if
        end subroutine
    end subroutine

    subroutine apply_tensor_operator(op,x,y,eliminate)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        logical, intent(in), optional :: eliminate
        logical :: mask
        integer :: n,m,nl,mq,ws,colour,ce,e,a,dir,q
        mask=.true.
        if (present(eliminate)) mask=eliminate
        n=op%ni(1); m=op%no(1)
        nl=product(op%ni); mq=product(op%no)
        ws=max(n,m)**3
        y=0
        !$omp parallel default(shared) private(colour,ce,e,a,dir,q)
        block
            real(dp), allocatable :: u(:),v(:),tmp(:),grad(:,:),flux(:,:), &
                t1(:),t2(:),b(:,:,:),bt(:,:,:)
            real(dp) :: g(3,3),mass,det,identity(1,1)
            allocate(u(nl),v(nl),tmp(max(nl,mq)),grad(mq,4),flux(mq,4), &
                t1(ws),t2(ws),b(m,n,3),bt(n,m,3))
            identity=1
            do colour=1,op%ncolours
                !$omp do schedule(static)
                do ce=op%colour_ptr(colour),op%colour_ptr(colour+1)-1
                    e=op%colour_cells(ce)
                    do a=1,nl
                        u(a)=x(op%mesh%dof(a,e))
                        if (mask) then
                            if (op%mesh%boundary(op%mesh%dof(a,e))) u(a)=0
                        end if
                    end do
                    do dir=1,op%mesh%dim+1
                        b(:,:,1)=op%line%b(:,:,1)
                        b(:,:,2)=op%line%b(:,:,1)
                        b(:,:,3)=op%line%b(:,:,1)
                        if (dir<=op%mesh%dim) b(:,:,dir)=op%line%d(:,:,1)
                        if (op%mesh%dim==2) then
                            call tensor_contract(op%ni,op%no,b(:,:,1),b(:,:,2), &
                                identity,u,grad(:,dir),t1,t2)
                        else
                            call tensor_contract(op%ni,op%no,b(:,:,1),b(:,:,2), &
                                b(:,:,3),u,grad(:,dir),t1,t2)
                        end if
                    end do
                    flux=0
                    do q=1,mq
                        if (op%precompute) then
                            g=op%metric(:,:,q,e); mass=op%mass(q,e)
                        else
                            call metric_at(op%jac(:,:,q,e),op%coef(:,:,q,e), &
                                op%weight(q),op%reaction(q,e),g,mass,det)
                        end if
                        do dir=1,op%mesh%dim
                            flux(q,dir)=dot_product(g(dir,1:op%mesh%dim), &
                                grad(q,1:op%mesh%dim))
                        end do
                        flux(q,op%mesh%dim+1)=mass*grad(q,op%mesh%dim+1)
                    end do
                    v=0
                    do dir=1,op%mesh%dim+1
                        bt(:,:,1)=transpose(op%line%b(:,:,1))
                        bt(:,:,2)=transpose(op%line%b(:,:,1))
                        bt(:,:,3)=transpose(op%line%b(:,:,1))
                        if (dir<=op%mesh%dim) &
                            bt(:,:,dir)=transpose(op%line%d(:,:,1))
                        if (op%mesh%dim==2) then
                            call tensor_contract(op%no,op%ni,bt(:,:,1),bt(:,:,2), &
                                identity,flux(:,dir),tmp,t1,t2)
                        else
                            call tensor_contract(op%no,op%ni,bt(:,:,1),bt(:,:,2), &
                                bt(:,:,3),flux(:,dir),tmp,t1,t2)
                        end if
                        v=v+tmp(1:nl)
                    end do
                    do a=1,nl
                        y(op%mesh%dof(a,e))=y(op%mesh%dof(a,e))+v(a)
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        if (mask) then
            where(op%mesh%boundary) y=x
        end if
    end subroutine

    subroutine assemble_tensor_operator(op,matrix)
        type(tensor_operator_t), intent(in) :: op
        type(csr_matrix_t), intent(out) :: matrix
        real(dp), allocatable :: basis(:,:),gradient(:,:,:),element(:,:)
        real(dp) :: g(3,3),mass,det,val
        integer :: a,b,q,e,i,j,k,l,dir,ia(3),iq(3),n,m,nl,mq,entry
        n=op%ni(1); m=op%no(1)
        nl=product(op%ni); mq=product(op%no)
        allocate(basis(mq,nl),gradient(3,mq,nl),element(nl,nl))
        gradient=0
        ! Deliberately unfactored quadrature oracle, independent of MF schedules.
        do a=1,nl
            ia=[modulo(a-1,n)+1,modulo((a-1)/n,n)+1,(a-1)/n**2+1]
            do q=1,mq
                iq=[modulo(q-1,m)+1,modulo((q-1)/m,m)+1,(q-1)/m**2+1]
                val=1
                do dir=1,op%mesh%dim
                    val=val*op%line%b(iq(dir),ia(dir),1)
                end do
                basis(q,a)=val
                do dir=1,op%mesh%dim
                    val=1
                    do l=1,op%mesh%dim
                        if (l==dir) then
                            val=val*op%line%d(iq(l),ia(l),1)
                        else
                            val=val*op%line%b(iq(l),ia(l),1)
                        end if
                    end do
                    gradient(dir,q,a)=val
                end do
            end do
        end do
        call csr_build_pattern(op%mesh%ndof,op%mesh%dof,.false.,matrix)
        do e=1,op%mesh%nelem
            element=0
            do q=1,mq
                call metric_at(op%jac(:,:,q,e),op%coef(:,:,q,e), &
                    op%weight(q),op%reaction(q,e),g,mass,det)
                do b=1,nl
                    do a=1,nl
                        element(a,b)=element(a,b)+mass*basis(q,a)*basis(q,b)+ &
                            dot_product(gradient(:,q,a),matmul(g,gradient(:,q,b)))
                    end do
                end do
            end do
            do b=1,nl
                j=op%mesh%dof(b,e)
                do a=1,nl
                    i=op%mesh%dof(a,e)
                    entry=csr_find_entry(matrix,i,j)
                    matrix%val(entry)=matrix%val(entry)+element(a,b)
                end do
            end do
        end do
        do i=1,matrix%n
            do k=matrix%rowptr(i),matrix%rowptr(i+1)-1
                j=matrix%colidx(k)
                if (op%mesh%boundary(i) .or. op%mesh%boundary(j)) then
                    matrix%val(k)=0
                    if (i==j) matrix%val(k)=1
                end if
            end do
        end do
    end subroutine

    subroutine tensor_diagonal(op,diagonal)
        type(tensor_operator_t), intent(in) :: op
        real(dp), allocatable, intent(out) :: diagonal(:)
        integer :: n,m,nl,mq,ws,colour,ce,e,a,r,s,d,q
        n=op%ni(1); m=op%no(1)
        nl=product(op%ni); mq=product(op%no); ws=max(n,m)**3
        allocate(diagonal(op%mesh%ndof))
        diagonal=0
        !$omp parallel default(shared) private(colour,ce,e,a,r,s,d,q)
        block
            real(dp), allocatable :: t1(:),t2(:),table(:,:,:), &
                local(:),term(:),values(:)
            real(dp) :: g(3,3),mass,det,identity(1,1)
            allocate(t1(ws),t2(ws),table(n,m,3),local(nl),term(nl),values(mq))
            identity=1
            do colour=1,op%ncolours
                !$omp do schedule(static)
                do ce=op%colour_ptr(colour),op%colour_ptr(colour+1)-1
                    e=op%colour_cells(ce)
                    local=0
                    do r=1,op%mesh%dim+1
                        do s=1,op%mesh%dim+1
                            if (max(r,s)==op%mesh%dim+1 .and. r/=s) cycle
                            do d=1,3
                                table(:,:,d)=transpose(op%line%b(:,:,1)**2)
                                if (r==d .and. s==d) then
                                    table(:,:,d)=transpose(op%line%d(:,:,1)**2)
                                else if (r==d .or. s==d) then
                                    table(:,:,d)=transpose(op%line%b(:,:,1)* &
                                        op%line%d(:,:,1))
                                end if
                            end do
                            do q=1,mq
                                call metric_at(op%jac(:,:,q,e),op%coef(:,:,q,e), &
                                    op%weight(q),op%reaction(q,e),g,mass,det)
                                if (r==op%mesh%dim+1) then
                                    values(q)=mass
                                else
                                    values(q)=g(r,s)
                                end if
                            end do
                            if (op%mesh%dim==2) then
                                call tensor_contract(op%no,op%ni, &
                                    table(:,:,1),table(:,:,2),identity, &
                                    values,term,t1,t2)
                            else
                                call tensor_contract(op%no,op%ni, &
                                    table(:,:,1),table(:,:,2),table(:,:,3), &
                                    values,term,t1,t2)
                            end if
                            local=local+term
                        end do
                    end do
                    do a=1,nl
                        diagonal(op%mesh%dof(a,e))= &
                            diagonal(op%mesh%dof(a,e))+local(a)
                    end do
                end do
                !$omp end do
            end do
        end block
        !$omp end parallel
        where(op%mesh%boundary) diagonal=1
    end subroutine

    subroutine tensor_load(op,forcing,rhs)
        type(tensor_operator_t), intent(in) :: op
        procedure(scalar_function) :: forcing
        real(dp), allocatable, intent(out) :: rhs(:)
        real(dp), allocatable :: values(:),local(:),t1(:),t2(:),bt(:,:)
        real(dp) :: identity(1,1)
        integer :: e,q,a,ws
        ws=maxval([op%ni,op%no])**3
        allocate(rhs(op%mesh%ndof),values(product(op%no)), &
            local(product(op%ni)),t1(ws),t2(ws))
        bt=transpose(op%line%b(:,:,1))
        identity=1
        rhs=0
        do e=1,op%mesh%nelem
            do q=1,product(op%no)
                values(q)=op%volume(q,e)*forcing(op%points(:,q,e))
            end do
            if (op%mesh%dim==2) then
                call tensor_contract(op%no,op%ni,bt,bt,identity,values,local,t1,t2)
            else
                call tensor_contract(op%no,op%ni,bt,bt,bt,values,local,t1,t2)
            end if
            do a=1,product(op%ni)
                rhs(op%mesh%dof(a,e))=rhs(op%mesh%dof(a,e))+local(a)
            end do
        end do
        where(op%mesh%boundary) rhs=0
    end subroutine

    function tensor_l2_error(op,x,exact) result(error)
        type(tensor_operator_t), intent(in) :: op
        real(dp), intent(in) :: x(:)
        procedure(scalar_function) :: exact
        real(dp) :: error,identity(1,1)
        real(dp), allocatable :: values(:),local(:),t1(:),t2(:),b(:,:)
        integer :: e,q,ws
        ws=maxval([op%ni,op%no])**3
        allocate(values(product(op%no)),local(product(op%ni)),t1(ws),t2(ws))
        b=op%line%b(:,:,1)
        identity=1
        error=0
        do e=1,op%mesh%nelem
            local=x(op%mesh%dof(:,e))
            if (op%mesh%dim==2) then
                call tensor_contract(op%ni,op%no,b,b,identity,local,values,t1,t2)
            else
                call tensor_contract(op%ni,op%no,b,b,b,local,values,t1,t2)
            end if
            do q=1,product(op%no)
                error=error+op%volume(q,e)*(values(q)-exact(op%points(:,q,e)))**2
            end do
        end do
        error=sqrt(error)
    end function
end module
