program test_tetra_vector_fallback_jets
    use fortfem_kinds, only: dp
    use check, only: check_condition, check_summary
    use fortfem_tetra_nedelec_arbitrary_order, only: &
        tetra_nedelec_first_kind_t, initialize_tetra_nedelec_first_kind, &
        tetra_nedelec_dof_count, evaluate_tetra_nedelec_first_kind, &
        evaluate_tetra_nedelec_first_kind_jvp, evaluate_tetra_nedelec_first_kind_vjp
    use fortfem_tetra_rt_arbitrary_order, only: tetra_rt_t, initialize_tetra_rt, &
        tetra_rt_dof_count, evaluate_tetra_rt, evaluate_tetra_rt_jvp, evaluate_tetra_rt_vjp
    implicit none
    type(tetra_nedelec_first_kind_t) :: nedelec
    type(tetra_rt_t) :: rt
    real(dp), parameter :: step=2e-6_dp
    real(dp), parameter :: points(3,8)=reshape([ &
        0._dp,0._dp,0._dp, 1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, &
        0._dp,0._dp,1._dp, .5_dp,.5_dp,0._dp, .25_dp,0._dp,.25_dp, &
        0._dp,.25_dp,.25_dp, .19_dp,.17_dp,.23_dp], [3,8])
    real(dp), allocatable :: values(:,:), curls(:,:), plus(:,:), twice(:,:)
    real(dp), allocatable :: curl_plus(:,:), curl_twice(:,:), values_dot(:,:), curls_dot(:,:)
    real(dp), allocatable :: divergences(:), div_plus(:), div_twice(:), div_dot(:)
    real(dp), allocatable :: values_bar(:,:), curls_bar(:,:), div_bar(:)
    real(dp) :: direction(3), point_bar(3), lhs, rhs
    integer :: degree, family, count, point_id, status, k
    do degree=5,6
        do family=1,2
            if(family==1)then
                call initialize_tetra_nedelec_first_kind(degree,nedelec,status)
                count=tetra_nedelec_dof_count(nedelec)
            else
                call initialize_tetra_rt(degree,rt,status)
                count=tetra_rt_dof_count(rt)
            end if
            call check_condition(status==0,'Unbounded fallback basis initializes at degree five/six')
            if(status/=0)cycle
            allocate(values(3,count),curls(3,count),plus(3,count),twice(3,count), &
                curl_plus(3,count),curl_twice(3,count),values_dot(3,count),curls_dot(3,count), &
                divergences(count),div_plus(count),div_twice(count),div_dot(count), &
                values_bar(3,count),curls_bar(3,count),div_bar(count))
            do k=1,count
                values_bar(:,k)=[.013_dp*k,-.009_dp*k,.004_dp*k]
                curls_bar(:,k)=[-.007_dp*k,.003_dp*k,.011_dp*k]
                div_bar(k)=.005_dp*k
            end do
            do point_id=1,8
                ! Inward direction keeps one-sided differences inside the reference tetrahedron.
                direction=[.21_dp,.17_dp,.13_dp]-points(:,point_id)
                if(family==1)then
                    call evaluate_tetra_nedelec_first_kind(nedelec,points(:,point_id),values,curls,status)
                    call check_condition(status==0,'Nedelec values/curls include vertices and zero axes')
                    call evaluate_tetra_nedelec_first_kind(nedelec,points(:,point_id)+step*direction, &
                        plus,curl_plus,status)
                    call evaluate_tetra_nedelec_first_kind(nedelec,points(:,point_id)+2*step*direction, &
                        twice,curl_twice,status)
                    call evaluate_tetra_nedelec_first_kind_jvp(nedelec,points(:,point_id),direction, &
                        values_dot,curls_dot,status)
                    call check_condition(status==0,'Nedelec tangent includes boundary points')
                    call check_fd(values_dot,(-3*values+4*plus-twice)/(2*step),'Nedelec value tangent')
                    call check_fd(curls_dot,(-3*curls+4*curl_plus-curl_twice)/(2*step), &
                        'Nedelec curl tangent')
                    call evaluate_tetra_nedelec_first_kind_vjp(nedelec,points(:,point_id), &
                        values_bar,curls_bar,point_bar,status)
                    lhs=sum(values_bar*values_dot)+sum(curls_bar*curls_dot)
                else
                    call evaluate_tetra_rt(rt,points(:,point_id),values,divergences,status)
                    call check_condition(status==0,'RT values/divergence include vertices and zero axes')
                    call evaluate_tetra_rt(rt,points(:,point_id)+step*direction,plus,div_plus,status)
                    call evaluate_tetra_rt(rt,points(:,point_id)+2*step*direction,twice,div_twice,status)
                    call evaluate_tetra_rt_jvp(rt,points(:,point_id),direction,values_dot,div_dot,status)
                    call check_condition(status==0,'RT tangent includes boundary points')
                    call check_fd(values_dot,(-3*values+4*plus-twice)/(2*step),'RT value tangent')
                    call check_fd(reshape(div_dot,[1,count]), &
                        reshape((-3*divergences+4*div_plus-div_twice)/(2*step),[1,count]), &
                        'RT divergence tangent')
                    call evaluate_tetra_rt_vjp(rt,points(:,point_id),values_bar,div_bar,point_bar,status)
                    lhs=sum(values_bar*values_dot)+sum(div_bar*div_dot)
                end if
                rhs=dot_product(point_bar,direction)
                call check_condition(status==0,'Fallback reverse product succeeds')
                call check_condition(abs(lhs-rhs)<2e-11_dp*max(1._dp,abs(lhs),abs(rhs)), &
                    'Fallback reverse product is adjoint to its directional tangent')
            end do
            deallocate(values,curls,plus,twice,curl_plus,curl_twice,values_dot,curls_dot, &
                divergences,div_plus,div_twice,div_dot,values_bar,curls_bar,div_bar)
        end do
    end do
    call check_summary('Tetrahedral monomial and Koornwinder fallback jets')
contains
    subroutine check_fd(tangent,finite_difference,label)
        real(dp), intent(in) :: tangent(:,:),finite_difference(:,:)
        character(*), intent(in) :: label
        real(dp) :: scale,error
        scale=max(1._dp,maxval(abs(tangent)),maxval(abs(finite_difference)))
        error=maxval(abs(tangent-finite_difference))/scale
        if(error>=2e-6_dp)print *,degree,family,point_id,label,'relative FD error',error
        call check_condition(error<2e-6_dp,label//' matches independent inward finite differences')
    end subroutine check_fd
end program test_tetra_vector_fallback_jets
