module crack_contact_module
    use globalParameters, only: dx, maxfam, friction_angle, &
        contact_damage_threshold, &
        intact_neighbor_damage_threshold, full_damage_threshold, &
        friction_regularization, &
        contact_normal_tolerance, damage_gradient_tolerance
        !contact_damage_threshold：节点被认为出现损伤的最低阈值；
        !intact_cnode_damage_threshold=0.05：损伤小于该值的邻居用于重构完好材料方向；
        !full_damage_threshold≈1：大于该值视为完全损伤，不计算摩擦；
        !friction_regularization：小滑移阶段的摩擦正则化长度。
    use Global_arrays, only: coord, disp, olddisp, pforce, nodefam, fail, &
        post_damage, kexi, bulk_modulus, m_volume, bond_idist, scr, &
        bond_omega, bond_fac, flag_node, flag_gradient_valid, &
        flag_pddo_kernel, flag_spatial_increment_gradient, &
        flag_deformation_inverse
    use scalars, only: totnode, totint, vol, step
    implicit none
    private

    integer, allocatable, save :: selected_index(:)
    real*8, allocatable, save :: crack_normal(:,:)
    logical, allocatable, save :: crack_normal_valid(:)
    logical, save :: contact_initialized = .false.

    public :: initialize_crack_contact
    public :: apply_crack_contact

contains

    !这个子程序主要是初始化分配数组
    subroutine initialize_crack_contact()
        implicit none
        integer :: selected_nnum
        integer :: node

        if (contact_initialized) return
        if (totnode <= 0 .or. totint <= 0) then
            error stop 'Crack contact: node counts must be initialized first'
        endif

        allocate(selected_index(totnode))
        allocate(crack_normal(totint,2))
        allocate(crack_normal_valid(totint))

        selected_index = 0
        !这是全局节点号到 flag_node 局部编号的映射。
        do selected_nnum = 1,size(flag_node)
            node = flag_node(selected_nnum)
            if (node >= 1 .and. node <= totnode) then
                selected_index(node) = selected_nnum
            endif
        enddo

        crack_normal = 0.0d0
        crack_normal_valid = .false.
        contact_initialized = .true.
    end subroutine initialize_crack_contact


    subroutine apply_crack_contact()
        implicit none

        integer :: friction_interaction_count
        integer :: valid_normal_count

        if (.not. contact_initialized) then
            error stop 'Crack contact: initialize_crack_contact was not called'
        endif

        call reconstruct_crack_normals(valid_normal_count)
        call apply_contact_friction(friction_interaction_count)

        if (mod(step,5000) == 0) then
            write(*,'(A,I8,A,I8)') &
                'friction interactions = ',friction_interaction_count, &
                ', valid crack normals = ',valid_normal_count
        endif
    end subroutine apply_crack_contact


    subroutine apply_contact_friction(friction_interaction_count)
        implicit none
        integer, intent(out) :: friction_interaction_count

        integer :: i, family_nnum, idx, j
        real*8 :: Y_vec(2), nlength
        real*8 :: idist, omega, fac
        real*8 :: tV_i, tV_j, normal_pressure
        real*8 :: friction_weight
        logical :: friction_i, friction_j

        friction_interaction_count = 0

        ! The material force already retains compressive tV at damaged
        ! endpoints.  Friction therefore traverses only active reference
        ! family interactions and adds no separate normal penalty force.
        ! This loop remains serial because the PDDO compensation scatters to
        ! neighboring nodes.
        do i = 1,totint
            if (post_damage(i) <= contact_damage_threshold) cycle

            do family_nnum = 1,maxfam
                idx = (i-1)*maxfam+family_nnum
                j = nodefam(idx)
                if (j == 0) exit
                if (j > totint) cycle
                if (fail(idx) == 0) cycle
                if (post_damage(j) <= contact_damage_threshold) cycle

                Y_vec = [coord(j,1)+disp(j,1)-coord(i,1)-disp(i,1), &
                    coord(j,2)+disp(j,2)-coord(i,2)-disp(i,2)]
                nlength = norm2(Y_vec)
                idist = bond_idist(idx)
                
                if (nlength <= 1.0d-15 .or. &
                    nlength >= idist) cycle
                if (.not. nodes_face_each_other(i,j,Y_vec)) cycle

                friction_i = node_can_friction(i,Y_vec)
                if (.not. friction_i) cycle
                friction_j = node_can_friction(j,-Y_vec)
                if (friction_j) then
                    friction_weight = 0.5d0
                else
                    friction_weight = 1.0d0
                endif

                omega = bond_omega(idx)
                fac = bond_fac(idx)
                tV_i = 2.0d0*bulk_modulus(i)*kexi(i)*omega* &
                    idist/m_volume(i)
                tV_j = 2.0d0*bulk_modulus(j)*kexi(j)*omega* &
                    idist/m_volume(j)

                ! Only the compressive portions of the two endpoint force
                ! states supply Coulomb pressure.  fac and scr match the
                ! material-force assembly of this same reference interaction.
                normal_pressure = (max(-tV_i,0.0d0) + &
                    max(-tV_j,0.0d0))*fac*scr(idx)
                if (normal_pressure <= 0.0d0) cycle

                call apply_directed_friction(i,j,Y_vec,normal_pressure, &
                    friction_weight)
                friction_interaction_count = friction_interaction_count+1
            enddo
        enddo
    end subroutine apply_contact_friction


    subroutine reconstruct_crack_normals(valid_normal_count)
        implicit none
        integer, intent(out) :: valid_normal_count
        integer :: i, family_nnum, idx, cnode, selected_nnum
        real*8 :: damage_gradient(2)
        real*8 :: material_kernel(2), spatial_kernel(2)
        real*8 :: damage_difference, gradient_norm

        crack_normal = 0.0d0
        crack_normal_valid = .false.
        valid_normal_count = 0

        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(i,family_nnum,idx,cnode,selected_nnum, &
        !$omp& damage_gradient,material_kernel,spatial_kernel, &
        !$omp& damage_difference,gradient_norm) &
        !$omp& reduction(+:valid_normal_count)
        do i = 1,totint
            if (post_damage(i) <= contact_damage_threshold) cycle
            selected_nnum = selected_index(i)
            if (selected_nnum <= 0) cycle
            if (.not. flag_gradient_valid(selected_nnum)) cycle

            damage_gradient = 0.0d0

            ! The gradient of the smoothed damage field points from intact
            ! material toward the center of the damage band.  Consequently,
            ! the gradients reconstructed on the two crack flanks point
            ! toward one another and provide consistently opposed normals.
            do family_nnum = 1,maxfam
                idx = (i-1)*maxfam+family_nnum
                cnode = nodefam(idx)
                if (cnode == 0) exit
                if (cnode > totint) cycle
                if (fail(idx) == 0) cycle

                material_kernel = &
                    flag_pddo_kernel(selected_nnum,family_nnum,:)

                ! grad_x(D) = F^{-T} grad_X(D).  Written as a row-vector
                ! product, this is grad_X(D)^T F^{-1}, consistent with the
                ! spatial-kernel transformation used by the friction update.
                spatial_kernel(1) = material_kernel(1)* &
                    flag_deformation_inverse(selected_nnum,1,1) + &
                    material_kernel(2)* &
                    flag_deformation_inverse(selected_nnum,2,1)
                spatial_kernel(2) = material_kernel(1)* &
                    flag_deformation_inverse(selected_nnum,1,2) + &
                    material_kernel(2)* &
                    flag_deformation_inverse(selected_nnum,2,2)

                damage_difference = post_damage(cnode)-post_damage(i)
                damage_gradient = damage_gradient + damage_difference* &
                    spatial_kernel*vol*bond_fac(idx)
            enddo

            gradient_norm = sqrt(dot_product(damage_gradient, &
                damage_gradient))
            if (gradient_norm > damage_gradient_tolerance) then
                crack_normal(i,:) = damage_gradient/gradient_norm
                crack_normal_valid(i) = .true.
                valid_normal_count = valid_normal_count+1
            endif
        enddo
        !$omp end parallel do
    end subroutine reconstruct_crack_normals


    logical function nodes_face_each_other(i,j,Y_vec)
        implicit none
        integer, intent(in) :: i, j
        real*8, intent(in) :: Y_vec(2)

        real*8 :: nlength
        real*8 :: direction_ij(2)
        logical :: valid_i, valid_j
        logical :: faces_i, faces_j

        nodes_face_each_other = .false.
        nlength = sqrt(dot_product(Y_vec,Y_vec))
        if (nlength <= 1.0d-15) return

        direction_ij = Y_vec/nlength
        valid_i = crack_normal_valid(i)
        valid_j = crack_normal_valid(j)

        faces_i = valid_i .and. dot_product(crack_normal(i,:), &
            direction_ij) > contact_normal_tolerance
        faces_j = valid_j .and. dot_product(crack_normal(j,:), &
            -direction_ij) > contact_normal_tolerance

        if (valid_i .and. valid_j) then
            ! Both reconstructed outward normals must point toward the other
            ! node and must oppose one another.
            nodes_face_each_other = faces_i .and. faces_j .and. &
                dot_product(crack_normal(i,:),crack_normal(j,:)) < 0.0d0
        elseif (valid_i) then
            ! A fully damaged interior node may not possess a reconstructable
            ! normal; in that case the reliable normal of the opposite,
            ! partially damaged surface node controls the geometric test.
            nodes_face_each_other = faces_i .and. &
                post_damage(j) >= full_damage_threshold
        elseif (valid_j) then
            nodes_face_each_other = faces_j .and. &
                post_damage(i) >= full_damage_threshold
        endif
    end function nodes_face_each_other


    logical function node_can_friction(i,Y_vec)
        implicit none
        integer, intent(in) :: i
        real*8, intent(in) :: Y_vec(2)
        real*8 :: nlength

        node_can_friction = .false.
        if (post_damage(i) <= contact_damage_threshold) return
        if (.not. crack_normal_valid(i)) return
        if (selected_index(i) <= 0) return
        if (.not. flag_gradient_valid(selected_index(i))) return

        nlength = sqrt(dot_product(Y_vec, &
            Y_vec))
        if (nlength <= 1.0d-15) return
        node_can_friction = dot_product(crack_normal(i,:), &
            Y_vec/nlength) > 0.0d0
    end function node_can_friction


    subroutine apply_directed_friction(i,j,Y_vec,normal_pressure, &
        friction_weight)
        implicit none
        integer, intent(in) :: i, j
        real*8, intent(in) :: Y_vec(2)
        real*8, intent(in) :: normal_pressure
        real*8, intent(in) :: friction_weight

        integer :: selected_nnum
        integer :: family_nnum, idx, cnode
        real*8 :: increment_i(2), increment_j(2)
        real*8 :: corrected_increment(2), tangent_increment(2)
        real*8 :: friction_force(2)
        real*8 :: material_kernel(2), spatial_kernel(2)
        real*8 :: tangential_norm
        real*8 :: friction_factor, friction_magnitude, alpha_k

        selected_nnum = selected_index(i)
        if (selected_nnum <= 0) return

        if (normal_pressure <= 0.0d0) return

        increment_i = disp(i,:)-olddisp(i,:)
        increment_j = disp(j,:)-olddisp(j,:)
        corrected_increment = increment_j-increment_i - &
            matmul(flag_spatial_increment_gradient(selected_nnum,:,:), &
            Y_vec)
        tangent_increment = corrected_increment - &
            dot_product(corrected_increment,crack_normal(i,:))* &
            crack_normal(i,:)

        tangential_norm = sqrt(dot_product(tangent_increment, &
            tangent_increment))
        if (tangential_norm <= 1.0d-30) return

        !friction_factor = min(tangential_norm/friction_regularization,1.0d0)
        friction_magnitude = tan(friction_angle)*normal_pressure! * friction_factor
        friction_force = -friction_magnitude*tangent_increment/ &
            tangential_norm

        ! Direct equal/opposite friction contribution.
        pforce(j,:) = pforce(j,:) + &
            friction_weight*friction_force*vol
        pforce(i,:) = pforce(i,:) - &
            friction_weight*friction_force*vol

        ! PDDO material-cnodehood compensation.  The sum of every update
        ! is zero, while the first-order reproduction identity compensates
        ! the moment of the non-central direct friction pair.
        do family_nnum = 1,maxfam
            idx = (i-1)*maxfam+family_nnum
            cnode = nodefam(idx)
            if (cnode == 0) exit
            if (fail(idx) == 0) cycle

            material_kernel = &
                flag_pddo_kernel(selected_nnum,family_nnum,:)
            spatial_kernel(1) = material_kernel(1)* &
                flag_deformation_inverse(selected_nnum,1,1) + &
                material_kernel(2)* &
                flag_deformation_inverse(selected_nnum,2,1)
            spatial_kernel(2) = material_kernel(1)* &
                flag_deformation_inverse(selected_nnum,1,2) + &
                material_kernel(2)* &
                flag_deformation_inverse(selected_nnum,2,2)

            alpha_k = dot_product(spatial_kernel,Y_vec)* &
                vol*bond_fac(idx)
            pforce(cnode,:) = pforce(cnode,:) - &
                friction_weight*friction_force*alpha_k*vol
            pforce(i,:) = pforce(i,:) + &
                friction_weight*friction_force*alpha_k*vol
        enddo
    end subroutine apply_directed_friction

end module crack_contact_module


!module crack_contact_module
!    use globalParameters, only: dx, maxfam, friction_angle, &
!        contact_damage_threshold, &
!        intact_neighbor_damage_threshold, full_damage_threshold, &
!        friction_regularization, &
!        contact_normal_tolerance
!        !contact_damage_threshold：节点被认为出现损伤的最低阈值；
!        !intact_cnode_damage_threshold=0.05：损伤小于该值的邻居用于重构完好材料方向；
!        !full_damage_threshold≈1：大于该值视为完全损伤，不计算摩擦；
!        !friction_regularization：小滑移阶段的摩擦正则化长度。
!    use Global_arrays, only: coord, disp, olddisp, pforce, nodefam, fail, &
!        node_damage, kexi, bulk_modulus, m_volume, bond_idist, scr, &
!        bond_omega, bond_fac, selected_node, &
!        selected_gradient_valid, selected_pddo_kernel, &
!        selected_spatial_increment_gradient, selected_deformation_inverse
!    use scalars, only: totnode, totint, vol, step
!    implicit none
!    private
!
!    integer, allocatable, save :: selected_index(:)
!    real*8, allocatable, save :: crack_normal(:,:)
!    logical, allocatable, save :: crack_normal_valid(:)
!    logical, save :: contact_initialized = .false.
!
!    public :: initialize_crack_contact
!    public :: apply_crack_contact
!
!contains
!
!    !这个子程序主要是初始化分配数组
!    subroutine initialize_crack_contact()
!        implicit none
!        integer :: selected_nnum
!        integer :: node
!
!        if (contact_initialized) return
!        if (totnode <= 0 .or. totint <= 0) then
!            error stop 'Crack contact: node counts must be initialized first'
!        endif
!
!        allocate(selected_index(totnode))
!        allocate(crack_normal(totint,2))
!        allocate(crack_normal_valid(totint))
!
!        selected_index = 0
!        !这是全局节点号到 selected_node 局部编号的映射。
!        do selected_nnum = 1,size(selected_node)
!            node = selected_node(selected_nnum)
!            if (node >= 1 .and. node <= totnode) then
!                selected_index(node) = selected_nnum
!            endif
!        enddo
!
!        crack_normal = 0.0d0
!        crack_normal_valid = .false.
!        contact_initialized = .true.
!    end subroutine initialize_crack_contact
!
!
!    subroutine apply_crack_contact()
!        implicit none
!
!        integer :: friction_interaction_count
!        integer :: valid_normal_count
!
!        if (.not. contact_initialized) then
!            error stop 'Crack contact: initialize_crack_contact was not called'
!        endif
!
!        call reconstruct_crack_normals(valid_normal_count)
!        call apply_contact_friction(friction_interaction_count)
!
!        if (mod(step,5000) == 0) then
!            write(*,'(A,I8,A,I8)') &
!                'friction interactions = ',friction_interaction_count, &
!                ', valid crack normals = ',valid_normal_count
!        endif
!    end subroutine apply_crack_contact
!
!
!    subroutine apply_contact_friction(friction_interaction_count)
!        implicit none
!        integer, intent(out) :: friction_interaction_count
!
!        integer :: i, family_nnum, idx, j
!        real*8 :: Y_vec(2), nlength
!        real*8 :: idist, omega, fac
!        real*8 :: tV_i, tV_j, normal_pressure
!        real*8 :: friction_weight
!        logical :: friction_i, friction_j
!
!        friction_interaction_count = 0
!
!        ! The material force already retains compressive tV at damaged
!        ! endpoints.  Friction therefore traverses only active reference
!        ! family interactions and adds no separate normal penalty force.
!        ! This loop remains serial because the PDDO compensation scatters to
!        ! neighboring nodes.
!        do i = 1,totint
!            if (node_damage(i) <= contact_damage_threshold) cycle
!
!            do family_nnum = 1,maxfam
!                idx = (i-1)*maxfam+family_nnum
!                j = nodefam(idx)
!                if (j == 0) exit
!                if (j > totint) cycle
!                if (fail(idx) == 0) cycle
!                if (node_damage(j) <= contact_damage_threshold) cycle
!
!                Y_vec = [coord(j,1)+disp(j,1)-coord(i,1)-disp(i,1), &
!                    coord(j,2)+disp(j,2)-coord(i,2)-disp(i,2)]
!                nlength = norm2(Y_vec)
!                if (nlength <= 1.0d-15 .or. &
!                    nlength >= dx) cycle
!                if (.not. nodes_face_each_other(i,j,Y_vec)) cycle
!
!                friction_i = node_can_friction(i,Y_vec)
!                if (.not. friction_i) cycle
!                friction_j = node_can_friction(j,-Y_vec)
!                if (friction_j) then
!                    friction_weight = 0.5d0
!                else
!                    friction_weight = 1.0d0
!                endif
!
!                idist = bond_idist(idx)
!                omega = bond_omega(idx)
!                fac = bond_fac(idx)
!                tV_i = 2.0d0*bulk_modulus(i)*kexi(i)*omega* &
!                    idist/m_volume(i)
!                tV_j = 2.0d0*bulk_modulus(j)*kexi(j)*omega* &
!                    idist/m_volume(j)
!
!                ! Only the compressive portions of the two endpoint force
!                ! states supply Coulomb pressure.  fac and scr match the
!                ! material-force assembly of this same reference interaction.
!                normal_pressure = (max(-tV_i,0.0d0) + &
!                    max(-tV_j,0.0d0))*fac*scr(idx)
!                if (normal_pressure <= 0.0d0) cycle
!
!                call apply_directed_friction(i,j,Y_vec,normal_pressure, &
!                    friction_weight)
!                friction_interaction_count = friction_interaction_count+1
!            enddo
!        enddo
!    end subroutine apply_contact_friction
!
!
!    subroutine reconstruct_crack_normals(valid_normal_count)
!        implicit none
!        integer, intent(out) :: valid_normal_count
!        integer :: i, j, idx, cnode
!        real*8 :: Y_vec(2)
!        real*8 :: nlength, normal_norm
!
!        crack_normal = 0.0d0
!        crack_normal_valid = .false.
!        valid_normal_count = 0
!
!        !$omp parallel do default(shared) schedule(static) &
!        !$omp& private(i,j,idx,cnode,Y_vec,nlength,normal_norm) &
!        !$omp& reduction(+:valid_normal_count)
!        do i = 1,totint
!            if (node_damage(i) <= contact_damage_threshold) cycle
!
!            do j = 1,maxfam
!                idx = (i-1)*maxfam+j
!                cnode = nodefam(idx)
!                if (cnode == 0) exit
!                if (cnode > totint) cycle
!                if (fail(idx) == 0) cycle
!                if (node_damage(cnode) > &
!                    intact_neighbor_damage_threshold) cycle
!
!                Y_vec = [coord(cnode,1)+disp(cnode,1)-coord(i,1)-disp(i,1), &
!                            coord(cnode,2)+disp(cnode,2)-coord(i,2)-disp(i,2)]
!                
!                nlength = norm2(Y_vec)
!                if (nlength <= 1.0d-15) cycle
!
!                crack_normal(i,:) = crack_normal(i,:) - Y_vec/nlength
!            enddo
!
!            normal_norm = sqrt(dot_product(crack_normal(i,:), &
!                crack_normal(i,:)))
!            if (normal_norm > contact_normal_tolerance) then
!                crack_normal(i,:) = crack_normal(i,:)/normal_norm
!                crack_normal_valid(i) = .true.
!                valid_normal_count = valid_normal_count+1
!            endif
!        enddo
!        !$omp end parallel do
!    end subroutine reconstruct_crack_normals
!
!
!    logical function nodes_face_each_other(i,j,Y_vec)
!        implicit none
!        integer, intent(in) :: i, j
!        real*8, intent(in) :: Y_vec(2)
!
!        real*8 :: nlength
!        real*8 :: direction_ij(2)
!        logical :: valid_i, valid_j
!        logical :: faces_i, faces_j
!
!        nodes_face_each_other = .false.
!        nlength = sqrt(dot_product(Y_vec,Y_vec))
!        if (nlength <= 1.0d-15) return
!
!        direction_ij = Y_vec/nlength
!        valid_i = crack_normal_valid(i)
!        valid_j = crack_normal_valid(j)
!
!        faces_i = valid_i .and. dot_product(crack_normal(i,:), &
!            direction_ij) > contact_normal_tolerance
!        faces_j = valid_j .and. dot_product(crack_normal(j,:), &
!            -direction_ij) > contact_normal_tolerance
!
!        if (valid_i .and. valid_j) then
!            ! Both reconstructed outward normals must point toward the other
!            ! node and must oppose one another.
!            nodes_face_each_other = faces_i .and. faces_j .and. &
!                dot_product(crack_normal(i,:),crack_normal(j,:)) < 0.0d0
!        elseif (valid_i) then
!            ! A fully damaged interior node may not possess a reconstructable
!            ! normal; in that case the reliable normal of the opposite,
!            ! partially damaged surface node controls the geometric test.
!            nodes_face_each_other = faces_i .and. &
!                node_damage(j) >= full_damage_threshold
!        elseif (valid_j) then
!            nodes_face_each_other = faces_j .and. &
!                node_damage(i) >= full_damage_threshold
!        endif
!    end function nodes_face_each_other
!
!
!    logical function node_can_friction(i,Y_vec)
!        implicit none
!        integer, intent(in) :: i
!        real*8, intent(in) :: Y_vec(2)
!        real*8 :: nlength
!
!        node_can_friction = .false.
!        if (node_damage(i) <= contact_damage_threshold) return
!        if (node_damage(i) >= full_damage_threshold) return
!        if (.not. crack_normal_valid(i)) return
!        if (selected_index(i) <= 0) return
!        if (.not. selected_gradient_valid(selected_index(i))) return
!
!        nlength = sqrt(dot_product(Y_vec, &
!            Y_vec))
!        if (nlength <= 1.0d-15) return
!        node_can_friction = dot_product(crack_normal(i,:), &
!            Y_vec/nlength) > 0.0d0
!    end function node_can_friction
!
!
!    subroutine apply_directed_friction(i,j,Y_vec,normal_pressure, &
!        friction_weight)
!        implicit none
!        integer, intent(in) :: i, j
!        real*8, intent(in) :: Y_vec(2)
!        real*8, intent(in) :: normal_pressure
!        real*8, intent(in) :: friction_weight
!
!        integer :: selected_nnum
!        integer :: family_nnum, idx, cnode
!        real*8 :: increment_i(2), increment_j(2)
!        real*8 :: corrected_increment(2), tangent_increment(2)
!        real*8 :: friction_force(2)
!        real*8 :: material_kernel(2), spatial_kernel(2)
!        real*8 :: tangential_norm
!        real*8 :: friction_factor, friction_magnitude, alpha_k
!
!        selected_nnum = selected_index(i)
!        if (selected_nnum <= 0) return
!
!        if (normal_pressure <= 0.0d0) return
!
!        increment_i = disp(i,:)-olddisp(i,:)
!        increment_j = disp(j,:)-olddisp(j,:)
!        corrected_increment = increment_j-increment_i - &
!            matmul(selected_spatial_increment_gradient(selected_nnum,:,:), &
!            Y_vec)
!        tangent_increment = corrected_increment - &
!            dot_product(corrected_increment,crack_normal(i,:))* &
!            crack_normal(i,:)
!
!        tangential_norm = sqrt(dot_product(tangent_increment, &
!            tangent_increment))
!        if (tangential_norm <= 1.0d-30) return
!
!        !friction_factor = min(tangential_norm/friction_regularization,1.0d0)
!        friction_magnitude = tan(friction_angle)*normal_pressure !* friction_factor
!        friction_force = -friction_magnitude*tangent_increment/ &
!            tangential_norm
!
!        ! Direct equal/opposite friction contribution.
!        pforce(j,:) = pforce(j,:) + &
!            friction_weight*friction_force*vol
!        pforce(i,:) = pforce(i,:) - &
!            friction_weight*friction_force*vol
!
!        ! PDDO material-cnodehood compensation.  The sum of every update
!        ! is zero, while the first-order reproduction identity compensates
!        ! the moment of the non-central direct friction pair.
!        do family_nnum = 1,maxfam
!            idx = (i-1)*maxfam+family_nnum
!            cnode = nodefam(idx)
!            if (cnode == 0) exit
!            if (fail(idx) == 0) cycle
!
!            material_kernel = &
!                selected_pddo_kernel(selected_nnum,family_nnum,:)
!            spatial_kernel(1) = material_kernel(1)* &
!                selected_deformation_inverse(selected_nnum,1,1) + &
!                material_kernel(2)* &
!                selected_deformation_inverse(selected_nnum,2,1)
!            spatial_kernel(2) = material_kernel(1)* &
!                selected_deformation_inverse(selected_nnum,1,2) + &
!                material_kernel(2)* &
!                selected_deformation_inverse(selected_nnum,2,2)
!
!            alpha_k = dot_product(spatial_kernel,Y_vec)* &
!                vol*bond_fac(idx)
!            pforce(cnode,:) = pforce(cnode,:) - &
!                friction_weight*friction_force*alpha_k*vol
!            pforce(i,:) = pforce(i,:) + &
!                friction_weight*friction_force*alpha_k*vol
!        enddo
!    end subroutine apply_directed_friction
!
!end module crack_contact_module
