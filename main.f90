program mixed_mode_uniaxial_compression_2d
    use globalParameters
    use Global_arrays
    use scalars
    use omp_lib
    use displacement_gradient_module, only: build_flag_node_list_from_flag, &
        calculate_selected_node_displacement_gradients
    use crack_contact_module, only: initialize_crack_contact, &
        apply_crack_contact
    implicit none
    character(len=64) :: filename

    ! ================================================================
    ! Allocation and initialization
    ! ================================================================
    allocate(grid_node(ndivx,ndivy+2*nbnd))
    allocate(coord(ntotmax,2),disp(ntotmax,2),olddisp(ntotmax,2))
    allocate(velhalf(ntotmax,2),velhalfold(ntotmax,2))
    allocate(pforce(ntotmax,2),pforceold(ntotmax,2),bforce(ntotmax,2))
    allocate(massvec(ntotmax,2))
    allocate(nodefam(ntotmax*maxfam),fail(ntotmax*maxfam))
    allocate(col(ntotmax),flag(ntotmax,1))
    allocate(bond_idist(ntotmax*maxfam),bond_omega(ntotmax*maxfam))
    allocate(bond_fac(ntotmax*maxfam),scr(ntotmax*maxfam))
    allocate(emod(ntotmax),pratio(ntotmax),tensile_strength(ntotmax))
    allocate(shear_strength(ntotmax),bulk_modulus(ntotmax))
    allocate(shear_modulus(ntotmax),G1(ntotmax),G2(ntotmax),bc(ntotmax))
    allocate(kexi(ntotmax),m_volume(ntotmax),vol_all(ntotmax))
    allocate(node_damage(ntotmax),post_damage(ntotmax))
    allocate(damage_tension(ntotmax),damage_shear(ntotmax))
    allocate(failure_index_I(ntotmax),failure_index_II(ntotmax))
    allocate(damage_initiation_lambda(ntotmax))
    allocate(equivalent_normal_strain_output(ntotmax))
    allocate(equivalent_shear_strain_output(ntotmax))
    allocate(principal_stress_1_output(ntotmax))
    allocate(principal_stress_3_output(ntotmax))
    allocate(damage_activated(ntotmax))
    allocate(tension_activated(ntotmax),shear_activated(ntotmax))
    allocate(stress_xx(ntotmax),stress_yy(ntotmax),stress_xy(ntotmax))

    grid_node = 0
    coord = 0.0d0
    disp = 0.0d0
    olddisp = 0.0d0
    velhalf = 0.0d0
    velhalfold = 0.0d0
    pforce = 0.0d0
    pforceold = 0.0d0
    bforce = 0.0d0
    massvec = 0.0d0
    nodefam = 0
    fail = 1
    col = 1
    flag = 0
    bond_idist = 0.0d0
    bond_omega = 0.0d0
    bond_fac = 0.0d0
    scr = 0.0d0
    emod = 0.0d0
    pratio = 0.0d0
    tensile_strength = 0.0d0
    shear_strength = 0.0d0
    bulk_modulus = 0.0d0
    shear_modulus = 0.0d0
    G1 = 0.0d0
    G2 = 0.0d0
    bc = 0.0d0
    kexi = 0.0d0
    m_volume = 0.0d0
    vol_all = 0.0d0
    node_damage = 0.0d0
    post_damage = 0.0d0
    damage_tension = 0.0d0
    damage_shear = 0.0d0
    failure_index_I = 0.0d0
    failure_index_II = 0.0d0
    damage_initiation_lambda = 0.0d0
    equivalent_normal_strain_output = 0.0d0
    equivalent_shear_strain_output = 0.0d0
    principal_stress_1_output = 0.0d0
    principal_stress_3_output = 0.0d0
    damage_activated = .false.
    tension_activated = .false.
    shear_activated = .false.
    stress_xx = 0.0d0
    stress_yy = 0.0d0
    stress_xy = 0.0d0
    vol = particle_volume
    radij = particle_radius
    thick = thickness
    cfri = tan(friction_angle)

    if (damage_mode < 1 .or. damage_mode > 3) then
        error stop 'damage_mode must be 1, 2 or 3'
    endif

    ! ================================================================
    ! Geometry: internal plate with a node-deletion flaw
    ! ================================================================
    nnum = 0
    do iy = 1,ndivy
        coordy = (dble(iy)-0.5d0*(dble(ndivy)+1.0d0))*dx
        do ix = 1,ndivx
            coordx = (dble(ix)-0.5d0*(dble(ndivx)+1.0d0))*dx
            along_flaw = coordx*cos(flaw_angle) + coordy*sin(flaw_angle)
            across_flaw = -coordx*sin(flaw_angle) + coordy*cos(flaw_angle)

            ! Points inside the finite-width flaw are not generated.
            if (abs(along_flaw) <= 0.5d0*flaw_length .and. &
                abs(across_flaw) <= 0.5d0*flaw_width) cycle

            nnum = nnum + 1
            coord(nnum,1) = coordx
            coord(nnum,2) = coordy
            grid_node(ix,iy+nbnd) = nnum
        enddo
    enddo
    totint = nnum

    ! Three lower virtual boundary layers.
    do iy = 1,nbnd
        coordy = -0.5d0*plate_height - (dble(iy)-0.5d0)*dx
        do ix = 1,ndivx
            coordx = (dble(ix)-0.5d0*(dble(ndivx)+1.0d0))*dx
            nnum = nnum + 1
            coord(nnum,1) = coordx
            coord(nnum,2) = coordy
            grid_node(ix,nbnd-iy+1) = nnum
        enddo
    enddo
    totbottom = nnum

    ! Three upper virtual boundary layers.
    do iy = 1,nbnd
        coordy = 0.5d0*plate_height + (dble(iy)-0.5d0)*dx
        do ix = 1,ndivx
            coordx = (dble(ix)-0.5d0*(dble(ndivx)+1.0d0))*dx
            nnum = nnum + 1
            coord(nnum,1) = coordx
            coord(nnum,2) = coordy
            grid_node(ix,ndivy+nbnd+iy) = nnum
        enddo
    enddo

    tottop = nnum
    totnode = nnum

    if (totnode > ntotmax) error stop 'Generated node count exceeds ntotmax'

    ! ================================================================
    ! Assign initial rock material properties
    ! ================================================================

    do i = 1, totnode
        emod(i) = young_modulus
        pratio(i) = poisson_ratio
        shear_strength(i) = cohesion
        tensile_strength(i) = tensile_strength_value

        shear_modulus(i) = emod(i) / (2.0d0 * (1.0d0 + pratio(i)))
        bulk_modulus(i) = emod(i) / (2.0d0 * (1.0d0 - pratio(i)))
        G1(i) = fracture_energy_I
        G2(i) = fracture_energy_II
        bc(i) = 9.0d0 * emod(i) / (pi * thick * delta**3)
    enddo
    ! ================================================================
    ! Family search, partial-volume factor, weighted volume and shape
    ! ================================================================
    do jj = 1,ndivy+2*nbnd
        do ii = 1,ndivx
            i = grid_node(ii,jj)
            if (i == 0) cycle

            do dj = -4,4
                if (jj + dj < 1 .or. jj + dj > ndivy+2*nbnd) cycle
                
                do di = -4,4
                    if (ii + di < 1 .or. ii + di > ndivx) cycle
                    
                    cnode = grid_node(ii + di,jj + dj)
                    if (cnode <= i) cycle
                
                    temp = (coord(cnode,1) - coord(i,1))**2 + &
                        (coord(cnode,2) - coord(i,2))**2

                    if (temp <= delta**2) then
                        if (col(i) > maxfam .or. col(cnode) > maxfam) then
                            error stop 'maxfam is too small for the selected horizon'
                        endif
                    
                        nodefam((i-1)*maxfam+col(i)) = cnode
                        nodefam((cnode-1)*maxfam+col(cnode)) = i
                        col(i) = col(i) + 1
                        col(cnode) = col(cnode) + 1
                    endif
                enddo
            enddo
        enddo
    enddo

    ! ================================================================
    ! Break bonds that cross the finite initial-flaw centerline
    ! ================================================================
    do i = 1,totint
        along_i = coord(i,1)*cos(flaw_angle) + &
            coord(i,2)*sin(flaw_angle)
        across_i = -coord(i,1)*sin(flaw_angle) + &
            coord(i,2)*cos(flaw_angle)

        do j = 1,maxfam
            idx = (i-1)*maxfam+j
            cnode = nodefam(idx)
            if (cnode == 0) exit
            if (cnode > totint) cycle

            along_j = coord(cnode,1)*cos(flaw_angle) + &
                coord(cnode,2)*sin(flaw_angle)
            across_j = -coord(cnode,1)*sin(flaw_angle) + &
                coord(cnode,2)*cos(flaw_angle)

            if (across_i*across_j < 0.0d0) then
                crack_cross_parameter = across_i/(across_i-across_j)
                crack_cross_along = along_i + crack_cross_parameter * &
                    (along_j-along_i)

                if (crack_cross_parameter > 0.0d0 .and. &
                    crack_cross_parameter < 1.0d0 .and. &
                    abs(crack_cross_along) <= 0.5d0*flaw_length) then
                    fail(idx) = 0
                endif
            endif
        enddo
    enddo
    ! ================================================================
    ! 标记节点，为后续计算所有节点的PDDO微分算子做铺垫
    ! ================================================================
    do i = 1,totnode
        flag(i,1) = 1
    enddo
    call build_flag_node_list_from_flag(flag_node)
    allocate(flag_gradient(size(flag_node),2,2))
    allocate(flag_gradient_valid(size(flag_node)))
    allocate(flag_pddo_kernel(size(flag_node),maxfam,2))
    allocate(flag_spatial_increment_gradient(size(flag_node),2,2))
    allocate(flag_deformation_inverse(size(flag_node),2,2))
    flag_gradient = 0.0d0
    flag_gradient_valid = .false.
    flag_pddo_kernel = 0.0d0
    flag_spatial_increment_gradient = 0.0d0
    flag_deformation_inverse = 0.0d0
    call initialize_crack_contact()
    
    ! ================================================================
    ! 计算表面修正系数
    ! ================================================================

    vol_all = 0.0d0

    do i = 1,totnode
        do j = 1,maxfam
            idx = (i-1)*maxfam+j
            cnode = nodefam(idx)
            if (cnode == 0) exit
        
            bond_idist(idx) = dsqrt((coord(cnode,1) - coord(i,1))**2 + (coord(cnode,2) - coord(i,2))**2)

            if (bond_idist(idx) <= delta - radij) then
                bond_fac(idx) = 1.0d0
            elseif (bond_idist(idx) <= delta + radij) then
                bond_fac(idx) = (delta + radij - bond_idist(idx)) / (2.0d0 * radij)
            else
                bond_fac(idx) = 0.0d0
            endif

            bond_omega(idx) = dexp(-4.0d0 * (bond_idist(idx)**2) / (delta**2))

            vol_all(i) = vol_all(i) + vol * bond_fac(idx)

        enddo
    enddo
    
    ! 2D surface correction coefficient
    do i = 1,totnode
        do j = 1,maxfam
            idx = (i-1)*maxfam+j
            cnode = nodefam(idx)
            if (cnode == 0) exit

            scr(idx) = (2.0d0 * dx * pi * delta**2) / &
                (vol_all(i) + vol_all(cnode))

        enddo
    enddo

    ! ================================================================
    ! Calculate weighted volume
    ! ================================================================

    do i = 1,totnode
        do j = 1,maxfam
            idx = (i-1)*maxfam+j
            cnode = nodefam(idx)
            if (cnode == 0) exit

            idist = bond_idist(idx)

            m_volume(i) = m_volume(i) + bond_omega(idx) * &
                bond_fac(idx) * vol * idist**2

        enddo
    enddo

    ! ================================================================
    ! Fictitious mass used by adaptive dynamic relaxation
    ! ================================================================
    do i = 1, totnode
        massvec(i,1) = 0.25d0 * dt**2 * pi * delta**2 * thick * bc(i) / dx * 5.0d0
        massvec(i,2) = 0.25d0 * dt**2 * pi * delta**2 * thick * bc(i) / dx * 5.0d0
    enddo

    ! ================================================================
    ! Quasi-static loading loop
    ! ================================================================
    section_area = plate_width * thickness
    open(newunit=stress_history_unit, file='stress_strain_history.txt', &
        status='replace', action='write')

    do step = 1,nsteps

        write(*,*) "tt=", step
        
        ! Symmetric prescribed compression.  Boundary y velocities remain
        ! zero because these are displacement-controlled DOFs; only their x
        ! DOFs participate in the subsequent ADR displacement iteration.
        do i = totint+1, totbottom
            disp(i,2) = disp(i,2) + load_increment*dt
            velhalf(i,2) = 0.0d0
            velhalfold(i,2) = 0.0d0
        enddo

        do i = totbottom+1, tottop
            disp(i,2) = disp(i,2) - load_increment*dt
            velhalf(i,2) = 0.0d0
            velhalfold(i,2) = 0.0d0
        enddo

        pforce = 0.0d0
        kexi = 0.0d0
    
        !$omp parallel private(i,j,idx,cnode,Y_vec,nlength,idist,omega,fac)
        !$omp do schedule(auto)
        do i = 1,totnode
            do j = 1,maxfam
                idx = (i-1)*maxfam+j
                cnode = nodefam(idx)
                if (cnode == 0) exit
                if (fail(idx) == 0) cycle

                Y_vec = [coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1), &
                            coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)]

                nlength = norm2(Y_vec)
                idist = bond_idist(idx)
                omega = bond_omega(idx)
                fac = bond_fac(idx)

                kexi(i) = kexi(i) + (2.0d0 / m_volume(i)) * omega * &
                    (nlength - idist) * idist * fac * vol

            enddo
        enddo
        !$omp end parallel
        
        call calculate_selected_node_displacement_gradients(flag_node, &
            flag_gradient, flag_gradient_valid, flag_pddo_kernel, &
            flag_spatial_increment_gradient, flag_deformation_inverse)

        ! Reconstruct the undamaged (effective) plane-stress tensor at every
        ! PDDO node.  Damage initiation must be driven by this trial stress,
        ! not by the already degraded force state.
        stress_xx = 0.0d0
        stress_yy = 0.0d0
        stress_xy = 0.0d0
        failure_index_I = 0.0d0
        failure_index_II = 0.0d0
        equivalent_normal_strain_output = 0.0d0
        equivalent_shear_strain_output = 0.0d0
        principal_stress_1_output = 0.0d0
        principal_stress_3_output = 0.0d0
        
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(flag_nnum,i,epsilon_xx,epsilon_yy,epsilon_xy, &
        !$omp& mean_strain,strain_radius,mean_stress,stress_radius, &
        !$omp& principal_strain_1,principal_strain_3, &
        !$omp& principal_stress_1,principal_stress_3, &
        !$omp& equivalent_normal_strain,equivalent_shear_strain, &
        !$omp& normal_initial_strain,shear_initial_strain, &
        !$omp& normal_failure_strain,shear_failure_strain, &
        !$omp& rankine_index,mohr_coulomb_index,mc_numerator, &
        !$omp& normalized_normal_strain,normalized_shear_strain, &
        !$omp& normalized_normal_initial,normalized_shear_initial, &
        !$omp& normal_mode_measure,shear_mode_measure,lambda,lambda_p, &
        !$omp& cos_theta,sin_theta,damage_trial,damage_old,damage_new, &
        !$omp& damage_increment,damage_lambda,damage_lambda_p,mode_weight_sum)
        do flag_nnum = 1, size(flag_node)
            if (.not. flag_gradient_valid(flag_nnum)) cycle
            
            i = flag_node(flag_nnum)
            
            !节点i的xx，yy，xy应变
            epsilon_xx = flag_gradient(flag_nnum,1,1)
            epsilon_yy = flag_gradient(flag_nnum,2,2)
            epsilon_xy = 0.5d0 * (flag_gradient(flag_nnum,1,2) + &
             flag_gradient(flag_nnum,2,1))

            stress_xx(i) = emod(i) / (1.0d0-pratio(i)**2) * &
                (epsilon_xx + pratio(i)*epsilon_yy)
            stress_yy(i) = emod(i) / (1.0d0-pratio(i)**2) * &
                (epsilon_yy + pratio(i)*epsilon_xx)
            stress_xy(i) = 2.0d0*shear_modulus(i)*epsilon_xy
        
            !根据节点i的应变计算最大主应变和最小主应变
            mean_strain = 0.5d0 * (epsilon_xx + epsilon_yy)
            strain_radius = dsqrt((0.5d0*(epsilon_xx-epsilon_yy))**2 + epsilon_xy**2)
            principal_strain_1 = mean_strain + strain_radius
            principal_strain_3 = mean_strain - strain_radius
        
            !等效张拉应变εeq设置为ε1，等效剪切应变γeq设置为ε1-ε3
            equivalent_normal_strain = principal_strain_1
            equivalent_shear_strain = principal_strain_1 - principal_strain_3

            ! Principal effective stresses.  Rankine and Mohr-Coulomb only
            ! initiate their respective damage mechanisms; strain controls
            ! the subsequent softening progress.
            mean_stress = 0.5d0 * (stress_xx(i) + stress_yy(i))
            stress_radius = dsqrt((0.5d0*(stress_xx(i)-stress_yy(i)))**2 + &
                stress_xy(i)**2)
            principal_stress_1 = mean_stress + stress_radius
            principal_stress_3 = mean_stress - stress_radius

            equivalent_normal_strain_output(i) = equivalent_normal_strain
            equivalent_shear_strain_output(i) = equivalent_shear_strain
            principal_stress_1_output(i) = principal_stress_1
            principal_stress_3_output(i) = principal_stress_3

            rankine_index = max(principal_stress_1,0.0d0) / &
                tensile_strength(i)

            ! Plane-stress, tension-positive Mohr-Coulomb index:
            ! (sigma1-sigma3)+(sigma1+sigma3)sin(phi)
            ! ------------------------------------------------- >= 1.
            !                  2 c cos(phi)
            mc_numerator = (principal_stress_1-principal_stress_3) + &
                (principal_stress_1+principal_stress_3)*sin(friction_angle)
            mohr_coulomb_index = max(mc_numerator,0.0d0) / &
                (2.0d0*shear_strength(i)*cos(friction_angle))

            failure_index_I(i) = rankine_index
            failure_index_II(i) = mohr_coulomb_index

            ! ============================================================
            ! Restored original strain-based mixed-mode damage criterion.
            ! Stress and failure indices above are diagnostic outputs only.
            ! ============================================================
            normal_initial_strain = tensile_strength(i) / emod(i)
            normal_failure_strain = G1(i) / (tensile_strength(i) * delta)
            shear_initial_strain = shear_strength(i) / shear_modulus(i)
            shear_failure_strain = G2(i) / (shear_strength(i) * delta)

            if (normal_failure_strain <= normal_initial_strain .or. &
                shear_failure_strain <= shear_initial_strain) then
                error stop 'Damage: failure strain must exceed initiation strain'
            endif

            normalized_normal_strain = equivalent_normal_strain / &
                normal_failure_strain
            normalized_shear_strain = equivalent_shear_strain / &
                shear_failure_strain
            normalized_normal_initial = normal_initial_strain / &
                normal_failure_strain
            normalized_shear_initial = shear_initial_strain / &
                shear_failure_strain
            normal_mode_measure = max(normalized_normal_strain,0.0d0)
            shear_mode_measure = max(normalized_shear_strain,0.0d0)

            select case (damage_mode)
            case (1)
                lambda = normal_mode_measure
                lambda_p = normalized_normal_initial
            case (2)
                lambda = shear_mode_measure
                lambda_p = normalized_shear_initial
            case (3)
                lambda = dsqrt(normal_mode_measure**2 + &
                    shear_mode_measure**2)
                if (lambda > tiny_value) then
                    cos_theta = shear_mode_measure/lambda
                    sin_theta = normal_mode_measure/lambda
                    lambda_p = normalized_normal_initial * &
                        normalized_shear_initial / dsqrt( &
                        normalized_normal_initial**2*cos_theta**2 + &
                        normalized_shear_initial**2*sin_theta**2)
                else
                    lambda_p = 0.0d0
                endif
            case default
                error stop 'damage_mode must be 1, 2 or 3'
            end select

            ! When lambda is only round-off noise, the mode angle is
            ! undefined and lambda_p was set to zero above.  Do not pass that
            ! state into the softening equation: lambda_p=0 would otherwise
            ! produce damage_trial=1 for any tiny positive lambda.
            if (lambda <= tiny_value) then
                damage_trial = 0.0d0
            elseif (lambda <= lambda_p) then
                damage_trial = 0.0d0
            elseif (lambda <= 1.0d0) then
                damage_trial = 1.0d0 - (1.0d0-lambda)*lambda_p / &
                    (lambda*(1.0d0-lambda_p))
            else
                damage_trial = 1.0d0
            endif

            damage_old = node_damage(i)
            damage_new = max(damage_old, &
                max(0.0d0,min(1.0d0,damage_trial)))
            damage_increment = damage_new-damage_old
            node_damage(i) = damage_new

            ! The former stress-initiated unified-lambda implementation is
            ! retained below for reference and is intentionally inactive.
            !
            ! normal_failure_strain = G1(i) / &
            !     (tensile_strength(i) * delta)
            ! shear_failure_strain = G2(i) / &
            !     (shear_strength(i) * delta)
            ! normalized_normal_strain = 0.0d0
            ! normalized_shear_strain = 0.0d0
            ! if ((damage_mode == 1 .or. damage_mode == 3) .and. &
            !     normal_failure_strain > tiny_value) then
            !     normalized_normal_strain = equivalent_normal_strain / &
            !         normal_failure_strain
            ! endif
            ! if ((damage_mode == 2 .or. damage_mode == 3) .and. &
            !     shear_failure_strain > tiny_value) then
            !     normalized_shear_strain = equivalent_shear_strain / &
            !         shear_failure_strain
            ! endif
            ! damage_lambda = dsqrt(normalized_normal_strain**2 + &
            !     normalized_shear_strain**2)
            ! if ((damage_mode == 1 .or. damage_mode == 3) .and. &
            !     .not. tension_activated(i) .and. &
            !     rankine_index >= 1.0d0) then
            !     tension_activated(i) = .true.
            ! endif
            ! if ((damage_mode == 2 .or. damage_mode == 3) .and. &
            !     .not. shear_activated(i) .and. &
            !     mohr_coulomb_index >= 1.0d0) then
            !     shear_activated(i) = .true.
            ! endif
            ! if (.not. damage_activated(i) .and. &
            !     (tension_activated(i) .or. shear_activated(i))) then
            !     damage_activated(i) = .true.
            !     damage_initiation_lambda(i) = damage_lambda
            ! endif
            ! if (damage_activated(i)) then
            !     damage_lambda_p = damage_initiation_lambda(i)
            !     if (damage_lambda >= 1.0d0) then
            !         damage_trial = 1.0d0
            !     elseif (damage_lambda <= damage_lambda_p .or. &
            !         damage_lambda_p >= 1.0d0-tiny_value) then
            !         damage_trial = 0.0d0
            !     else
            !         damage_trial = 1.0d0 - &
            !             (1.0d0-damage_lambda)*damage_lambda_p / &
            !             (damage_lambda*(1.0d0-damage_lambda_p))
            !     endif
            !     node_damage(i) = max(node_damage(i), &
            !         max(0.0d0,min(1.0d0,damage_trial)))
            ! endif
            ! mode_weight_sum = normalized_normal_strain**2 + &
            !     normalized_shear_strain**2
            ! if (mode_weight_sum > tiny_value) then
            !     damage_tension(i) = node_damage(i) * &
            !         normalized_normal_strain**2 / mode_weight_sum
            !     damage_shear(i) = node_damage(i) * &
            !         normalized_shear_strain**2 / mode_weight_sum
            ! else
            !     damage_tension(i) = 0.0d0
            !     damage_shear(i) = 0.0d0
            ! endif
        enddo
        !$omp end parallel do
        
        ! 后处理节点损伤取其同材料近场点本构损伤的算术平均值。
        ! node_damage 保留为本构变量，避免空间平均反过来改变键力。
        do i = 1, totnode
            damage_sum = 0.0d0
            nearfield_count = 0

            do j = 1, maxfam
                idx = (i-1)*maxfam+j
                cnode = nodefam(idx)
                if (cnode == 0) exit

                ! Do not smooth across a pre-broken material bond.
                if (fail(idx) == 0) cycle

                damage_sum = damage_sum+node_damage(cnode)
                nearfield_count = nearfield_count+1
            enddo

            if (nearfield_count > 0) then
                post_damage(i) = damage_sum/dble(nearfield_count)
            else
                post_damage(i) = node_damage(i)
            endif
        enddo
        

        if (size(flag_node) > 0 .and. mod(step,5000) == 0) then
        write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,I8)') &
            'max total/tensile/shear damage = ', &
            maxval(node_damage(flag_node)), ', ', &
            maxval(damage_tension(flag_node)), ', ', &
            maxval(damage_shear(flag_node)), ', damaged nodes = ', &
            count(node_damage(flag_node) > 0.0d0)
        endif
    
    
        !$omp parallel private(i,j,idx,cnode,Y_vec,idist,nlength,omega,fac, &
        !$omp& damage_i,damage_j,damage_tension_i,damage_tension_j, &
        !$omp& tV_i,tS_i,tF_i,tV_j,tS_j,tF_j, &
        !$omp& t_i,t_j,dforce)
        !$omp do schedule(static)
        do i = 1,totnode
            do j = 1,maxfam
                idx = (i-1)*maxfam+j
                cnode = nodefam(idx)
                if (cnode == 0) exit
                if (fail(idx) == 0) cycle
            

                Y_vec = [coord(cnode,1) + disp(cnode,1) - coord(i,1) - disp(i,1), &
                            coord(cnode,2) + disp(cnode,2) - coord(i,2) - disp(i,2)]

                nlength = norm2(Y_vec)
                idist = bond_idist(idx)
                omega = bond_omega(idx)
                fac = bond_fac(idx)

                if (nlength <= 1.0d-15) cycle

                damage_i = node_damage(i)
                damage_j = node_damage(cnode)
                damage_tension_i = damage_tension(i)
                damage_tension_j = damage_tension(cnode)
                
                ! Volumetric and deviatoric scalar force states at endpoint i.
                tV_i = 2.0d0 * bulk_modulus(i) * kexi(i) * omega * &
                    idist / m_volume(i)
                tS_i = 8.0d0 * shear_modulus(i) / m_volume(i) * omega * &
                    ((nlength-idist)-0.5d0*kexi(i)*idist)

                ! Endpoint cnode is evaluated independently.
                tV_j = 2.0d0 * bulk_modulus(cnode) * kexi(cnode) * &
                    omega * idist / m_volume(cnode)
                tS_j = 8.0d0 * shear_modulus(cnode) / &
                    m_volume(cnode) * omega * &
                    ((nlength-idist)-0.5d0*kexi(cnode)*idist)

                ! Requested split damage force law:
                ! kexi > 0: T=(1-D_I)T_V+(1-D_total)T_S
                ! kexi < 0: T=T_V+(1-D_total)T_S+D_total*T_F
                if (kexi(i) > 0.0d0) then
                    t_i = (1.0d0-damage_i)*(tV_i+tS_i)*(Y_vec/nlength)
                else
                    t_i = (tV_i+tS_i*(1.0d0-damage_i))*(Y_vec/nlength)
                endif
                
                if (kexi(cnode) > 0.0d0) then
                    t_j = (1.0d0-damage_j)*(tV_j+tS_j)*(-Y_vec/nlength)
                else
                    t_j = (tV_j+tS_j*(1.0d0-damage_j))*(-Y_vec/nlength)
                endif
             

                ! Scalar force states act along the current bond direction.
                dforce = (t_i-t_j) * vol * fac * scr(idx)
                pforce(i,:) = pforce(i,:) + dforce
            enddo
        enddo
        !$omp end parallel

        ! Contact/friction uses post_damage only; node_damage remains the
        ! constitutive variable used by the material force state above.
        call apply_crack_contact()

        ! ================================================================
        ! Global stress-strain response from virtual-boundary reactions
        ! ================================================================
        bottom_reaction = sum(pforce(totint+1:totbottom,2)) * &
            particle_volume
        top_reaction = sum(pforce(totbottom+1:tottop,2)) * &
            particle_volume

        axial_stress_mpa = 0.5d0 * (abs(bottom_reaction) + &
            abs(top_reaction)) / section_area

        axial_strain_percent = 100.0d0 * ( &
            sum(disp(totint+1:totbottom,2)) / &
            dble(totbottom-totint) - &
            sum(disp(totbottom+1:tottop,2)) / &
            dble(tottop-totbottom)) / plate_height

        boundary_unbalance = abs(abs(top_reaction) - &
            abs(bottom_reaction)) / max(0.5d0 * &
            (abs(top_reaction)+abs(bottom_reaction)),tiny_value)

        write(stress_history_unit,'(I10,5(1X,ES18.10E3))') &
            step, axial_strain_percent, axial_stress_mpa, &
            top_reaction, bottom_reaction, boundary_unbalance

        olddisp = disp
        ! Adaptive dynamic relaxation coefficient
        cn = 0.0d0
        cn1 = 0.0d0
        cn2 = 0.0d0

        !$omp parallel private(i) reduction(+:cn1,cn2)
        !$omp do schedule(static)
        do i = 1,totnode
            ! The x direction is free for both internal and virtual-boundary
            ! nodes, so it participates in the ADR damping estimate.
            if (velhalfold(i,1) .ne. 0.0d0) then
                cn1 = cn1 - disp(i,1) * disp(i,1) * &
                    (pforce(i,1) / massvec(i,1) - pforceold(i,1) / massvec(i,1)) / &
                    (dt * velhalfold(i,1))
            endif
            cn2 = cn2 + disp(i,1) * disp(i,1)

            ! Only internal-node y DOFs are free.  The y displacements of
            ! the upper and lower virtual boundaries are prescribed loading
            ! DOFs and must not enter either cn1 or cn2.
            if (i <= totint) then
                if (velhalfold(i,2) .ne. 0.0d0) then
                    cn1 = cn1 - disp(i,2) * disp(i,2) * &
                        (pforce(i,2) / massvec(i,2) - &
                         pforceold(i,2) / massvec(i,2)) / &
                        (dt * velhalfold(i,2))
                endif
                cn2 = cn2 + disp(i,2) * disp(i,2)
            endif
        enddo
        !$omp end parallel
        
        if (cn2 .ne. 0.0d0) then
            if (cn1 / cn2 > 0.0d0) then
                cn = 2.0d0 / dt * dsqrt(cn1 / cn2)
            else
                cn = 0.0d0
            endif
        else
            cn = 0.0d0
        endif

        if (cn > 2.0d0 / dt) then
            cn = 1.9d0 / dt
        endif
        
        !内部区域，x和y方向自由更新
        do i = 1,totint
            if (step == 1) then
                velhalf(i,1) = dt / massvec(i,1) * &
                    (pforce(i,1) + bforce(i,1)) / 2.0d0
                velhalf(i,2) = dt / massvec(i,2) * &
                    (pforce(i,2) + bforce(i,2)) / 2.0d0
            else
                velhalf(i,1) = ((2.0d0 - cn * dt) * velhalfold(i,1) + &
                    2.0d0 * dt / massvec(i,1) * &
                    (pforce(i,1) + bforce(i,1))) / (2.0d0 + cn * dt)

                velhalf(i,2) = ((2.0d0 - cn * dt) * velhalfold(i,2) + &
                    2.0d0 * dt / massvec(i,2) * &
                    (pforce(i,2) + bforce(i,2))) / (2.0d0 + cn * dt)
            endif

            disp(i,1) = disp(i,1) + velhalf(i,1) * dt
            disp(i,2) = disp(i,2) + velhalf(i,2) * dt

            velhalfold(i,1) = velhalf(i,1)
            velhalfold(i,2) = velhalf(i,2)

            pforceold(i,1) = pforce(i,1)
            pforceold(i,2) = pforce(i,2)
        enddo
        
        !底部边界：x方向自由迭代，y方向保持本步规定的向上位移
        do i = totint + 1,totbottom
            if (step == 1) then
                velhalf(i,1) = dt / massvec(i,1) * &
                    (pforce(i,1) + bforce(i,1)) / 2.0d0
            else
                velhalf(i,1) = ((2.0d0 - cn * dt) * velhalfold(i,1) + &
                    2.0d0 * dt / massvec(i,1) * &
                    (pforce(i,1) + bforce(i,1))) / (2.0d0 + cn * dt)
            endif

            disp(i,1) = disp(i,1) + velhalf(i,1) * dt
            velhalfold(i,1) = velhalf(i,1)
            pforceold(i,1) = pforce(i,1)
        enddo
        
        !顶部边界：x方向自由迭代，y方向保持本步规定的向下位移
        do i = totbottom+1, tottop
            if (step == 1) then
                velhalf(i,1) = dt / massvec(i,1) * &
                    (pforce(i,1) + bforce(i,1)) / 2.0d0
            else
                velhalf(i,1) = ((2.0d0 - cn * dt) * velhalfold(i,1) + &
                    2.0d0 * dt / massvec(i,1) * &
                    (pforce(i,1) + bforce(i,1))) / (2.0d0 + cn * dt)
            endif

            disp(i,1) = disp(i,1) + velhalf(i,1) * dt
            velhalfold(i,1) = velhalf(i,1)
            pforceold(i,1) = pforce(i,1)
        enddo

        ! ================================================================
        ! 输出位移和不可逆拉伸-剪切联合损伤。
        ! ================================================================
        ! With 0.85 mm/500000 steps, retain the original output cadence and
        ! the two additional snapshots used by the earlier post-processing.
        if (mod(step,5000) == 0) then
            write(filename,'("output_",I0,".txt")') step
            open(65, file = filename)
                    do i = 1, totnode
                        write(65,'(7(ES24.16E3,1X))') coord(i,1), coord(i,2), &
                            disp(i,1), disp(i,2), node_damage(i), &
                            damage_tension(i), damage_shear(i)
                    enddo
            close(65)

            ! Keep the original five-column field output unchanged.  Write
            ! damage diagnostics to a separate file with the same step number.
            write(filename,'("diagnostic_",I0,".txt")') step
            open(newunit=diagnostic_unit, file=trim(filename), &
                status='replace', action='write')
            write(diagnostic_unit,'(A)') &
                'x y ux uy Dtotal DI DII FI FII epsIeq gammaIIeq sigma1 sigma3 kexi'
            do i = 1, totnode
                write(diagnostic_unit,'(14(ES24.16E3,1X))') &
                    coord(i,1), coord(i,2), disp(i,1), disp(i,2), &
                    node_damage(i), damage_tension(i), damage_shear(i), &
                    failure_index_I(i), failure_index_II(i), &
                    equivalent_normal_strain_output(i), &
                    equivalent_shear_strain_output(i), &
                    principal_stress_1_output(i), &
                    principal_stress_3_output(i), kexi(i)
            enddo
            close(diagnostic_unit)
        endif
        
    enddo

    close(stress_history_unit)

    write(*,'(A)') 'Simulation completed.'

end program mixed_mode_uniaxial_compression_2d

