module displacement_gradient_module
    use globalParameters, only: dx, maxfam, &
        exclude_fully_damaged_pddo_neighbors, &
        pddo_full_damage_threshold
    use Global_arrays, only: coord, disp, olddisp, nodefam, fail, bond_omega, &
        bond_fac, flag, node_damage, bulk_modulus, shear_modulus
    use scalars, only: totnode, vol
    implicit none
    private

    !理论上，M 维、最高总阶数 N 的基函数数量为组合数 C(M+N,N)。
    integer, parameter :: n_pddo = 6
    integer, parameter :: bicgstab_max_iterations = 100
    real*8, parameter :: bicgstab_relative_tolerance = 1.0d-12

    public :: build_flag_node_list_from_flag
    public :: calculate_selected_node_displacement_gradients

    contains

    ! 直接读取主程序 flag(:,1)
    subroutine build_flag_node_list_from_flag(flag_node)
        implicit none

        integer, allocatable, intent(out) :: flag_node(:)

        integer :: i
        integer :: totnode_of_flag
        integer :: flag_nnum

        call check_global_gradient_arrays()

        totnode_of_flag = count(flag(1:totnode,1) == 1)
        allocate(flag_node(totnode_of_flag))

        flag_nnum = 0
        do i = 1,totnode
            if (flag(i,1) == 1) then
                flag_nnum = flag_nnum + 1
                flag_node(flag_nnum) = i
            endif
        enddo
    end subroutine build_flag_node_list_from_flag

    ! Calculate N=2 PDDO displacement gradients only at selected global
    ! node ids. This avoids looping over every node at every output step.
    !
    ! displacement_gradient(k,component,direction):
    !   (:,1,1) = du_x/dx, (:,1,2) = du_x/dy
    !   (:,2,1) = du_y/dx, (:,2,2) = du_y/dy
    !
    ! When exclude_fully_damaged_pddo_neighbors is enabled, the same
    ! active-neighbor set must be used to construct A and to accumulate
    ! the displacement gradient.
    subroutine calculate_selected_node_displacement_gradients(flag_node, &
        displacement_gradient, gradient_valid, pddo_kernel, &
        spatial_increment_gradient, deformation_inverse)
        implicit none

        integer, intent(in) :: flag_node(:)
        real*8, intent(out) :: displacement_gradient(:,:,:)
        logical, intent(out) :: gradient_valid(:)
        real*8, intent(out) :: pddo_kernel(:,:,:)
        real*8, intent(out) :: spatial_increment_gradient(:,:,:)
        real*8, intent(out) :: deformation_inverse(:,:,:)

        real*8 :: A_matrix(n_pddo,n_pddo)
        real*8 :: b_array(n_pddo,2)
        real*8 :: aa_matrix(n_pddo,2)
        real*8 :: monomial(n_pddo)
        real*8 :: relative_displacement(2)
        real*8 :: relative_increment(2)
        real*8 :: increment_gradient(2,2)
        real*8 :: deformation(2,2), inverse_local(2,2)
        real*8 :: determinant
        real*8 :: relative_x
        real*8 :: relative_y
        real*8 :: g_10
        real*8 :: g_01
        integer :: flag_nnum
        integer :: totnode_of_flag
        integer :: i
        integer :: j
        integer :: idx
        integer :: cnode
        integer :: row_index
        integer :: column_index
        logical :: matrix_is_nonsingular

        call check_global_gradient_arrays()
        totnode_of_flag = size(flag_node)

        if (size(displacement_gradient,1) < totnode_of_flag .or. &
            size(displacement_gradient,2) < 2 .or. &
            size(displacement_gradient,3) < 2) then
            error stop 'PDDO gradient: displacement_gradient has an invalid shape'
        endif
        if (size(gradient_valid) < totnode_of_flag) then
            error stop 'PDDO gradient: gradient_valid is too small'
        endif
        if (size(pddo_kernel,1) < totnode_of_flag .or. &
            size(pddo_kernel,2) < maxfam .or. size(pddo_kernel,3) < 2) then
            error stop 'PDDO gradient: pddo_kernel has an invalid shape'
        endif
        if (size(spatial_increment_gradient,1) < totnode_of_flag .or. &
            size(spatial_increment_gradient,2) < 2 .or. &
            size(spatial_increment_gradient,3) < 2) then
            error stop 'PDDO gradient: spatial increment gradient has an invalid shape'
        endif
        if (size(deformation_inverse,1) < totnode_of_flag .or. &
            size(deformation_inverse,2) < 2 .or. &
            size(deformation_inverse,3) < 2) then
            error stop 'PDDO gradient: deformation inverse has an invalid shape'
        endif

        displacement_gradient = 0.0d0
        gradient_valid = .false.
        pddo_kernel = 0.0d0
        spatial_increment_gradient = 0.0d0
        deformation_inverse = 0.0d0

        ! Use dimensionless bond coordinates xi/dx to keep the moment
        ! matrix well-conditioned. The 1/dx factors restore derivatives
        ! with respect to the physical coordinates.
        ! Basis order: [1, q_x, q_y, q_x^2, q_x*q_y, q_y^2].
        b_array = 0.0d0
        b_array(2,1) = 1.0d0 / dx
        b_array(3,2) = 1.0d0 / dx

        ! Each selected center node owns a separate moment matrix, solver
        ! state, and output slice. Parallelizing this outer loop therefore
        ! requires no reduction or atomic update.
        !$omp parallel do default(shared) schedule(static) &
        !$omp& private(flag_nnum,i,j,idx,cnode,row_index,column_index, &
        !$omp& A_matrix,aa_matrix,monomial,relative_displacement, &
        !$omp& relative_increment,increment_gradient,deformation, &
        !$omp& inverse_local,determinant,relative_x,relative_y,g_10,g_01, &
        !$omp& matrix_is_nonsingular)
        do flag_nnum = 1,totnode_of_flag
            i = flag_node(flag_nnum)

            if (i < 1 .or. i > totnode) cycle

            A_matrix = 0.0d0
            increment_gradient = 0.0d0

            ! The MATLAB reference includes the center point in its family.
            ! At xi = 0 only the constant-constant moment receives a term.
            A_matrix(1,1) = vol

            ! 该循环是计算矩阵A的数值
            do j = 1,maxfam
                idx = (i-1)*maxfam + j
                cnode = nodefam(idx)

                if (cnode == 0) exit
                if (fail(idx) == 0) cycle
                if (exclude_fully_damaged_pddo_neighbors .and. &
                    node_damage(cnode) >= pddo_full_damage_threshold) cycle

                relative_x = (coord(cnode,1) - coord(i,1)) / dx
                relative_y = (coord(cnode,2) - coord(i,2)) / dx

                monomial(1) = 1.0d0
                monomial(2) = relative_x
                monomial(3) = relative_y
                monomial(4) = relative_x * relative_x
                monomial(5) = relative_x * relative_y
                monomial(6) = relative_y * relative_y

                do row_index = 1,n_pddo
                    do column_index = 1,n_pddo
                        A_matrix(row_index,column_index) = &
                            A_matrix(row_index,column_index) + &
                            bond_omega(idx) * monomial(row_index) * &
                            monomial(column_index) * vol * bond_fac(idx)
                    enddo
                enddo
            enddo

            ! 矩阵求解输出x方向和y方向导数核函数
            call solve_6x6_bicgstab_two_rhs(A_matrix, b_array, &
                aa_matrix, matrix_is_nonsingular)

            if (.not. matrix_is_nonsingular) cycle

            do j = 1,maxfam
                idx = (i-1)*maxfam + j
                cnode = nodefam(idx)

                if (cnode == 0) exit
                if (fail(idx) == 0) cycle
                if (exclude_fully_damaged_pddo_neighbors .and. &
                    node_damage(cnode) >= pddo_full_damage_threshold) cycle

                relative_x = (coord(cnode,1) - coord(i,1)) / dx
                relative_y = (coord(cnode,2) - coord(i,2)) / dx

                monomial(1) = 1.0d0
                monomial(2) = relative_x
                monomial(3) = relative_y
                monomial(4) = relative_x * relative_x
                monomial(5) = relative_x * relative_y
                monomial(6) = relative_y * relative_y

                g_10 = bond_omega(idx) * &
                    dot_product(monomial, aa_matrix(:,1))
                g_01 = bond_omega(idx) * &
                    dot_product(monomial, aa_matrix(:,2))
                pddo_kernel(flag_nnum,j,1) = g_10
                pddo_kernel(flag_nnum,j,2) = g_01

                relative_displacement(:) = disp(cnode,:) - disp(i,:)
                relative_increment(:) = &
                    (disp(cnode,:)-olddisp(cnode,:)) - &
                    (disp(i,:)-olddisp(i,:))

                displacement_gradient(flag_nnum,1,1) = &
                    displacement_gradient(flag_nnum,1,1) + &
                    relative_displacement(1) * g_10 * vol * bond_fac(idx)
                displacement_gradient(flag_nnum,1,2) = &
                    displacement_gradient(flag_nnum,1,2) + &
                    relative_displacement(1) * g_01 * vol * bond_fac(idx)
                displacement_gradient(flag_nnum,2,1) = &
                    displacement_gradient(flag_nnum,2,1) + &
                    relative_displacement(2) * g_10 * vol * bond_fac(idx)
                displacement_gradient(flag_nnum,2,2) = &
                    displacement_gradient(flag_nnum,2,2) + &
                    relative_displacement(2) * g_01 * vol * bond_fac(idx)

                increment_gradient(:,1) = increment_gradient(:,1) + &
                    relative_increment*g_10*vol*bond_fac(idx)
                increment_gradient(:,2) = increment_gradient(:,2) + &
                    relative_increment*g_01*vol*bond_fac(idx)
            enddo

            deformation = displacement_gradient(flag_nnum,:,:)
            deformation(1,1) = deformation(1,1) + 1.0d0
            deformation(2,2) = deformation(2,2) + 1.0d0
            determinant = deformation(1,1)*deformation(2,2) - &
                deformation(1,2)*deformation(2,1)
            if (abs(determinant) <= 1.0d-12) cycle

            inverse_local(1,1) = deformation(2,2)/determinant
            inverse_local(1,2) = -deformation(1,2)/determinant
            inverse_local(2,1) = -deformation(2,1)/determinant
            inverse_local(2,2) = deformation(1,1)/determinant
            deformation_inverse(flag_nnum,:,:) = inverse_local
            spatial_increment_gradient(flag_nnum,:,:) = &
                matmul(increment_gradient,inverse_local)
            gradient_valid(flag_nnum) = .true.
        enddo
        !$omp end parallel do
    end subroutine calculate_selected_node_displacement_gradients


    subroutine check_global_gradient_arrays()
        implicit none

        if (totnode < 1) error stop 'PDDO gradient: totnode must be positive'
        if (.not. allocated(coord)) error stop 'PDDO gradient: coord is not allocated'
        if (.not. allocated(disp)) error stop 'PDDO gradient: disp is not allocated'
        if (.not. allocated(olddisp)) error stop 'PDDO gradient: olddisp is not allocated'
        if (.not. allocated(nodefam)) error stop 'PDDO gradient: nodefam is not allocated'
        if (.not. allocated(fail)) error stop 'PDDO gradient: fail is not allocated'
        if (.not. allocated(bond_omega)) error stop 'PDDO gradient: bond_omega is not allocated'
        if (.not. allocated(bond_fac)) error stop 'PDDO gradient: bond_fac is not allocated'
        if (.not. allocated(flag)) error stop 'PDDO gradient: flag is not allocated'
        if (.not. allocated(node_damage)) error stop 'PDDO gradient: node_damage is not allocated'
        if (size(flag,1) < totnode .or. size(flag,2) < 1) then
            error stop 'PDDO gradient: flag has an invalid shape'
        endif
        if (size(node_damage) < totnode) then
            error stop 'PDDO gradient: node_damage is too small'
        endif
    end subroutine check_global_gradient_arrays


    ! Solve the fixed 6-by-6 PDDO systems with the BiCGSTAB iteration
    ! presented in Appendix B of the supplied PDF. The two columns of B
    ! are solved independently because they represent D^[1,0] and D^[0,1].
    subroutine solve_6x6_bicgstab_two_rhs(matrix_a, matrix_b, solution, success)
        implicit none

        real*8, intent(in) :: matrix_a(n_pddo,n_pddo)
        real*8, intent(in) :: matrix_b(n_pddo,2)
        real*8, intent(out) :: solution(n_pddo,2)
        logical, intent(out) :: success

        real*8 :: residual(n_pddo)
        real*8 :: shadow_residual(n_pddo)
        real*8 :: search_direction(n_pddo)
        real*8 :: matrix_search(n_pddo)
        real*8 :: intermediate_residual(n_pddo)
        real*8 :: matrix_intermediate(n_pddo)
        real*8 :: right_vector(n_pddo)
        real*8 :: rho_previous
        real*8 :: rho_current
        real*8 :: alpha
        real*8 :: beta
        real*8 :: omega_previous
        real*8 :: omega_current
        real*8 :: denominator
        real*8 :: right_norm
        real*8 :: residual_norm
        real*8 :: breakdown_scale
        integer :: rhs_index
        integer :: iteration
        logical :: rhs_converged

        solution = 0.0d0
        success = .true.

        do rhs_index = 1,2
            right_vector = matrix_b(:,rhs_index)
            right_norm = sqrt(dot_product(right_vector,right_vector))

            if (right_norm <= tiny(1.0d0)) cycle

            ! PDF steps (1)-(5): x0=0, r0=b-A*x0, r_tilde=r0,
            ! rho0=alpha=omega0=1, and v0=p0=0.
            solution(:,rhs_index) = 0.0d0
            residual = right_vector
            shadow_residual = residual
            search_direction = 0.0d0
            matrix_search = 0.0d0
            rho_previous = 1.0d0
            alpha = 1.0d0
            omega_previous = 1.0d0
            rhs_converged = .false.

            do iteration = 1,bicgstab_max_iterations
                ! PDF steps (6)-(10).
                rho_current = dot_product(shadow_residual,residual)
                breakdown_scale = epsilon(1.0d0) * &
                    sqrt(dot_product(shadow_residual,shadow_residual)) * &
                    sqrt(dot_product(residual,residual))
                if (abs(rho_current) <= max(tiny(1.0d0),breakdown_scale)) exit

                beta = (rho_current/rho_previous) * (alpha/omega_previous)
                search_direction = residual + beta * &
                    (search_direction - omega_previous*matrix_search)
                matrix_search = matmul(matrix_a,search_direction)

                denominator = dot_product(shadow_residual,matrix_search)
                breakdown_scale = epsilon(1.0d0) * &
                    sqrt(dot_product(shadow_residual,shadow_residual)) * &
                    sqrt(dot_product(matrix_search,matrix_search))
                if (abs(denominator) <= max(tiny(1.0d0),breakdown_scale)) exit

                alpha = rho_current / denominator

                ! PDF step (11). When s is already converged, the omega
                ! correction is unnecessary and t may be the zero vector.
                intermediate_residual = residual - alpha*matrix_search
                residual_norm = sqrt(dot_product(intermediate_residual, &
                    intermediate_residual))
                if (residual_norm/right_norm <= bicgstab_relative_tolerance) then
                    solution(:,rhs_index) = solution(:,rhs_index) + &
                        alpha*search_direction
                    rhs_converged = .true.
                    exit
                endif

                ! PDF steps (12)-(15).
                matrix_intermediate = matmul(matrix_a,intermediate_residual)
                denominator = dot_product(matrix_intermediate,matrix_intermediate)
                if (denominator <= tiny(1.0d0)) exit

                omega_current = dot_product(matrix_intermediate, &
                    intermediate_residual) / denominator
                if (abs(omega_current) <= epsilon(1.0d0)) exit

                solution(:,rhs_index) = solution(:,rhs_index) + &
                    alpha*search_direction + omega_current*intermediate_residual
                residual = intermediate_residual - omega_current*matrix_intermediate

                ! PDF step (16): relative residual stopping criterion.
                residual_norm = sqrt(dot_product(residual,residual))
                if (residual_norm/right_norm <= bicgstab_relative_tolerance) then
                    rhs_converged = .true.
                    exit
                endif

                rho_previous = rho_current
                omega_previous = omega_current
            enddo

            if (.not. rhs_converged) then
                solution = 0.0d0
                success = .false.
                return
            endif
        enddo
    end subroutine solve_6x6_bicgstab_two_rhs

end module displacement_gradient_module
