module globalParameters
    implicit none

    ! Figure 11(a) single-flaw specimen. With dx=0.4 mm, 191 columns
    ! give a centered 76.0 mm point span (76.2 mm nominal specimen width),
    ! while 381 rows give the nominal 152.4 mm specimen height.
    integer, parameter :: ndivx = 95
    integer, parameter :: ndivy = 190
    integer, parameter :: nbnd = 3
    integer, parameter :: maxfam = 32

    integer, parameter :: nsteps = 700000
    integer, parameter :: ntotmax = ndivx*ndivy + 2*nbnd*ndivx

    real*8, parameter :: pi = 3.1415926535897932384626433832795d0
    real*8, parameter :: plate_width = 76.0d-3
    real*8, parameter :: plate_height = 152.0d-3
    real*8, parameter :: dx = 0.8d-3
    real*8, parameter :: delta = 3.015*dx
    real*8, parameter :: thickness = dx
    real*8, parameter :: particle_volume = dx*dx*thickness
    real*8, parameter :: particle_radius = 0.5d0*dx

    real*8, parameter :: young_modulus = 5.964d9
    real*8, parameter :: poisson_ratio = 0.15d0
    real*8, parameter :: tensile_strength_value = 2.0d6
    real*8, parameter :: cohesion = 6.0d6                       ! 调小参数开裂早，并且在裂纹扩展过程中会有新的起裂发生。
    real*8, parameter :: friction_angle = 28.0d0*pi/180.0d0
    real*8, parameter :: fracture_energy_I = 300.0d0
    real*8, parameter :: fracture_energy_II = 432.0d0           ! 调大参数裂纹延申的慢，并且在裂纹慢慢扩展过程中会有新的起裂发生，如果调小了，裂纹可能一瞬间扩展过去了，就没有新的裂纹起裂了。

    ! User-selected damage calculation mode:
    ! 1 = tensile only, 2 = shear only, 3 = mixed tensile-shear.
    integer, parameter :: damage_mode = 3

    ! Experimental PDDO option. When true, a fully damaged cnode is
    ! excluded from both the moment matrix and displacement-gradient sum
    ! of an otherwise intact center node. Set to false to restore the
    ! original fixed-neighborhood PDDO calculation.
    logical, parameter :: exclude_fully_damaged_pddo_neighbors = .false.
    real*8, parameter :: pddo_full_damage_threshold = 1.0d0 - 1.0d-12

    ! Crack-surface friction is activated by the smoothed post_damage field.
    real*8, parameter :: contact_damage_threshold = 0.05d0     ! 当前节点大于该值时即认定为损伤，存在接触摩擦
    real*8, parameter :: intact_neighbor_damage_threshold = 0.01d0      !当前邻域点大于该值时即认定为不构成当前节点的外法线方向
    real*8, parameter :: full_damage_threshold = 1.0d0 - 1.0d-12
    real*8, parameter :: friction_regularization = 0.001d0*dx
    real*8, parameter :: contact_normal_tolerance = 1.0d-12
    real*8, parameter :: damage_gradient_tolerance = 1.0d-10

    real*8, parameter :: flaw_length = 12.0d-3
    real*8, parameter :: flaw_width = 1.20d-3
    real*8, parameter :: flaw_angle = 45.0d0*pi/180.0d0

    ! Symmetric compression: the prescribed total shortening is split
    ! equally between the lower and upper virtual boundary layers.
    real*8, parameter :: total_compression_displacement = 0.85d-3
    real*8, parameter :: final_boundary_displacement = &
        0.5d0*total_compression_displacement
    real*8, parameter :: load_increment = &
        final_boundary_displacement/dble(nsteps)

    ! Adaptive dynamic relaxation uses pseudo-time and a fictitious mass.
    ! No physical density is used.
    real*8, parameter :: dt = 1.0d0
    real*8, parameter :: tiny_value = 1.0d-18
end module globalParameters


module Global_arrays
    implicit none

    integer, allocatable :: grid_node(:,:), col(:), flag_node(:)
    integer, allocatable :: nodefam(:), fail(:), flag(:,:)
    real*8, allocatable :: coord(:,:), disp(:,:), olddisp(:,:)
    real*8, allocatable :: velhalf(:,:), velhalfold(:,:)
    real*8, allocatable :: pforce(:,:), pforceold(:,:)
    real*8, allocatable :: bforce(:,:), massvec(:,:)
    real*8, allocatable :: bond_idist(:), bond_omega(:)
    real*8, allocatable :: bond_fac(:), scr(:)
    real*8, allocatable :: emod(:), pratio(:), tensile_strength(:)
    real*8, allocatable :: shear_strength(:), G1(:), G2(:), bc(:)
    real*8, allocatable :: bulk_modulus(:), shear_modulus(:)
    real*8, allocatable :: kexi(:), m_volume(:), vol_all(:)
    real*8, allocatable :: node_damage(:), post_damage(:)
    real*8, allocatable :: damage_tension(:), damage_shear(:)
    real*8, allocatable :: failure_index_I(:), failure_index_II(:)
    real*8, allocatable :: damage_initiation_lambda(:)
    real*8, allocatable :: equivalent_normal_strain_output(:)
    real*8, allocatable :: equivalent_shear_strain_output(:)
    real*8, allocatable :: principal_stress_1_output(:)
    real*8, allocatable :: principal_stress_3_output(:)
    real*8, allocatable :: flag_gradient(:,:,:)
    real*8, allocatable :: flag_pddo_kernel(:,:,:)
    real*8, allocatable :: flag_spatial_increment_gradient(:,:,:)
    real*8, allocatable :: flag_deformation_inverse(:,:,:)
    real*8, allocatable :: stress_xx(:), stress_yy(:), stress_xy(:)
    logical, allocatable :: damage_activated(:)
    logical, allocatable :: tension_activated(:), shear_activated(:)
    logical, allocatable :: flag_gradient_valid(:)
end module Global_arrays


module scalars
    implicit none

    integer :: i, j, ix, iy, ii, jj, di, dj, idx, cnode
    integer :: nnum, totnode = 0, totint, totbottom, tottop
    integer :: step, flag_nnum, stress_history_unit, diagnostic_unit
    integer :: nearfield_count

    real*8 :: vol = 0.0d0
    real*8 :: coordx, coordy, along_flaw, across_flaw
    real*8 :: along_i, along_j, across_i, across_j
    real*8 :: crack_cross_parameter, crack_cross_along
    real*8 :: temp, idist, omega, fac, radij, thick, cfri, nlength
    real*8 :: epsilon_xx, epsilon_yy, epsilon_xy
    real*8 :: mean_strain, strain_radius, mean_stress, stress_radius
    real*8 :: principal_strain_1, principal_strain_3
    real*8 :: principal_stress_1, principal_stress_3
    real*8 :: equivalent_normal_strain, equivalent_shear_strain
    real*8 :: normal_initial_strain, shear_initial_strain
    real*8 :: normal_failure_strain, shear_failure_strain
    real*8 :: normalized_normal_strain, normalized_shear_strain
    real*8 :: normalized_normal_initial, normalized_shear_initial
    real*8 :: normal_mode_measure, shear_mode_measure
    real*8 :: lambda, lambda_p, cos_theta, sin_theta
    real*8 :: damage_lambda, damage_lambda_p, mode_weight_sum
    real*8 :: rankine_index, mohr_coulomb_index, mc_numerator
    real*8 :: damage_trial, damage_old, damage_new, damage_increment
    real*8 :: damage_sum
    real*8 :: damage_i, damage_j, damage_tension_i, damage_tension_j
    real*8 :: tV_i, tS_i, tF_i, tV_j, tS_j, tF_j
    real*8 :: cn, cn1, cn2
    real*8 :: bottom_reaction, top_reaction, section_area
    real*8 :: axial_stress_mpa, axial_strain_percent
    real*8 :: boundary_unbalance
    real*8 :: Y_vec(2), t_i(2), t_j(2), dforce(2)
end module scalars
