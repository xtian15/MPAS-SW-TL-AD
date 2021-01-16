module mpas_sw_module
  implicit none
  
  integer, parameter :: RKIND  = selected_real_kind(12)
  integer, parameter :: StrKIND = 512

  ! ========== constants ==========
  real (kind=RKIND), parameter :: pii     = 3.141592653589793   !< Constant: Pi
  real (kind=RKIND), parameter :: a       = 6371229.0           !< Constant: Spherical Earth radius [m]
  real (kind=RKIND), parameter :: omega   = 7.29212e-5          !< Constant: Angular rotation rate of the Earth [s-1]
  real (kind=RKIND), parameter :: gravity = 9.80616             !< Constant: Acceleration due to gravity [m s-2]
  real (kind=RKIND), parameter :: rgas    = 287.0               !< Constant: Gas constant for dry air [J kg-1 K-1]
  real (kind=RKIND), parameter :: rv      = 461.6               !< Constant: Gas constant for water vapor [J kg-1 K-1]
  real (kind=RKIND), parameter :: rvord   = rv/rgas             !
  real (kind=RKIND), parameter :: cp      = 7.*rgas/2.          !< Constant: Specific heat of dry air at constant pressure [J kg-1 K-1]
  real (kind=RKIND), parameter :: cv      = cp - rgas           !< Constant: Specific heat of dry air at constant volume [J kg-1 K-1]
  real (kind=RKIND), parameter :: cvpm    = -cv/cp              !
  real (kind=RKIND), parameter :: prandtl = 1.0                 !< Constant: Prandtl number
  
  ! ========== end constants ==========
  
  ! ========== configurations ==========
  integer :: config_test_case = 5
  character(len=StrKIND) :: config_time_integration = 'RK4'
  real(kind=RKIND) :: config_dt = 150.0
  character(len=StrKIND) :: config_start_time = '2000-01-01_00:00:00'
  character(len=StrKIND) :: config_stop_time = 'none'
  character(len=StrKIND) :: config_run_duration = 'none'
  integer :: config_stats_interval = 100
  logical :: config_h_ScaleWithMesh = .false.
  real(kind=RKIND) :: config_h_mom_eddy_visc2 = 0.0
  real(kind=RKIND) :: config_h_mom_eddy_visc4 = 0.0
  real(kind=RKIND) :: config_h_tracer_eddy_diff2 = 0.0
  real(kind=RKIND) :: config_h_tracer_eddy_diff4 = 0.0
  integer :: config_thickness_adv_order = 2
  integer :: config_tracer_adv_order = 2
  logical :: config_positive_definite = .false.
  logical :: config_monotonic = .false.
  logical :: config_wind_stress = .false.
  logical :: config_bottom_drag = .false.
  real(kind=RKIND) :: config_apvm_upwinding = 0.5
  integer :: config_num_halos = 2
  
  ! ========== dimensions ==========
  integer :: nVertLevels=1
  integer :: nCells, nEdges, nVertices, nTracers
  integer :: cellStart, cellEnd
  integer :: vertexStart, vertexEnd
  integer :: edgeStart, edgeEnd
  integer :: vertexDegree
  integer :: nCellsSolve, nEdgesSolve
  integer :: maxEdges, maxEdges2
  integer :: cellSolveStart, cellSolveEnd
  integer :: edgeSolveStart, edgeSolveEnd
  integer :: R3

  ! ========== state ==========
  real(kind=RKIND), dimension(:,:,:), allocatable :: u
  real(kind=RKIND), dimension(:,:,:), allocatable :: h
  real(kind=RKIND), dimension(:,:,:,:), allocatable :: tracers
  
  ! ========== state TLM ==========
  real(kind=RKIND), dimension(:,:,:), allocatable :: u_tl
  real(kind=RKIND), dimension(:,:,:), allocatable :: h_tl
  real(kind=RKIND), dimension(:,:,:,:), allocatable :: tracers_tl

  ! ========== state ADJ ==========
  real(kind=RKIND), dimension(:,:,:), allocatable :: u_ad
  real(kind=RKIND), dimension(:,:,:), allocatable :: h_ad
  real(kind=RKIND), dimension(:,:,:,:), allocatable :: tracers_ad

  ! ========== top level ==========
  logical :: on_a_sphere
  logical :: is_periodic
  real(kind=RKIND) :: x_period, y_period

  ! ========== mesh ==========
  real(kind=RKIND), dimension(:,:), allocatable :: u_src
  real(kind=RKIND), dimension(:), allocatable :: h_s
  real(kind=RKIND), dimension(:), allocatable :: fCell
  real(kind=RKIND), dimension(:), allocatable :: fVertex
  real(kind=RKIND), dimension(:), allocatable :: fEdge
  real(kind=RKIND), dimension(:), allocatable :: areaCell
  real(kind=RKIND), dimension(:), allocatable :: invAreaCell
  real(kind=RKIND), dimension(:), allocatable :: dvEdge
  real(kind=RKIND), dimension(:), allocatable :: dcEdge
  real(kind=RKIND), dimension(:,:), allocatable :: weightsOnEdge
  real(kind=RKIND), dimension(:,:), allocatable :: kiteAreasOnVertex
  real(kind=RKIND), dimension(:), allocatable :: meshDensity
  real(kind=RKIND), dimension(:), allocatable :: areaTriangle
  integer, dimension(:,:), allocatable :: boundaryCell
  integer, dimension(:,:), allocatable :: boundaryEdge
  integer, dimension(:,:), allocatable :: advCells
  integer, dimension(:,:), allocatable :: cellsOnCell
  integer, dimension(:,:), allocatable :: cellsOnEdge
  integer, dimension(:,:), allocatable :: cellsOnVertex
  integer, dimension(:,:), allocatable :: edgesOnCell
  integer, dimension(:,:), allocatable :: edgesOnEdge
  integer, dimension(:,:), allocatable :: edgesOnVertex
  integer, dimension(:,:), allocatable :: verticesOnCell
  integer, dimension(:,:), allocatable :: verticesOnEdge
  integer, dimension(:), allocatable :: nEdgesOnCell
  integer, dimension(:), allocatable :: nEdgesOnEdge
  real(kind=RKIND), dimension(:), allocatable :: latCell
  real(kind=RKIND), dimension(:), allocatable :: lonCell
  real(kind=RKIND), dimension(:), allocatable :: xCell, yCell, zCell
  real(kind=RKIND), dimension(:), allocatable :: meshScalingDel2
  real(kind=RKIND), dimension(:), allocatable :: meshScalingDel4
  real(kind=RKIND), dimension(:,:,:), allocatable :: deriv_two
  real(kind=RKIND), dimension(:), allocatable :: latEdge
  real(kind=RKIND), dimension(:), allocatable :: lonEdge
  real(kind=RKIND), dimension(:), allocatable :: xEdge, yEdge, zEdge
  real(kind=RKIND), dimension(:), allocatable :: angleEdge
  real(kind=RKIND), dimension(:), allocatable :: latVertex
  real(kind=RKIND), dimension(:), allocatable :: lonVertex
  real(kind=RKIND), dimension(:), allocatable :: xVertex, yVertex, zVertex
  real(kind=RKIND) :: sphere_radius
  real(kind=RKIND), dimension(:,:), allocatable :: defc_a
  real(kind=RKIND), dimension(:,:), allocatable :: defc_b
  real(kind=RKIND), dimension(:,:), allocatable :: localVerticalUnitVectors
  real(kind=RKIND), dimension(:,:), allocatable :: edgeNormalVectors
  real(kind=RKIND), dimension(:,:,:), allocatable :: cellTangentPlane
  real(kind=RKIND), dimension(:,:,:), allocatable :: coeffs_reconstruct

  ! ========== diag ==========
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructX
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructY
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructZ
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructZonal
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructMeridional
  real(kind=RKIND), dimension(:,:), allocatable :: v
  real(kind=RKIND), dimension(:,:), allocatable :: h_edge
  real(kind=RKIND), dimension(:,:), allocatable :: h_vertex
  real(kind=RKIND), dimension(:,:), allocatable :: pv_edge
  real(kind=RKIND), dimension(:,:), allocatable :: pv_vertex
  real(kind=RKIND), dimension(:,:), allocatable :: pv_cell
  real(kind=RKIND), dimension(:,:), allocatable :: circulation
  real(kind=RKIND), dimension(:,:), allocatable :: vorticity
  real(kind=RKIND), dimension(:,:), allocatable :: divergence
  real(kind=RKIND), dimension(:,:), allocatable :: ke
  real(kind=RKIND), dimension(:,:), allocatable :: vh
  real(kind=RKIND), dimension(:,:), allocatable :: vorticity_cell
  real(kind=RKIND), dimension(:,:), allocatable :: gradPVn
  real(kind=RKIND), dimension(:,:), allocatable :: gradPVt

  ! ========== diag TLM ==========
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructX_tl
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructY_tl
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructZ_tl
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructZonal_tl
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructMeridional_tl
  real(kind=RKIND), dimension(:,:), allocatable :: v_tl
  real(kind=RKIND), dimension(:,:), allocatable :: h_edge_tl
  real(kind=RKIND), dimension(:,:), allocatable :: h_vertex_tl
  real(kind=RKIND), dimension(:,:), allocatable :: pv_edge_tl
  real(kind=RKIND), dimension(:,:), allocatable :: pv_vertex_tl
  real(kind=RKIND), dimension(:,:), allocatable :: pv_cell_tl
  real(kind=RKIND), dimension(:,:), allocatable :: circulation_tl
  real(kind=RKIND), dimension(:,:), allocatable :: vorticity_tl
  real(kind=RKIND), dimension(:,:), allocatable :: divergence_tl
  real(kind=RKIND), dimension(:,:), allocatable :: ke_tl
  real(kind=RKIND), dimension(:,:), allocatable :: vh_tl
  real(kind=RKIND), dimension(:,:), allocatable :: vorticity_cell_tl
  real(kind=RKIND), dimension(:,:), allocatable :: gradPVn_tl
  real(kind=RKIND), dimension(:,:), allocatable :: gradPVt_tl

  ! ========== diag ADJ ==========
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructX_ad
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructY_ad
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructZ_ad
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructZonal_ad
  real(kind=RKIND), dimension(:,:), allocatable :: uReconstructMeridional_ad
  real(kind=RKIND), dimension(:,:), allocatable :: v_ad
  real(kind=RKIND), dimension(:,:), allocatable :: h_edge_ad
  real(kind=RKIND), dimension(:,:), allocatable :: h_vertex_ad
  real(kind=RKIND), dimension(:,:), allocatable :: pv_edge_ad
  real(kind=RKIND), dimension(:,:), allocatable :: pv_vertex_ad
  real(kind=RKIND), dimension(:,:), allocatable :: pv_cell_ad
  real(kind=RKIND), dimension(:,:), allocatable :: circulation_ad
  real(kind=RKIND), dimension(:,:), allocatable :: vorticity_ad
  real(kind=RKIND), dimension(:,:), allocatable :: divergence_ad
  real(kind=RKIND), dimension(:,:), allocatable :: ke_ad
  real(kind=RKIND), dimension(:,:), allocatable :: vh_ad
  real(kind=RKIND), dimension(:,:), allocatable :: vorticity_cell_ad
  real(kind=RKIND), dimension(:,:), allocatable :: gradPVn_ad
  real(kind=RKIND), dimension(:,:), allocatable :: gradPVt_ad

  ! ========== tend ==========
  real(kind=RKIND), dimension(:,:), allocatable :: tend_u
  real(kind=RKIND), dimension(:,:), allocatable :: tend_h
  real(kind=RKIND), dimension(:,:,:), allocatable :: tend_tracer

  ! ========== tend TLM ==========
  real(kind=RKIND), dimension(:,:), allocatable :: tend_u_tl
  real(kind=RKIND), dimension(:,:), allocatable :: tend_h_tl
  real(kind=RKIND), dimension(:,:,:), allocatable :: tend_tracer_tl
  
  ! ========== tend ADJ ==========
  real(kind=RKIND), dimension(:,:), allocatable :: tend_u_ad
  real(kind=RKIND), dimension(:,:), allocatable :: tend_h_ad
  real(kind=RKIND), dimension(:,:,:), allocatable :: tend_tracer_ad
  
  ! ========== end FWD storing ==========

  ! ===== transient indices =====
  integer :: rk_step
  real(kind=RKIND) :: LHS, RHS, pt

CONTAINS
  
  subroutine var_allocation()
    implicit none

    ! ===== allocating state =====
    allocate(u(2,nVertLevels,nEdges))
    u=0
    allocate(h(2,nVertLevels,nCells))
    h=0
    allocate(tracers(2,nTracers,nVertLevels,nCells))
    tracers=0

    ! ===== allocating mesh =====
    allocate(u_src(nVertLevels, nEdges))
    u_src=0
    allocate(h_s(nCells))
    h_s=0
    allocate(fCell(nCells))
    fCell=0
    allocate(fVertex(nVertices))
    fVertex=0
    allocate(fEdge(nEdges))
    fEdge=0
    allocate(areaTriangle(nVertices))
    areaTriangle=0
    allocate(areaCell(nCells))
    areaCell=0
    allocate(dvEdge(nEdges))
    dvEdge=0
    allocate(dcEdge(nEdges))
    dcEdge=0
    allocate(weightsOnEdge(maxEdges2,nEdges))
    weightsOnEdge=0
    allocate(kiteAreasOnVertex(vertexDegree,nVertices))
    kiteAreasOnVertex=0
    allocate(meshDensity(nCells))
    meshDensity=0
    allocate(boundaryCell(nVertLevels,nCells))
    boundaryCell=0
    allocate(boundaryEdge(nVertLevels,nEdges))
    boundaryEdge=0
    allocate(advCells(21, nCells))
    advCells=0
    allocate(cellsOnCell(maxEdges,nCells))
    cellsOnCell=0
    allocate(cellsOnEdge(2,nEdges))
    cellsOnEdge=0
    allocate(cellsOnVertex(vertexDegree,nVertices))
    cellsOnVertex=0
    allocate(edgesOnCell(maxEdges,nCells))
    edgesOnCell=0
    allocate(edgesOnEdge(maxEdges2,nEdges))
    edgesOnEdge=0
    allocate(edgesOnVertex(vertexDegree,nVertices))
    edgesOnVertex=0
    allocate(verticesOnCell(maxEdges,nCells))
    verticesOnCell=0
    allocate(verticesOnEdge(2,nEdges))
    verticesOnEdge=0
    allocate(nEdgesOnCell(nCells))
    nEdgesOnCell=0
    allocate(nEdgesOnEdge(nEdges))
    nEdgesOnEdge=0
    allocate(latCell(nCells))
    latCell=0
    allocate(lonCell(nCells))
    lonCell=0
    allocate(xCell(nCells), yCell(nCells), zCell(nCells))
    xCell=0
    yCell=0
    zCell=0
    allocate(meshScalingDel2(nEdges))
    meshScalingDel2=0
    allocate(meshScalingDel4(nEdges))
    meshScalingDel4=0
    allocate(deriv_two(15,2,nEdges))
    deriv_two=0
    allocate(latEdge(nEdges))
    latEdge=0
    allocate(lonEdge(nEdges))
    lonEdge=0
    allocate(xEdge(nEdges), yEdge(nEdges), zEdge(nEdges))
    xEdge=0
    yEdge=0
    zEdge=0
    allocate(angleEdge(nEdges))
    angleEdge=0
    allocate(latVertex(nVertices))
    latVertex=0
    allocate(lonVertex(nVertices))
    lonVertex=0
    allocate(xVertex(nVertices), yVertex(nVertices), zVertex(nVertices))
    xVertex=0
    yVertex=0
    zVertex=0
    allocate(defc_a(maxEdges,nCells))
    defc_a=0
    allocate(defc_b(maxEdges,nCells))
    defc_b=0
    allocate(localVerticalUnitVectors(R3,nCells))
    localVerticalUnitVectors=0
    allocate(edgeNormalVectors(R3,nEdges))
    edgeNormalVectors=0
    allocate(cellTangentPlane(R3,2,nCells))
    cellTangentPlane=0
    allocate(coeffs_reconstruct(R3,maxEdges,nCells))
    coeffs_reconstruct=0

    ! ===== allocating diag =====
    allocate(uReconstructX(nVertLevels,nCells))
    allocate(uReconstructY(nVertLevels,nCells))
    allocate(uReconstructZ(nVertLevels,nCells))
    allocate(uReconstructZonal(nVertLevels,nCells))
    allocate(uReconstructMeridional(nVertLevels,nCells))
    allocate(v(nVertLevels,nEdges))
    allocate(h_edge(nVertLevels,nEdges))
    allocate(h_vertex(nVertLevels,nVertices))
    allocate(pv_edge(nVertLevels,nEdges))
    allocate(pv_vertex(nVertLevels,nVertices))
    allocate(pv_cell(nVertLevels,nCells))
    allocate(circulation(nVertLevels,nVertices))
    allocate(vorticity(nVertLevels,nVertices))
    allocate(divergence(nVertLevels,nCells))
    allocate(ke(nVertLevels,nCells))
    allocate(vh(nVertLevels,nEdges))
    allocate(vorticity_cell(nVertLevels,nCells))
    allocate(gradPVn(nVertLevels,nEdges))
    allocate(gradPVt(nVertLevels,nEdges))

    uReconstructX=0
    uReconstructY=0
    uReconstructZ=0
    uReconstructZonal=0
    uReconstructMeridional=0
    v=0
    h_edge=0
    h_vertex=0
    pv_edge=0
    pv_vertex=0
    pv_cell=0
    circulation=0
    vorticity=0
    divergence=0
    ke=0
    vh=0
    vorticity_cell=0
    gradPVn=0
    gradPVt=0

    ! ===== allocating tend =====
    allocate(tend_u(nVertLevels,nEdges))
    allocate(tend_h(nVertLevels,nCells))
    allocate(tend_tracer(nTracers,nVertLevels,nCells))
    tend_u=0
    tend_h=0
    tend_tracer=0

  end subroutine var_allocation  

  subroutine var_allocation_tlm()
    implicit none

    ! ===== allocating state =====
    allocate(u_tl(2,nVertLevels,nEdges))
    u_tl=0
    allocate(h_tl(2,nVertLevels,nCells))
    h_tl=0
    allocate(tracers_tl(2,nTracers,nVertLevels,nCells))
    tracers_tl=0

    ! ===== allocating diag =====
    allocate(uReconstructX_tl(nVertLevels,nCells))
    allocate(uReconstructY_tl(nVertLevels,nCells))
    allocate(uReconstructZ_tl(nVertLevels,nCells))
    allocate(uReconstructZonal_tl(nVertLevels,nCells))
    allocate(uReconstructMeridional_tl(nVertLevels,nCells))
    allocate(v_tl(nVertLevels,nEdges))
    allocate(h_edge_tl(nVertLevels,nEdges))
    allocate(h_vertex_tl(nVertLevels,nVertices))
    allocate(pv_edge_tl(nVertLevels,nEdges))
    allocate(pv_vertex_tl(nVertLevels,nVertices))
    allocate(pv_cell_tl(nVertLevels,nCells))
    allocate(circulation_tl(nVertLevels,nVertices))
    allocate(vorticity_tl(nVertLevels,nVertices))
    allocate(divergence_tl(nVertLevels,nCells))
    allocate(ke_tl(nVertLevels,nCells))
    allocate(vh_tl(nVertLevels,nEdges))
    allocate(vorticity_cell_tl(nVertLevels,nCells))
    allocate(gradPVn_tl(nVertLevels,nEdges))
    allocate(gradPVt_tl(nVertLevels,nEdges))

    uReconstructX_tl=0
    uReconstructY_tl=0
    uReconstructZ_tl=0
    uReconstructZonal_tl=0
    uReconstructMeridional_tl=0
    v_tl=0
    h_edge_tl=0
    h_vertex_tl=0
    pv_edge_tl=0
    pv_vertex_tl=0
    pv_cell_tl=0
    circulation_tl=0
    vorticity_tl=0
    divergence_tl=0
    ke_tl=0
    vh_tl=0
    vorticity_cell_tl=0
    gradPVn_tl=0
    gradPVt_tl=0

    ! ===== allocating tend =====
    allocate(tend_u_tl(nVertLevels,nEdges))
    allocate(tend_h_tl(nVertLevels,nCells))
    allocate(tend_tracer_tl(nTracers,nVertLevels,nCells))
    tend_u_tl=0
    tend_h_tl=0
    tend_tracer_tl=0

  end subroutine var_allocation_tlm

  subroutine var_allocation_adj()
    implicit none

    integer :: dynamics_split

    ! ===== allocating state =====
    allocate(u_ad(2,nVertLevels,nEdges))
    u_ad=0
    allocate(h_ad(2,nVertLevels,nCells))
    h_ad=0
    allocate(tracers_ad(2,nTracers,nVertLevels,nCells))
    tracers_ad=0

    ! ===== allocating diag =====
    allocate(uReconstructX_ad(nVertLevels,nCells))
    allocate(uReconstructY_ad(nVertLevels,nCells))
    allocate(uReconstructZ_ad(nVertLevels,nCells))
    allocate(uReconstructZonal_ad(nVertLevels,nCells))
    allocate(uReconstructMeridional_ad(nVertLevels,nCells))
    allocate(v_ad(nVertLevels,nEdges))
    allocate(h_edge_ad(nVertLevels,nEdges))
    allocate(h_vertex_ad(nVertLevels,nVertices))
    allocate(pv_edge_ad(nVertLevels,nEdges))
    allocate(pv_vertex_ad(nVertLevels,nVertices))
    allocate(pv_cell_ad(nVertLevels,nCells))
    allocate(circulation_ad(nVertLevels,nVertices))
    allocate(vorticity_ad(nVertLevels,nVertices))
    allocate(divergence_ad(nVertLevels,nCells))
    allocate(ke_ad(nVertLevels,nCells))
    allocate(vh_ad(nVertLevels,nEdges))
    allocate(vorticity_cell_ad(nVertLevels,nCells))
    allocate(gradPVn_ad(nVertLevels,nEdges))
    allocate(gradPVt_ad(nVertLevels,nEdges))

    uReconstructX_ad=0
    uReconstructY_ad=0
    uReconstructZ_ad=0
    uReconstructZonal_ad=0
    uReconstructMeridional_ad=0
    v_ad=0
    h_edge_ad=0
    h_vertex_ad=0
    pv_edge_ad=0
    pv_vertex_ad=0
    pv_cell_ad=0
    circulation_ad=0
    vorticity_ad=0
    divergence_ad=0
    ke_ad=0
    vh_ad=0
    vorticity_cell_ad=0
    gradPVn_ad=0
    gradPVt_ad=0

    ! ===== allocating tend =====
    allocate(tend_u_ad(nVertLevels,nEdges))
    allocate(tend_h_ad(nVertLevels,nCells))
    allocate(tend_tracer_ad(nTracers,nVertLevels,nCells))
    tend_u_ad=0
    tend_h_ad=0
    tend_tracer_ad=0

  end subroutine var_allocation_adj  

  subroutine sw_mpas_init_block()
    implicit none

    call sw_compute_solve_diagnostics(config_dt, nCells, nEdges, nVertices,       &
         maxEdges, maxEdges2, nVertLevels, vertexDegree,                          &
         config_monotonic, config_thickness_adv_order, config_apvm_upwinding,     &
         h(1,:,:), u(1,:,:), v, vh, h_edge, h_vertex, circulation,                &
         vorticity, divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell,  &
         gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge,         &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell,   &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell,      &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call compute_mesh_scaling(nCells, nEdges, config_h_ScaleWithMesh, &
         cellsOnEdge, meshDensity, meshScalingDel2, meshScalingDel4)

    call mpas_initialize_vectors(nCells, nEdges, maxEdges, R3, &
         verticesOnEdge, cellsOnEdge, edgesOnCell, xCell, yCell, zCell, &
         xEdge, yEdge, zEdge, localVerticalUnitVectors, edgeNormalVectors, &
         cellTangentPlane, on_a_sphere, is_periodic, x_period, y_period)

    call mpas_init_reconstruct(nCells, nEdges, maxEdges, R3, edgesOnCell, &
         nEdgesOnCell, xCell, yCell, zCell, xEdge, yEdge, zEdge, &
         edgeNormalVectors, cellTangentPlane, coeffs_reconstruct, &
         is_periodic, x_period, y_period)

    call mpas_reconstruct_2d(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), uReconstructX, uReconstructY, uReconstructZ,&
         uReconstructZonal, uReconstructMeridional)

  end subroutine sw_mpas_init_block

  subroutine sw_rk4()
    implicit none

    real(kind=RKIND) :: dt
    real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

    ! ----- local variables -----
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis
    real(kind=RKIND), dimension(nTracers,nVertLevels,nCells) :: tracersProvis
    ! ----- end local variables -----

    dt=config_dt
    ! ----- restarting combo -----

    call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(1,:,:), u(1,:,:), v, vh, h_edge, h_vertex, circulation, vorticity, &
         divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
         gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call mpas_reconstruct_2d(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), uReconstructX, uReconstructY, uReconstructZ,&
         uReconstructZonal, uReconstructMeridional)
    ! ----- end restarting combo -----
    
    uProvis=u(1,:,:)
    hProvis=h(1,:,:)
    tracersProvis=tracers(1,:,:,:)
    
    u(2,:,:)=u(1,:,:)
    h(2,:,:)=h(1,:,:)
    tracers(2,:,:,:)=tracers(1,:,:,:)

    rk_weights(1) = dt/6.
    rk_weights(2) = dt/3.
    rk_weights(3) = dt/3.
    rk_weights(4) = dt/6.

    rk_substep_weights(1) = dt/2.
    rk_substep_weights(2) = dt/2.
    rk_substep_weights(3) = dt
    rk_substep_weights(4) = 0.

    do rk_step=1, 4
       call sw_compute_tend(nCells, nEdges, nVertices, nVertLevels, &
            maxEdges, maxEdges2, vertexDegree, config_bottom_drag, &
            config_wind_stress,config_h_mom_eddy_visc2, &
            config_h_mom_eddy_visc4, hProvis, uProvis, v, h_edge, circulation, &
            vorticity, divergence, ke, pv_edge, vh, weightsOnEdge, &
            kiteAreasOnVertex, cellsOnEdge, cellsOnVertex, verticesOnEdge, &
            nEdgesOnCell, edgesOnCell, nEdgesOnEdge, edgesOnEdge, edgesOnVertex, &
            dcEdge, dvEdge, areaCell, areaTriangle, h_s, fVertex, fEdge, u_src, &
            meshScalingDel2, meshScalingDel4, tend_h, tend_u)
       ! call sw_compute_scalar_tend() ! omitted in this version
       ! call sw_enforce_boundary_edge() ! global model, no lateral boundaries

       if (rk_step<4) then
          uProvis(:,:) = u(1,:,:) + rk_substep_weights(rk_step) * tend_u(:,:)
          hProvis(:,:) = h(1,:,:) + rk_substep_weights(rk_step) * tend_h(:,:)
          ! tracersProvis = ...
       end if

       if (config_test_case == 1) then
          uProvis(:,:) = u(1,:,:) ! For case 1, wind field should be fixed
       end if

       call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
            maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
            config_thickness_adv_order, config_apvm_upwinding, &
            hProvis, uProvis, v, vh, h_edge, h_vertex, circulation, vorticity, &
            divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
            gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
            cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
            nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
            areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

       u(2,:,:)=u(2,:,:) + rk_weights(rk_step) * tend_u
       h(2,:,:)=h(2,:,:) + rk_weights(rk_step) * tend_h
       ! tracersProvis = ...
    end do ! rk_step=1, 4

    if (config_test_case == 1) u(2,:,:)=u(1,:,:)
    call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(2,:,:), u(2,:,:), v, vh, h_edge, h_vertex, circulation, vorticity, &
         divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
         gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call mpas_reconstruct_2d(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(2,:,:), uReconstructX, uReconstructY, uReconstructZ,&
         uReconstructZonal, uReconstructMeridional)

    u(1,:,:)=u(2,:,:)
    h(1,:,:)=h(2,:,:)

    tracers(1,:,:,:)=tracers(2,:,:,:)
  end subroutine sw_rk4
  
  subroutine sw_rk4_tlm()
    implicit none

    real(kind=RKIND) :: dt
    real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

    ! ----- local variables -----
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis_tl
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis_tl
    real(kind=RKIND), dimension(nTracers,nVertLevels,nCells) :: tracersProvis
    real(kind=RKIND), dimension(nTracers,nVertLevels,nCells) :: tracersProvis_tl
    ! ----- end local variables -----
    
    dt=config_dt

    ! ----- restarting -----
    call sw_compute_solve_diagnostics_tlm(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(1,:,:), h_tl(1,:,:), u(1,:,:), u_tl(1,:,:), v, v_tl, vh, vh_tl, &
         h_edge, h_edge_tl, h_vertex, h_vertex_tl, circulation, circulation_tl, &
         vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
         pv_edge, pv_edge_tl, pv_vertex, pv_vertex_tl, pv_cell, pv_cell_tl, &
         vorticity_cell, vorticity_cell_tl, gradPVn, gradPVn_tl, gradPVt, gradPVt_tl, &
         weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call mpas_reconstruct_2d_tlm(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), u_tl(1,:,:), &
         uReconstructX, uReconstructX_tl, uReconstructY, uReconstructY_tl, &
         uReconstructZ, uReconstructZ_tl, uReconstructZonal, uReconstructZonal_tl, &
         uReconstructMeridional, uReconstructMeridional_tl)
    ! ----- end restarting -----

    uProvis_tl=u_tl(1,:,:)
    uProvis=u(1,:,:)
    hProvis_tl=h_tl(1,:,:)
    hProvis=h(1,:,:)
    tracersProvis_tl=tracers_tl(1,:,:,:)
    tracersProvis=tracers(1,:,:,:)

    u_tl(2,:,:)=u_tl(1,:,:)
    u(2,:,:)=u(1,:,:)
    h_tl(2,:,:)=h_tl(1,:,:)
    h(2,:,:)=h(1,:,:)
    tracers_tl(2,:,:,:)=tracers_tl(1,:,:,:)
    tracers(2,:,:,:)=tracers(1,:,:,:)

    rk_weights(1) = dt/6.
    rk_weights(2) = dt/3.
    rk_weights(3) = dt/3.
    rk_weights(4) = dt/6.

    rk_substep_weights(1) = dt/2.
    rk_substep_weights(2) = dt/2.
    rk_substep_weights(3) = dt
    rk_substep_weights(4) = 0.

    do rk_step=1, 4
       call sw_compute_tend_tlm(nCells, nEdges, nVertices, nVertLevels, &
            maxEdges, maxEdges2, vertexDegree, config_bottom_drag, &
            config_wind_stress,config_h_mom_eddy_visc2, &
            config_h_mom_eddy_visc4, hProvis, hProvis_tl, uProvis, uProvis_tl, &
            v, v_tl, h_edge, h_edge_tl, circulation, circulation_tl, &
            vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
            pv_edge, pv_edge_tl, vh, vh_tl, weightsOnEdge, kiteAreasOnVertex, &
            cellsOnEdge, cellsOnVertex, verticesOnEdge, nEdgesOnCell, &
            edgesOnCell, nEdgesOnEdge, edgesOnEdge, edgesOnVertex, &
            dcEdge, dvEdge, areaCell, areaTriangle, h_s, fVertex, fEdge, u_src, &
            meshScalingDel2, meshScalingDel4, &
            tend_h, tend_h_tl, tend_u, tend_u_tl)
       ! call sw_compute_scalar_tend() ! omitted in this version
       ! call sw_enforce_boundary_edge() ! global model, no lateral boundaries

       if (rk_step<4) then
          uProvis_tl(:,:) = u_tl(1,:,:) + rk_substep_weights(rk_step) * tend_u_tl(:,:)
          uProvis(:,:) = u(1,:,:) + rk_substep_weights(rk_step) * tend_u(:,:)
          hProvis_tl(:,:) = h_tl(1,:,:) + rk_substep_weights(rk_step) * tend_h_tl(:,:)
          hProvis(:,:) = h(1,:,:) + rk_substep_weights(rk_step) * tend_h(:,:)
          ! tracersProvis = ...
       end if

       if (config_test_case == 1) then
          uProvis_tl(:,:) = u_tl(1,:,:)
          uProvis(:,:) = u(1,:,:) ! For case 1, wind field should be fixed
       end if

       call sw_compute_solve_diagnostics_tlm(dt, nCells, nEdges, nVertices, &
            maxEdges, maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
            config_thickness_adv_order, config_apvm_upwinding, &
            hProvis, hProvis_tl, uProvis, uProvis_tl, v, v_tl, vh, vh_tl, &
            h_edge, h_edge_tl, h_vertex, h_vertex_tl, circulation, circulation_tl, &
            vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
            pv_edge, pv_edge_tl, pv_vertex, pv_vertex_tl, pv_cell, pv_cell_tl, &
            vorticity_cell, vorticity_cell_tl, gradPVn, gradPVn_tl, &
            gradPVt, gradPVt_tl, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
            cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
            nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
            areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

       u_tl(2,:,:)=u_tl(2,:,:) + rk_weights(rk_step) * tend_u_tl
       u(2,:,:)=u(2,:,:) + rk_weights(rk_step) * tend_u
       h_tl(2,:,:)=h_tl(2,:,:) + rk_weights(rk_step) * tend_h_tl
       h(2,:,:)=h(2,:,:) + rk_weights(rk_step) * tend_h
       ! tracersProvis = ...
    end do ! rk_step=1, 4

    if (config_test_case == 1) then
       u_tl(2,:,:)=u_tl(1,:,:)
       u(2,:,:)=u(1,:,:)
    end if

    call sw_compute_solve_diagnostics_tlm(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(2,:,:), h_tl(2,:,:), u(2,:,:), u_tl(2,:,:), v, v_tl, vh, vh_tl, &
         h_edge, h_edge_tl, h_vertex, h_vertex_tl, circulation, circulation_tl, &
         vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
         pv_edge, pv_edge_tl, pv_vertex, pv_vertex_tl, pv_cell, pv_cell_tl, &
         vorticity_cell, vorticity_cell_tl, gradPVn, gradPVn_tl, gradPVt, gradPVt_tl, &
         weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call mpas_reconstruct_2d_tlm(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(2,:,:), u_tl(2,:,:), &
         uReconstructX, uReconstructX_tl, uReconstructY, uReconstructY_tl, &
         uReconstructZ, uReconstructZ_tl, uReconstructZonal, uReconstructZonal_tl, &
         uReconstructMeridional, uReconstructMeridional_tl)

    u_tl(1,:,:)=u_tl(2,:,:)
    u(1,:,:)=u(2,:,:)
    h_tl(1,:,:)=h_tl(2,:,:)
    h(1,:,:)=h(2,:,:)

    tracers_tl(1,:,:,:)=tracers_tl(2,:,:,:)
    tracers(1,:,:,:)=tracers(2,:,:,:)
  end subroutine sw_rk4_tlm

  subroutine sw_rk4_adj()
    implicit none

    real(kind=RKIND) :: dt
    real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

    ! ----- local variables -----
    real(kind=RKIND), dimension(5,nVertLevels,nEdges) :: uProvis
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis_ad
    real(kind=RKIND), dimension(5,nVertLevels,nCells) :: hProvis
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis_ad
    !real(kind=RKIND), dimension(nTracers,nVertLevels,nCells) :: tracersProvis
    !real(kind=RKIND), dimension(nTracers,nVertLevels,nCells) :: tracersProvis_tl
    ! ----- end local variables -----

    uProvis_ad=0.
    hProvis_ad=0.
    
    dt=config_dt

    rk_weights(1) = dt/6.
    rk_weights(2) = dt/3.
    rk_weights(3) = dt/3.
    rk_weights(4) = dt/6.

    rk_substep_weights(1) = dt/2.
    rk_substep_weights(2) = dt/2.
    rk_substep_weights(3) = dt
    rk_substep_weights(4) = 0.

    ! ========== forward calculation ==========

    call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(1,:,:), u(1,:,:), v, vh, h_edge, h_vertex, circulation, vorticity, &
         divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
         gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call mpas_reconstruct_2d(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), uReconstructX, uReconstructY, uReconstructZ,&
         uReconstructZonal, uReconstructMeridional)

    uProvis(1,:,:)=u(1,:,:)
    hProvis(1,:,:)=h(1,:,:)
    ! tracersProvis=tracers(1,:,:,:)

    u(2,:,:)=u(1,:,:)
    h(2,:,:)=h(1,:,:)
    ! tracers(2,:,:,:)=tracers(1,:,:,:)

    do rk_step=1, 4
       call sw_compute_tend(nCells, nEdges, nVertices, nVertLevels, &
            maxEdges, maxEdges2, vertexDegree, config_bottom_drag, &
            config_wind_stress,config_h_mom_eddy_visc2, &
            config_h_mom_eddy_visc4, hProvis(rk_step,:,:), uProvis(rk_step,:,:), &
            v, h_edge, circulation, &
            vorticity, divergence, ke, pv_edge, vh, weightsOnEdge, &
            kiteAreasOnVertex, cellsOnEdge, cellsOnVertex, verticesOnEdge, &
            nEdgesOnCell, edgesOnCell, nEdgesOnEdge, edgesOnEdge, edgesOnVertex, &
            dcEdge, dvEdge, areaCell, areaTriangle, h_s, fVertex, fEdge, u_src, &
            meshScalingDel2, meshScalingDel4, tend_h, tend_u)
       ! call sw_compute_scalar_tend() ! omitted in this version
       ! call sw_enforce_boundary_edge() ! global model, no lateral boundaries

       if (rk_step<4) then
          uProvis(rk_step+1,:,:) = u(1,:,:) + rk_substep_weights(rk_step) * tend_u(:,:)
          hProvis(rk_step+1,:,:) = h(1,:,:) + rk_substep_weights(rk_step) * tend_h(:,:)
          ! tracersProvis = ...
       else
          uProvis(rk_step+1,:,:)=uProvis(rk_step,:,:)
          hProvis(rk_step+1,:,:)=hProvis(rk_step,:,:)
       end if

       if (config_test_case == 1) then
          uProvis(rk_step+1,:,:) = u(1,:,:) ! For case 1, wind field should be fixed
       end if

       call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
            maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
            config_thickness_adv_order, config_apvm_upwinding, &
            hProvis(rk_step+1,:,:), uProvis(rk_step+1,:,:), v, vh, h_edge, h_vertex, circulation, vorticity, &
            divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
            gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
            cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
            nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
            areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

       u(2,:,:)=u(2,:,:) + rk_weights(rk_step) * tend_u
       h(2,:,:)=h(2,:,:) + rk_weights(rk_step) * tend_h
       ! tracersProvis = ...
    end do ! rk_step=1, 4

    ! ========== end forward calculation ==========

    h_ad(2,:,:)=h_ad(1,:,:)
    h_ad(1,:,:)=0
    u_ad(2,:,:)=u_ad(1,:,:)
    u_ad(1,:,:)=0

    call mpas_reconstruct_2d_adj(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(2,:,:), u_ad(2,:,:), &
         uReconstructX, uReconstructX_ad, uReconstructY, uReconstructY_ad, &
         uReconstructZ, uReconstructZ_ad, uReconstructZonal, uReconstructZonal_ad, &
         uReconstructMeridional, uReconstructMeridional_ad)

    call sw_compute_solve_diagnostics_adj(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(2,:,:), h_ad(2,:,:), u(2,:,:), u_ad(2,:,:), v, v_ad, vh, vh_ad, &
         h_edge, h_edge_ad, h_vertex, h_vertex_ad, circulation, circulation_ad, &
         vorticity, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
         pv_edge, pv_edge_ad, pv_vertex, pv_vertex_ad, pv_cell, pv_cell_ad, &
         vorticity_cell, vorticity_cell_ad, gradPVn, gradPVn_ad, gradPVt, gradPVt_ad, &
         weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    if (config_test_case == 1) then
       u_ad(1,:,:)=u_ad(1,:,:)+u_ad(2,:,:)
       u_ad(2,:,:)=0.
    end if

    do rk_step=4, 1, -1
       tend_h_ad=tend_h_ad + rk_weights(rk_step)*h_ad(2,:,:)
       tend_u_ad=tend_u_ad + rk_weights(rk_step)*u_ad(2,:,:)

       call sw_compute_solve_diagnostics_adj(dt, nCells, nEdges, nVertices, &
            maxEdges, maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
            config_thickness_adv_order, config_apvm_upwinding, &
            hProvis(rk_step+1,:,:), hProvis_ad, uProvis(rk_step+1,:,:), uProvis_ad, &
            v, v_ad, vh, vh_ad, &
            h_edge, h_edge_ad, h_vertex, h_vertex_ad, circulation, circulation_ad, &
            vorticity, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
            pv_edge, pv_edge_ad, pv_vertex, pv_vertex_ad, pv_cell, pv_cell_ad, &
            vorticity_cell, vorticity_cell_ad, gradPVn, gradPVn_ad, &
            gradPVt, gradPVt_ad, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
            cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
            nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
            areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

       if (config_test_case==1) then
          u_ad(1,:,:)=u_ad(1,:,:)+uProvis_ad(:,:)
          uProvis_ad(:,:)=0.
       end if

       if (rk_step<4) then
          tend_h_ad(:,:)=tend_h_ad(:,:) + rk_substep_weights(rk_step)*hProvis_ad(:,:)
          h_ad(1,:,:)=h_ad(1,:,:)+hProvis_ad(:,:)
          hProvis_ad(:,:)=0.
          tend_u_ad(:,:)=tend_u_ad(:,:) + rk_substep_weights(rk_step)*uProvis_ad(:,:)
          u_ad(1,:,:)=u_ad(1,:,:)+uProvis_ad(:,:)
          uProvis_ad(:,:)=0.
       end if

       ! ----- forward calculation -----
       call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
            maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
            config_thickness_adv_order, config_apvm_upwinding, &
            hProvis(rk_step,:,:), uProvis(rk_step,:,:), v, vh, h_edge, h_vertex, circulation, &
            vorticity, &
            divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
            gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
            cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
            nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
            areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)
       ! ----- end forward calculation -----

       call sw_compute_tend_adj(nCells, nEdges, nVertices, nVertLevels, &
            maxEdges, maxEdges2, vertexDegree, config_bottom_drag, &
            config_wind_stress,config_h_mom_eddy_visc2, &
            config_h_mom_eddy_visc4, hProvis(rk_step,:,:), hProvis_ad, &
            uProvis(rk_step,:,:), uProvis_ad, &
            v, v_ad, h_edge, h_edge_ad, circulation, circulation_ad, &
            vorticity, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
            pv_edge, pv_edge_ad, vh, vh_ad, weightsOnEdge, kiteAreasOnVertex, &
            cellsOnEdge, cellsOnVertex, verticesOnEdge, nEdgesOnCell, &
            edgesOnCell, nEdgesOnEdge, edgesOnEdge, edgesOnVertex, &
            dcEdge, dvEdge, areaCell, areaTriangle, h_s, fVertex, fEdge, u_src, &
            meshScalingDel2, meshScalingDel4, &
            tend_h, tend_h_ad, tend_u, tend_u_ad)
    end do

    h_ad(1,:,:)=h_ad(1,:,:)+h_ad(2,:,:)
    h_ad(2,:,:)=0.
    u_ad(1,:,:)=u_ad(1,:,:)+u_ad(2,:,:)
    u_ad(2,:,:)=0.

    h_ad(1,:,:)=h_ad(1,:,:) + hProvis_ad
    hProvis_ad=0.
    u_ad(1,:,:)=u_ad(1,:,:) + uProvis_ad
    uProvis_ad=0.

    call mpas_reconstruct_2d_adj(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), u_ad(1,:,:), &
         uReconstructX, uReconstructX_ad, uReconstructY, uReconstructY_ad, &
         uReconstructZ, uReconstructZ_ad, uReconstructZonal, uReconstructZonal_ad, &
         uReconstructMeridional, uReconstructMeridional_ad)

    call sw_compute_solve_diagnostics_adj(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(1,:,:), h_ad(1,:,:), u(1,:,:), u_ad(1,:,:), v, v_ad, vh, vh_ad, &
         h_edge, h_edge_ad, h_vertex, h_vertex_ad, circulation, circulation_ad, &
         vorticity, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
         pv_edge, pv_edge_ad, pv_vertex, pv_vertex_ad, pv_cell, pv_cell_ad, &
         vorticity_cell, vorticity_cell_ad, gradPVn, gradPVn_ad, gradPVt, gradPVt_ad, &
         weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

  end subroutine sw_rk4_adj

  include "mpas_sw_testcases.f90"

  subroutine link_tlm_adj()
    implicit none

    real(kind=RKIND) :: dt
    real (kind=RKIND), dimension(4) :: rk_weights, rk_substep_weights

    real(kind=RKIND), dimension(5,nVertLevels,nEdges) :: uProvis_fw
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis_tl
    real(kind=RKIND), dimension(nVertLevels,nEdges) :: uProvis_ad
    real(kind=RKIND), dimension(5,nVertLevels,nCells) :: hProvis_fw
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis_tl
    real(kind=RKIND), dimension(nVertLevels,nCells) :: hProvis_ad

    pt=1.E-3
    dt=config_dt
    rk_weights(1) = dt/6.
    rk_weights(2) = dt/3.
    rk_weights(3) = dt/3.
    rk_weights(4) = dt/6.

    rk_substep_weights(1) = dt/2.
    rk_substep_weights(2) = dt/2.
    rk_substep_weights(3) = dt
    rk_substep_weights(4) = 0.

    ! ===== TLM =====
    call sw_compute_solve_diagnostics_tlm(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(1,:,:), h_tl(1,:,:), u(1,:,:), u_tl(1,:,:), v, v_tl, vh, vh_tl, &
         h_edge, h_edge_tl, h_vertex, h_vertex_tl, circulation, circulation_tl, &
         vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
         pv_edge, pv_edge_tl, pv_vertex, pv_vertex_tl, pv_cell, pv_cell_tl, &
         vorticity_cell, vorticity_cell_tl, gradPVn, gradPVn_tl, gradPVt, gradPVt_tl, &
         weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)

    call mpas_reconstruct_2d_tlm(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), u_tl(1,:,:), &
         uReconstructX, uReconstructX_tl, uReconstructY, uReconstructY_tl, &
         uReconstructZ, uReconstructZ_tl, uReconstructZonal, uReconstructZonal_tl, &
         uReconstructMeridional, uReconstructMeridional_tl)

    uProvis_tl=u_tl(1,:,:)
    uProvis=u(1,:,:)
    uProvis_fw(1,:,:)=u(1,:,:)
    hProvis_tl=h_tl(1,:,:)
    hProvis=h(1,:,:)
    hProvis_fw(1,:,:)=h(1,:,:)

    u_tl(2,:,:)=u_tl(1,:,:)
    u(2,:,:)=u(1,:,:)
    h_tl(2,:,:)=h_tl(1,:,:)
    h(2,:,:)=h(1,:,:)

    do rk_step=1, 1 !4
       call sw_compute_tend_tlm(nCells, nEdges, nVertices, nVertLevels, &
            maxEdges, maxEdges2, vertexDegree, config_bottom_drag, &
            config_wind_stress,config_h_mom_eddy_visc2, &
            config_h_mom_eddy_visc4, hProvis, hProvis_tl, uProvis, uProvis_tl, &
            v, v_tl, h_edge, h_edge_tl, circulation, circulation_tl, &
            vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
            pv_edge, pv_edge_tl, vh, vh_tl, weightsOnEdge, kiteAreasOnVertex, &
            cellsOnEdge, cellsOnVertex, verticesOnEdge, nEdgesOnCell, &
            edgesOnCell, nEdgesOnEdge, edgesOnEdge, edgesOnVertex, &
            dcEdge, dvEdge, areaCell, areaTriangle, h_s, fVertex, fEdge, u_src, &
            meshScalingDel2, meshScalingDel4, &
            tend_h, tend_h_tl, tend_u, tend_u_tl)

       if (rk_step<4) then
          uProvis_tl(:,:) = u_tl(1,:,:) + rk_substep_weights(rk_step) * tend_u_tl(:,:)
          uProvis(:,:) = u(1,:,:) + rk_substep_weights(rk_step) * tend_u(:,:)
          hProvis_tl(:,:) = h_tl(1,:,:) + rk_substep_weights(rk_step) * tend_h_tl(:,:)
          hProvis(:,:) = h(1,:,:) + rk_substep_weights(rk_step) * tend_h(:,:)
       end if
       
    end do
    
    ! ===== end TLM =====
    LHS= SUM(uProvis_tl**2) + SUM(hProvis_tl**2) + &
         SUM(tend_h_tl**2) + SUM(tend_u_tl**2) + &
         SUM(u_tl(2,:,:)**2) + SUM(h_tl(2,:,:)**2) + &
         SUM(v_tl**2) + &
         SUM(h_edge_tl**2) + &
         SUM(h_vertex_tl**2) + &
         SUM(circulation_tl**2) + &
         SUM(vorticity_tl**2) + &
         SUM(divergence_tl**2) + &
         SUM(ke_tl**2) + &
         SUM(pv_edge_tl**2) + &
         SUM(pv_vertex_tl**2) + &
         SUM(pv_cell_tl**2) + &
         SUM(vorticity_cell_tl**2) + &
         SUM(gradPVn_tl**2) + &
         SUM(gradPVt_tl**2) + &
         SUM(uReconstructX_tl**2) + &
         SUM(uReconstructY_tl**2) + &
         SUM(uReconstructZ_tl**2) + &
         SUM(uReconstructZonal_tl**2) + &
         SUM(uReconstructMeridional_tl**2)

    uProvis_ad=uProvis_tl
    hProvis_ad=hProvis_tl
    tend_h_ad=tend_h_tl
    tend_u_ad=tend_u_tl
    u_ad(2,:,:)=u_tl(2,:,:)
    h_ad(2,:,:)=h_tl(2,:,:)
    v_ad=v_tl
    h_edge_ad=h_edge_tl
    h_vertex_ad=h_vertex_tl
    circulation_ad=circulation_tl
    vorticity_ad=vorticity_tl
    divergence_ad=divergence_tl
    ke_ad=ke_tl
    pv_edge_ad=pv_edge_tl
    pv_vertex_ad=pv_vertex_tl
    pv_cell_ad=pv_cell_tl
    vorticity_cell_ad=vorticity_cell_tl
    gradPVn_ad=gradPVn_tl
    gradPVt_ad=gradPVt_tl
    uReconstructX_ad=uReconstructX_tl
    uReconstructY_ad=uReconstructY_tl
    uReconstructZ_ad=uReconstructZ_tl
    uReconstructZonal_ad=uReconstructZonal_tl
    uReconstructMeridional_ad=uReconstructMeridional_tl

    ! ===== ADJ =====

    do rk_step=1, 1, -1

       if (rk_step<4) then
          tend_h_ad(:,:)=tend_h_ad(:,:) + rk_substep_weights(rk_step)*hProvis_ad(:,:)
          h_ad(1,:,:)=h_ad(1,:,:)+hProvis_ad(:,:)
          hProvis_ad(:,:)=0.
          tend_u_ad(:,:)=tend_u_ad(:,:) + rk_substep_weights(rk_step)*uProvis_ad(:,:)
          u_ad(1,:,:)=u_ad(1,:,:)+uProvis_ad(:,:)
          uProvis_ad(:,:)=0.
       end if
       
       ! ----- forward calculation -----
       call sw_compute_solve_diagnostics(dt, nCells, nEdges, nVertices, maxEdges, &
            maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
            config_thickness_adv_order, config_apvm_upwinding, &
            hProvis_fw(rk_step,:,:), uProvis_fw(rk_step,:,:), v, vh, h_edge, h_vertex, circulation, &
            vorticity, &
            divergence, ke, pv_edge, pv_vertex, pv_cell, vorticity_cell, &
            gradPVn, gradPVt, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
            cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
            nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
            areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)
       ! ----- end forward calculation -----

       call sw_compute_tend_adj(nCells, nEdges, nVertices, nVertLevels, &
            maxEdges, maxEdges2, vertexDegree, config_bottom_drag, &
            config_wind_stress,config_h_mom_eddy_visc2, &
            config_h_mom_eddy_visc4, hProvis_fw(rk_step,:,:), hProvis_ad, &
            uProvis_fw(rk_step,:,:), uProvis_ad, &
            v, v_ad, h_edge, h_edge_ad, circulation, circulation_ad, &
            vorticity, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
            pv_edge, pv_edge_ad, vh, vh_ad, weightsOnEdge, kiteAreasOnVertex, &
            cellsOnEdge, cellsOnVertex, verticesOnEdge, nEdgesOnCell, &
            edgesOnCell, nEdgesOnEdge, edgesOnEdge, edgesOnVertex, &
            dcEdge, dvEdge, areaCell, areaTriangle, h_s, fVertex, fEdge, u_src, &
            meshScalingDel2, meshScalingDel4, &
            tend_h, tend_h_ad, tend_u, tend_u_ad)
    end do ! rk_step=4, 1, -1
    
    h_ad(1,:,:)=h_ad(1,:,:)+h_ad(2,:,:)
    h_ad(2,:,:)=0.
    u_ad(1,:,:)=u_ad(1,:,:)+u_ad(2,:,:)
    u_ad(2,:,:)=0.

    h_ad(1,:,:)=h_ad(1,:,:) + hProvis_ad
    hProvis_ad=0.
    u_ad(1,:,:)=u_ad(1,:,:) + uProvis_ad
    uProvis_ad=0.
    
    call mpas_reconstruct_2d_adj(nCells, nEdges, maxEdges, nVertLevels, R3, &
         on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
         coeffs_reconstruct, u(1,:,:), u_ad(1,:,:), &
         uReconstructX, uReconstructX_ad, uReconstructY, uReconstructY_ad, &
         uReconstructZ, uReconstructZ_ad, uReconstructZonal, uReconstructZonal_ad, &
         uReconstructMeridional, uReconstructMeridional_ad)
    
    call sw_compute_solve_diagnostics_adj(dt, nCells, nEdges, nVertices, maxEdges, &
         maxEdges2, nVertLevels, vertexDegree, config_monotonic, &
         config_thickness_adv_order, config_apvm_upwinding, &
         h(1,:,:), h_ad(1,:,:), u(1,:,:), u_ad(1,:,:), v, v_ad, vh, vh_ad, &
         h_edge, h_edge_ad, h_vertex, h_vertex_ad, circulation, circulation_ad, &
         vorticity, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
         pv_edge, pv_edge_ad, pv_vertex, pv_vertex_ad, pv_cell, pv_cell_ad, &
         vorticity_cell, vorticity_cell_ad, gradPVn, gradPVn_ad, gradPVt, gradPVt_ad, &
         weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
         cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
         nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
         areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)
    ! ===== end ADJ =====

    
    RHS=SUM(u_ad(1,:,:)*u(1,:,:)*pt) + SUM(h_ad(1,:,:)*h(1,:,:)*pt)

    print *, LHS
    print *, RHS
    print *, (LHS-RHS)/LHS
    print *, (LHS-RHS)/RHS
    STOP
    
  end subroutine link_tlm_adj
end module mpas_sw_module
