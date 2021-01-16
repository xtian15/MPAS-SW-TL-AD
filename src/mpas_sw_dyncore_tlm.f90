subroutine sw_compute_solve_diagnostics_tlm(dt, nCells, nEdges, nVertices, &
     maxEdges, maxEdges2, nVertLevels, vertexDegree, &
     config_monotonic, config_thickness_adv_order, config_apvm_upwinding, &
     h, h_tl, u, u_tl, v, v_tl, vh, vh_tl, h_edge, h_edge_tl, &
     h_vertex, h_vertex_tl, circulation, circulation_tl, &
     vorticity, vorticity_tl, divergence, divergence_tl, ke, ke_tl, &
     pv_edge, pv_edge_tl, pv_vertex, pv_vertex_tl, pv_cell, pv_cell_tl, &
     vorticity_cell, vorticity_cell_tl,  gradPVn, gradPVn_tl, &
     gradPVt, gradPVt_tl, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
     cellsOnVertex, cellsOnCell, verticesOnEdge, nEdgesOnCell, edgesOnCell, &
     nEdgesOnEdge, edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, &
     areaTriangle, fVertex, deriv_two, boundaryEdge, boundaryCell)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute diagnostic fields used in the tendency computations
  !
  ! Input: grid - grid metadata
  !
  ! Output: s - computed diagnostics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)
  logical :: config_monotonic
  integer :: config_thickness_adv_order
  real (kind=RKIND) :: config_apvm_upwinding
  real (kind=RKIND), intent(in) :: dt
  integer :: nCells, nEdges, nVertices, nVertLevels, vertexDegree
  integer :: maxEdges, maxEdges2

  real(kind=RKIND), dimension(nVertices) :: fVertex
  real(kind=RKIND), dimension(nEdges) :: dvEdge
  real(kind=RKIND), dimension(nEdges) :: dcEdge
  real(kind=RKIND), dimension(nCells) :: areaCell
  real(kind=RKIND), dimension(nVertices) :: areaTriangle
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: vh
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: vh_tl
  real(kind=RKIND), dimension(maxEdges2,nEdges) :: weightsOnEdge
  real(kind=RKIND), dimension(vertexDegree,nVertices) :: kiteAreasOnVertex
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge_tl
  real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: h
  real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: h_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: u
  real(kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: u_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: pv_vertex
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: pv_vertex_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: pv_cell
  real(kind=RKIND), dimension(nVertLevels,nCells) :: pv_cell_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVn
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVn_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVt
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVt_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: vorticity_cell
  real(kind=RKIND), dimension(nVertLevels,nCells) :: vorticity_cell_tl

  integer, dimension(2,nEdges) :: cellsOnEdge
  integer, dimension(vertexDegree,nVertices) :: cellsOnVertex
  integer, dimension(2,nEdges) :: verticesOnEdge
  integer, dimension(maxEdges,nCells) :: edgesOnCell
  integer, dimension(maxEdges2,nEdges) :: edgesOnEdge
  integer, dimension(vertexDegree,nVertices) :: edgesOnVertex
  integer, dimension(nVertLevels,nEdges) :: boundaryEdge
  integer, dimension(nVertLevels,nCells) :: boundaryCell
  integer, dimension(maxEdges,nCells) :: cellsOnCell

  integer, dimension(nCells) :: nEdgesOnCell
  integer, dimension(nEdges) :: nEdgesOnEdge
  real (kind=RKIND), dimension(15,2,nEdges) :: deriv_two

  ! ----- local variables -----
  real (kind=RKIND) :: r, h1, h2
  real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2
  real (kind=RKIND) :: d2fdx2_cell1_tl, d2fdx2_cell2_tl
  real (kind=RKIND) :: coef_3rd_order

  integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
  real (kind=RKIND) :: flux, vorticity_abs, workpv
  real (kind=RKIND) :: flux_tl, vorticity_abs_tl, workpv_tl
  ! ----- end local variables -----

  ! Find those cells that have an edge on the boundary
  boundaryCell(:,:) = 0
  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        if(boundaryEdge(k,iEdge).eq.1) then
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           boundaryCell(k,cell1) = 1
           boundaryCell(k,cell2) = 1
        endif
     enddo
  enddo

  ! Compute height on cell edges at velocity locations
  !   Namelist options control the order of accuracy of the reconstructed h_edge value
  coef_3rd_order = 0.
  if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
  if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

  if (config_thickness_adv_order == 2) then

     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        if (cell1 <= nCells .and. cell2 <= nCells) then
           do k = 1, nVertLevels
              h_edge_tl(k,iEdge) = 0.5 * (h_tl(k,cell1) + h_tl(k,cell2))
              h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
           end do
        end if
     end do

  else if (config_thickness_adv_order == 3) then
     
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        !-- if a cell not on the most outside ring of the halo
        if (cell1 <= nCells .and. cell2 <= nCells) then

           do k = 1, nVertLevels
              d2fdx2_cell1_tl = 0.0
              d2fdx2_cell1 = 0.0
              d2fdx2_cell2_tl = 0.0
              d2fdx2_cell2 = 0.0

              !-- if not a boundary cell
              if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then

                 d2fdx2_cell1_tl = deriv_two(1,1,iEdge) * h_tl(k,cell1)
                 d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                 d2fdx2_cell2_tl = deriv_two(1,2,iEdge) * h_tl(k,cell2)
                 d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                 !-- all edges of cell 1
                 do i = 1, nEdgesOnCell(cell1)
                    d2fdx2_cell1_tl = d2fdx2_cell1_tl + &
                         deriv_two(i+1,1,iEdge) * h_tl(k, cellsOnCell(i,cell1))
                    d2fdx2_cell1 = d2fdx2_cell1 + &
                         deriv_two(i+1,1,iEdge) * h(k, cellsOnCell(i,cell1))
                 end do

                 !-- all edges of cell 2
                 do i = 1, nEdgesOnCell(cell2)
                    d2fdx2_cell2_tl = d2fdx2_cell2_tl + &
                         deriv_two(i+1,2,iEdge) * h_tl(k, cellsOnCell(i,cell2))
                    d2fdx2_cell2 = d2fdx2_cell2 + &
                         deriv_two(i+1,2,iEdge) * h(k, cellsOnCell(i,cell2))
                 end do

              endif

              !-- if u > 0:
              if (u(k,iEdge) > 0) then
                 h_edge_tl(k,iEdge) =     &
                      0.5*(h_tl(k,cell1) + h_tl(k,cell2))      &
                      -(dcEdge(iEdge) **2) * (d2fdx2_cell1_tl + d2fdx2_cell2_tl) / 12.          &
                      -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1_tl - d2fdx2_cell2_tl)/12.
                 h_edge(k,iEdge) =     &
                      0.5*(h(k,cell1) + h(k,cell2))      &
                      -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                      -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
                 !-- else u <= 0:
              else
                 h_edge_tl(k,iEdge) =     &
                      0.5*(h_tl(k,cell1) + h_tl(k,cell2))      &
                      -(dcEdge(iEdge) **2) * (d2fdx2_cell1_tl + d2fdx2_cell2_tl) / 12.          &
                      +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1_tl - d2fdx2_cell2_tl)/12.
                 h_edge(k,iEdge) =     &
                      0.5*(h(k,cell1) + h(k,cell2))      &
                      -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
                      +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12.
              end if

           end do   ! do k
        end if      ! if (cell1 <=
     end do         ! do iEdge

  else  if (config_thickness_adv_order == 4) then

     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        !-- if a cell not on the most outside ring of the halo
        if (cell1 <= nCells .and. cell2 <= nCells) then

           do k = 1, nVertLevels

              d2fdx2_cell1_tl = 0.0
              d2fdx2_cell1 = 0.0
              d2fdx2_cell2_tl = 0.0
              d2fdx2_cell2 = 0.0

              !-- if not a boundary cell
              if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then

                 d2fdx2_cell1_tl = deriv_two(1,1,iEdge) * h_tl(k,cell1)
                 d2fdx2_cell1 = deriv_two(1,1,iEdge) * h(k,cell1)
                 d2fdx2_cell2_tl = deriv_two(1,2,iEdge) * h_tl(k,cell2)
                 d2fdx2_cell2 = deriv_two(1,2,iEdge) * h(k,cell2)

                 !-- all edges of cell 1
                 do i = 1, nEdgesOnCell(cell1)
                    d2fdx2_cell1_tl = d2fdx2_cell1_tl + &
                         deriv_two(i+1,1,iEdge) * h_tl(k, cellsOnCell(i,cell1))
                    d2fdx2_cell1 = d2fdx2_cell1 + &
                         deriv_two(i+1,1,iEdge) * h(k, cellsOnCell(i,cell1))
                 end do

                 !-- all edges of cell 2
                 do i = 1, nEdgesOnCell(cell2)
                    d2fdx2_cell2_tl = d2fdx2_cell2_tl + &
                         deriv_two(i+1,2,iEdge) * h_tl(k, cellsOnCell(i,cell2))
                    d2fdx2_cell2 = d2fdx2_cell2 + &
                         deriv_two(i+1,2,iEdge) * h(k, cellsOnCell(i,cell2))
                 end do

              endif

              h_edge_tl(k,iEdge) =   &
                   0.5*(h_tl(k,cell1) + h_tl(k,cell2))      &
                   -(dcEdge(iEdge) **2) * (d2fdx2_cell1_tl + d2fdx2_cell2_tl) / 12.
              h_edge(k,iEdge) =   &
                   0.5*(h(k,cell1) + h(k,cell2))      &
                   -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.

           end do   ! do k
        end if      ! if (cell1 <=
     end do         ! do iEdge

  endif   ! if(config_thickness_adv_order == 2)

  ! set the velocity in the nEdges+1 slot to zero, this is a dummy address
  !    used to when reading for edges that do not exist
  !u(:,nEdges+1) = 0.0

  ! Compute circulation and relative vorticity at each vertex
  circulation_tl(:,:) = 0.0
  circulation(:,:) = 0.0
  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        circulation_tl(k,verticesOnEdge(1,iEdge)) = circulation_tl(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u_tl(k,iEdge)
        circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)

        circulation_tl(k,verticesOnEdge(2,iEdge)) = circulation_tl(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u_tl(k,iEdge)
        circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
     end do
  end do
  do iVertex = 1, nVertices
     do k = 1, nVertLevels
        vorticity_tl(k,iVertex) = circulation_tl(k,iVertex) / areaTriangle(iVertex)
        vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
     end do
  end do

  ! Compute the divergence at each cell center
  divergence_tl(:,:) = 0.0
  divergence(:,:) = 0.0
  do iEdge = 1, nEdges
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     if (cell1 <= nCells) then
        do k = 1, nVertLevels
           divergence_tl(k,cell1) = divergence_tl(k,cell1) + u_tl(k,iEdge)*dvEdge(iEdge)
           divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
        enddo
     endif
     if(cell2 <= nCells) then
        do k = 1, nVertLevels
           divergence_tl(k,cell2) = divergence_tl(k,cell2) - u_tl(k,iEdge)*dvEdge(iEdge)
           divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
        enddo
     end if
  end do
  do iCell = 1, nCells
     r = 1.0 / areaCell(iCell)
     do k = 1, nVertLevels
        divergence_tl(k,iCell) = divergence_tl(k,iCell) * r
        divergence(k,iCell) = divergence(k,iCell) * r
     enddo
  enddo

  ! Compute kinetic energy in each cell
  ke_tl(:,:)=0.
  ke(:,:) = 0.0
  do iCell = 1, nCells
     do i = 1, nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        do k = 1, nVertLevels
           ke_tl(k,iCell) = ke_tl(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) &
                * (2.*u(k,iEdge)*u_tl(k,iEdge))
           ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
        end do
     end do
     do k = 1, nVertLevels
        ke_tl(k,iCell) = ke_tl(k,iCell) / areaCell(iCell)
        ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
     end do
  end do

  ! Compute v (tangential) velocities
  v_tl(:,:)=0.
  v(:,:) = 0.0
  do iEdge = 1,nEdges
     do i = 1, nEdgesOnEdge(iEdge)
        eoe = edgesOnEdge(i,iEdge)
        do k = 1,nVertLevels
           v_tl(k,iEdge) = v_tl(k,iEdge) + weightsOnEdge(i,iEdge) * u_tl(k, eoe)
           v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
        end do
     end do
  end do

  ! Compute mass fluxes tangential to each edge (i.e., through the faces of dual grid cells)
  !vh(:,:) = 0.0
  !do iEdge = 1, nEdges
  !   do j = 1, nEdgesOnEdge(iEdge)
  !      eoe = edgesOnEdge(j,iEdge)
  !      do k = 1, nVertLevels
  !         vh(k,iEdge) = vh(k,iEdge) + weightsOnEdge(j,iEdge) * u(k,eoe) * h_edge(k,eoe)
  !      end do
  !   end do
  !end do

  ! Compute height at vertices, pv at vertices, and average pv to edge locations
  !  ( this computes pv_vertex at all vertices bounding real cells and distance-1 ghost cells )
  do iVertex = 1,nVertices
     do k = 1, nVertLevels
        h_vertex_tl(k,iVertex) = 0.0
        h_vertex(k,iVertex) = 0.0
        do i = 1, vertexDegree
           h_vertex_tl(k,iVertex) = h_vertex_tl(k,iVertex) + h_tl(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
           h_vertex(k,iVertex) = h_vertex(k,iVertex) + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
        end do
        h_vertex_tl(k,iVertex) = h_vertex_tl(k,iVertex) / areaTriangle(iVertex)
        h_vertex(k,iVertex) = h_vertex(k,iVertex) / areaTriangle(iVertex)

        pv_vertex_tl(k,iVertex) = vorticity_tl(k,iVertex) / h_vertex(k,iVertex) &
             - h_vertex_tl(k,iVertex)*(fVertex(iVertex) + vorticity(k,iVertex)) &
             / (h_vertex(k,iVertex)**2)
        pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex(k,iVertex)
     end do
  end do

  ! Compute gradient of PV in the tangent direction
  !   ( this computes gradPVt at all edges bounding real cells and distance-1 ghost cells )
  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        gradPVt_tl(k,iEdge) = (pv_vertex_tl(k,verticesOnEdge(2,iEdge)) - pv_vertex_tl(k,verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
        gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
     enddo
  enddo

  ! Compute pv at the edges
  !   ( this computes pv_edge at all edges bounding real cells )
  pv_edge_tl(:,:) = 0.0
  pv_edge(:,:) = 0.0
  do iVertex = 1,nVertices
     do i = 1, vertexDegree
        iEdge = edgesOnVertex(i,iVertex)
        do k = 1, nVertLevels
           pv_edge_tl(k,iEdge) =  pv_edge_tl(k,iEdge)  + 0.5 * pv_vertex_tl(k,iVertex)
           pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
        end do
     end do
  end do

  ! Modify PV edge with upstream bias.
  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        pv_edge_tl(k,iEdge) = pv_edge_tl(k,iEdge) &
             - config_apvm_upwinding * v_tl(k,iEdge) * dt * gradPVt(k,iEdge) &
             - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt_tl(k,iEdge)
        pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
     enddo
  enddo

  ! Compute pv at cell centers
  !    ( this computes pv_cell for all real cells and distance-1 ghost cells )
  pv_cell_tl(:,:) = 0.0
  pv_cell(:,:) = 0.0
  vorticity_cell_tl(:,:) = 0.0
  vorticity_cell(:,:) = 0.0
  do iVertex = 1, nVertices
     do i = 1, vertexDegree
        iCell = cellsOnVertex(i,iVertex)
        if (iCell <= nCells) then
           do k = 1, nVertLevels
              pv_cell_tl(k,iCell) = pv_cell_tl(k,iCell) &
                   + kiteAreasOnVertex(i, iVertex) * pv_vertex_tl(k, iVertex) &
                   / areaCell(iCell)
              pv_cell(k,iCell) = pv_cell(k,iCell) &
                   + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) &
                   / areaCell(iCell)
              vorticity_cell_tl(k,iCell) = vorticity_cell_tl(k,iCell) &
                   + kiteAreasOnVertex(i, iVertex) * vorticity_tl(k, iVertex) &
                   / areaCell(iCell)
              vorticity_cell(k,iCell) = vorticity_cell(k,iCell) &
                   + kiteAreasOnVertex(i, iVertex) * vorticity(k, iVertex) &
                   / areaCell(iCell)
           enddo
        endif
     enddo
  enddo

  ! Compute gradient of PV in normal direction
  !   ( this computes gradPVn for all edges bounding real cells )
  gradPVn_tl(:,:) = 0.0
  gradPVn(:,:) = 0.0
  do iEdge = 1, nEdges
     if( cellsOnEdge(1,iEdge) <= nCells .and. cellsOnEdge(2,iEdge) <= nCells) then
        do k = 1, nVertLevels
           gradPVn_tl(k,iEdge) = (pv_cell_tl(k,cellsOnEdge(2,iEdge)) - pv_cell_tl(k,cellsOnEdge(1,iEdge))) / dcEdge(iEdge)
           gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / dcEdge(iEdge)
        enddo
     endif
  enddo

  ! Modify PV edge with upstream bias.
  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        pv_edge_tl(k,iEdge) = pv_edge_tl(k,iEdge) &
             - config_apvm_upwinding * u_tl(k,iEdge) * dt * gradPVn(k,iEdge) &
             - config_apvm_upwinding * u(k,iEdge) * dt * gradPVn_tl(k,iEdge)
        pv_edge(k,iEdge) = pv_edge(k,iEdge) &
             - config_apvm_upwinding * u(k,iEdge) * dt * gradPVn(k,iEdge)
     enddo
  enddo
  
end subroutine sw_compute_solve_diagnostics_tlm

subroutine sw_compute_tend_tlm(nCells, nEdges, nVertices, nVertLevels, &
     maxEdges, maxEdges2, vertexDegree, &
     config_bottom_drag, config_wind_stress, config_h_mom_eddy_visc2, &
     config_h_mom_eddy_visc4, h, h_tl, u, u_tl, v, v_tl, h_edge, h_edge_tl, &
     circulation, circulation_tl, vorticity, vorticity_tl, &
     divergence, divergence_tl, ke, ke_tl, pv_edge, pv_edge_tl, vh, vh_tl, &
     weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, cellsOnVertex, &
     verticesOnEdge, nEdgesOnCell, edgesOnCell, nEdgesOnEdge, &
     edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, areaTriangle, &
     h_s, fVertex, fEdge, u_src, meshScalingDel2, meshScalingDel4, &
     tend_h, tend_h_tl, tend_u, tend_u_tl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute height and normal wind tendencies, as well as diagnostic variables
  ! Input: s - current model state
  !        grid - grid metadata
  ! Output: tend - computed tendencies for prognostic variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)
  real (kind=RKIND), parameter :: gravity = 9.80616
  integer :: nCells, nEdges, nVertices, nVertLevels
  integer :: maxEdges, maxEdges2, vertexDegree
  real(kind=RKIND), dimension(nCells) :: h_s
  real(kind=RKIND), dimension(nVertices) :: fVertex
  real(kind=RKIND), dimension(nEdges) :: fEdge
  real(kind=RKIND), dimension(nEdges) :: dvEdge
  real(kind=RKIND), dimension(nEdges) :: dcEdge
  real(kind=RKIND), dimension(nCells) :: areaCell
  real(kind=RKIND), dimension(nVertices) :: areaTriangle
  real(kind=RKIND), dimension(nEdges) :: meshScalingDel2
  real(kind=RKIND), dimension(nEdges) :: meshScalingDel4
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: vh
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: vh_tl
  real(kind=RKIND), dimension(maxEdges2,nEdges) :: weightsOnEdge
  real(kind=RKIND), dimension(vertexDegree,nVertices) :: kiteAreasOnVertex
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: h
  real(kind=RKIND), dimension(nVertLevels,nCells) :: h_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: u
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: u_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: tend_h
  real(kind=RKIND), dimension(nVertLevels,nCells) :: tend_h_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: tend_u
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: tend_u_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke_tl
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge_tl
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence_tl
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex_tl

  integer, dimension(2,nEdges) :: cellsOnEdge
  integer, dimension(vertexDegree,nVertices) :: cellsOnVertex
  integer, dimension(2,nEdges) :: verticesOnEdge
  integer, dimension(maxEdges,nCells) :: edgesOnCell
  integer, dimension(maxEdges2,nEdges) :: edgesOnEdge
  integer, dimension(vertexDegree,nVertices) :: edgesOnVertex
  integer, dimension(nCells) :: nEdgesOnCell
  integer, dimension(nEdges) :: nEdgesOnEdge

  real (kind=RKIND), dimension(nVertLevels, nEdges) :: u_src

  logical :: config_wind_stress, config_bottom_drag
  real (kind=RKIND) :: config_h_mom_eddy_visc2, config_h_mom_eddy_visc4

  ! ----- local variables -----
  integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, timeLevel
  real (kind=RKIND) :: flux, vorticity_abs, workpv, q, upstream_bias
  real (kind=RKIND) :: flux_tl, vorticity_abs_tl, workpv_tl, q_tl, upstream_bias_tl
  real (kind=RKIND) :: ke_edge, ke_edge_tl
  real (kind=RKIND), parameter :: rho_ref = 1000.0
  real (kind=RKIND) :: r, u_diffusion, u_diffusion_tl
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence_tl
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u_tl
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation_tl, delsq_vorticity_tl
  ! ----- end local variables

  ! Compute height tendency for each cell
  tend_h_tl(:,:) = 0.0
  tend_h(:,:) = 0.0
  do iEdge = 1, nEdges
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     do k = 1, nVertLevels
        flux_tl = u_tl(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge) &
             + u(k,iEdge) * dvEdge(iEdge) * h_edge_tl(k,iEdge)
        flux = u(k,iEdge) * dvEdge(iEdge) * h_edge(k,iEdge)
        tend_h_tl(k,cell1) = tend_h_tl(k,cell1) - flux_tl
        tend_h(k,cell1) = tend_h(k,cell1) - flux
        tend_h_tl(k,cell2) = tend_h_tl(k,cell2) + flux_tl
        tend_h(k,cell2) = tend_h(k,cell2) + flux
     end do
  end do
  do iCell = 1, nCells
     do k = 1, nVertLevels
        tend_h_tl(k,iCell) = tend_h_tl(k,iCell) / areaCell(iCell)
        tend_h(k,iCell) = tend_h(k,iCell) / areaCell(iCell)
     end do
  end do

  ! Compute u (normal) velocity tendency for each edge (cell face)
  tend_u_tl(:,:) = 0.0
  tend_u(:,:) = 0.0
  do iEdge = 1, nEdges
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     vertex1 = verticesOnEdge(1,iEdge)
     vertex2 = verticesOnEdge(2,iEdge)

     do k = 1, nVertLevels
        q_tl=0.
        q = 0.0
        do j = 1, nEdgesOnEdge(iEdge)
           eoe = edgesOnEdge(j,iEdge)
           workpv_tl = 0.5 * (pv_edge_tl(k,iEdge) + pv_edge_tl(k,eoe))
           workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))

           q_tl = q_tl &
                + weightsOnEdge(j,iEdge) * u_tl(k,eoe) * workpv * h_edge(k,eoe) &
                + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv_tl * h_edge(k,eoe) &
                + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * h_edge_tl(k,eoe)
           q = q + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * h_edge(k,eoe)
        end do

        tend_u_tl(k,iEdge) = q_tl     &
             - (   ke_tl(k,cell2) - ke_tl(k,cell1) + &
             gravity * (h_tl(k,cell2) - h_tl(k,cell1))  ) / dcEdge(iEdge)
        tend_u(k,iEdge) = q     &
             - (   ke(k,cell2) - ke(k,cell1) + &
             gravity * (h(k,cell2) + h_s(cell2) - h(k,cell1) - h_s(cell1)) &
             ) / dcEdge(iEdge)
     end do
  end do


  ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
  !                    only valid for visc == constant
  if (config_h_mom_eddy_visc2 > 0.0) then
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)

        do k = 1, nVertLevels
           u_diffusion_tl =   ( divergence_tl(k,cell2) - divergence_tl(k,cell1) ) / dcEdge(iEdge) &
                -(vorticity_tl(k,vertex2) - vorticity_tl(k,vertex1) ) / dvEdge(iEdge)
           u_diffusion =   ( divergence(k,cell2)  -  divergence(k,cell1) ) / dcEdge(iEdge) &
                -(vorticity(k,vertex2)  - vorticity(k,vertex1) ) / dvEdge(iEdge)
           u_diffusion_tl = meshScalingDel2(iEdge) * config_h_mom_eddy_visc2 * u_diffusion_tl
           u_diffusion = meshScalingDel2(iEdge) * config_h_mom_eddy_visc2 * u_diffusion
           tend_u_tl(k,iEdge) = tend_u_tl(k,iEdge) + u_diffusion_tl
           tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion
        end do
     end do
  end if

  ! velocity tendency: del4 dissipation, -\nu_4 \nabla^4 u
  !   computed as \nabla^2 u = \nabla divergence + k \times \nabla vorticity
  !   applied recursively.
  !   strictly only valid for h_mom_eddy_visc4 == constant
  if (config_h_mom_eddy_visc4 > 0.0) then
     allocate(delsq_divergence_tl(nVertLevels, nCells+1))
     allocate(delsq_divergence(nVertLevels, nCells+1))
     allocate(delsq_u_tl(nVertLevels, nEdges+1))
     allocate(delsq_u(nVertLevels, nEdges+1))
     allocate(delsq_circulation_tl(nVertLevels, nVertices+1))
     allocate(delsq_circulation(nVertLevels, nVertices+1))
     allocate(delsq_vorticity_tl(nVertLevels, nVertices+1))
     allocate(delsq_vorticity(nVertLevels, nVertices+1))

     delsq_u_tl(:,:) = 0.0
     delsq_u(:,:) = 0.0

     ! Compute \nabla^2 u = \nabla divergence + k \times \nabla vorticity
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)

        do k = 1, nVertLevels

           delsq_u_tl(k,iEdge) = ( divergence_tl(k,cell2)-divergence_tl(k,cell1) ) /dcEdge(iEdge)  &
                -( vorticity_tl(k,vertex2) - vorticity_tl(k,vertex1)) / dvEdge(iEdge)
           delsq_u(k,iEdge) = ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                -( vorticity(k,vertex2) - vorticity(k,vertex1)) / dvEdge(iEdge)

        end do
     end do

     ! vorticity using \nabla^2 u
     delsq_circulation_tl(:,:) = 0.0
     delsq_circulation(:,:) = 0.0
     do iEdge = 1, nEdges
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)
        do k=1,nVertLevels
           delsq_circulation_tl(k,vertex1) = delsq_circulation_tl(k,vertex1) &
                - dcEdge(iEdge) * delsq_u_tl(k,iEdge)
           delsq_circulation(k,vertex1) = delsq_circulation(k,vertex1) &
                - dcEdge(iEdge) * delsq_u(k,iEdge)
           delsq_circulation_tl(k,vertex2) = delsq_circulation_tl(k,vertex2) &
                + dcEdge(iEdge) * delsq_u_tl(k,iEdge)
           delsq_circulation(k,vertex2) = delsq_circulation(k,vertex2) &
                + dcEdge(iEdge) * delsq_u(k,iEdge)
        end do
     end do
     do iVertex = 1, nVertices
        r = 1.0 / areaTriangle(iVertex)
        do k = 1, nVertLevels
           delsq_vorticity_tl(k,iVertex) = delsq_circulation_tl(k,iVertex) * r
           delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
        end do
     end do

     ! Divergence using \nabla^2 u
     delsq_divergence_tl(:,:) = 0.0
     delsq_divergence(:,:) = 0.0
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        do k=1,nVertLevels
           delsq_divergence_tl(k,cell1) = delsq_divergence_tl(k,cell1) &
                + delsq_u_tl(k,iEdge)*dvEdge(iEdge)
           delsq_divergence(k,cell1) = delsq_divergence(k,cell1) &
                + delsq_u(k,iEdge)*dvEdge(iEdge)
           delsq_divergence_tl(k,cell2) = delsq_divergence_tl(k,cell2) &
                - delsq_u_tl(k,iEdge)*dvEdge(iEdge)
           delsq_divergence(k,cell2) = delsq_divergence(k,cell2) &
                - delsq_u(k,iEdge)*dvEdge(iEdge)
        end do
     end do
     do iCell = 1,nCells
        r = 1.0 / areaCell(iCell)
        do k = 1, nVertLevels
           delsq_divergence_tl(k,iCell) = delsq_divergence_tl(k,iCell) * r
           delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
        end do
     end do

     ! Compute - \kappa \nabla^4 u
     ! as  \nabla div(\nabla^2 u) + k \times \nabla ( k \cross curl(\nabla^2 u) )
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)

        do k = 1, nVertLevels

           u_diffusion_tl = (  delsq_divergence_tl(k,cell2) &
                - delsq_divergence_tl(k,cell1) ) / dcEdge(iEdge)  &
                -(  delsq_vorticity_tl(k,vertex2) &
                - delsq_vorticity_tl(k,vertex1) ) / dvEdge(iEdge)
           u_diffusion = (  delsq_divergence(k,cell2) &
                - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                -(  delsq_vorticity(k,vertex2) &
                - delsq_vorticity(k,vertex1) ) / dvEdge(iEdge)

           u_diffusion_tl = meshScalingDel4(iEdge) * config_h_mom_eddy_visc4 * u_diffusion_tl
           u_diffusion = meshScalingDel4(iEdge) * config_h_mom_eddy_visc4 * u_diffusion

           tend_u_tl(k,iEdge) = tend_u_tl(k,iEdge) - u_diffusion_tl
           tend_u(k,iEdge) = tend_u(k,iEdge) - u_diffusion
        end do
     end do

     deallocate(delsq_divergence_tl)
     deallocate(delsq_divergence)
     deallocate(delsq_u_tl)
     deallocate(delsq_u)
     deallocate(delsq_circulation_tl)
     deallocate(delsq_circulation)
     deallocate(delsq_vorticity_tl)
     deallocate(delsq_vorticity)

  end if

  ! Compute u (velocity) tendency from wind stress (u_src)
  if(config_wind_stress) then
     do iEdge = 1, nEdges
        tend_u_tl(1,iEdge) =  tend_u_tl(1,iEdge) &
             - h_edge_tl(1,iEdge)*u_src(1,iEdge) / rho_ref / (h_edge(1,iEdge)**2)
        tend_u(1,iEdge) =  tend_u(1,iEdge) &
             + u_src(1,iEdge) / rho_ref / h_edge(1,iEdge)
     end do
  endif

  if (config_bottom_drag) then
     do iEdge = 1,  nEdges
        ! bottom drag is the same as POP:
        ! -c |u| u  where c is unitless and 1.0e-3.
        ! see POP Reference guide, section 3.4.4.
        ke_edge_tl = 0.5 * ( ke_tl(1,cellsOnEdge(1,iEdge)) &
             + ke_tl(1,cellsOnEdge(2,iEdge)))
        ke_edge = 0.5 * ( ke(1,cellsOnEdge(1,iEdge)) &
             + ke(1,cellsOnEdge(2,iEdge)))
        tend_u_tl(1,iEdge) = tend_u_tl(1,iEdge) &
             - 1.0e-3*u_tl(1,iEdge)*sqrt(2.0*ke_edge)/h_edge(1,iEdge) &
             - 0.5*1.0e-3*u(1,iEdge)*(2.0*ke_edge_tl)/sqrt(2.0*ke_edge)/h_edge(1,iEdge) &
             + 1.0e-3*u(1,iEdge)*sqrt(2.0*ke_edge)*h_edge_tl(1,iEdge)/(h_edge(1,iEdge)**2)
        tend_u(1,iEdge) = tend_u(1,iEdge) &
             - 1.0e-3*u(1,iEdge)*sqrt(2.0*ke_edge)/h_edge(1,iEdge)
     end do
  endif

end subroutine sw_compute_tend_tlm

subroutine mpas_reconstruct_2d_tlm(nCells, nEdges, maxEdges, nVertLevels, R3, &
     on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
     coeffs_reconstruct, u, u_tl, uReconstructX, uReconstructX_tl, &
     uReconstructY, uReconstructY_tl, uReconstructZ, uReconstructZ_tl, &
     uReconstructZonal, uReconstructZonal_tl, &
     uReconstructMeridional, uReconstructMeridional_tl)
  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)

  integer :: nCells, nEdges, nVertLevels, maxEdges, R3
  real (kind=RKIND), dimension(nVertLevels,nEdges) :: u, u_tl
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructX
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructX_tl
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructY
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructY_tl
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZ
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZ_tl
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZonal
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZonal_tl
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructMeridional
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructMeridional_tl
  !logical, optional, intent(in) :: includeHalos

  integer, dimension(maxEdges,nCells) :: edgesOnCell
  integer, dimension(nCells) :: nEdgesOnCell
  real(kind=RKIND), dimension(nCells) :: latCell, lonCell
  real (kind=RKIND), dimension(R3,maxEdges,nCells) :: coeffs_reconstruct
  logical :: on_a_sphere

  !f2py intent(in) nCells, nEdges, nVertLevels, maxEdges, R3, on_a_sphere
  !f2py intent(in) edgesOnCell, nEdgesOnCell, latCell, lonCell,
  !f2py intent(in) coeffs_reconstruct, u, u_tl
  !f2py intent(out) uReconstructX, uReconstructX_tl
  !f2py intent(out) uReconstructY, uReconstructY_tl
  !f2py intent(out) uReconstructZ, uReconstructZ_tl
  !f2py intent(out) uReconstructZonal, uReconstructZonal_tl
  !f2py intent(out) uReconstructMeridional, uReconstructMeridional_tl

  ! local variable
  logical :: includeHalosLocal
  integer :: iCell,iEdge, i
  real (kind=RKIND) :: clat, slat, clon, slon

  do iCell = 1, nCells
     uReconstructX_tl(:,iCell) = 0.0_RKIND
     uReconstructX(:,iCell) = 0.0_RKIND
     uReconstructY_tl(:,iCell) = 0.0_RKIND
     uReconstructY(:,iCell) = 0.0_RKIND
     uReconstructZ_tl(:,iCell) = 0.0_RKIND
     uReconstructZ(:,iCell) = 0.0_RKIND

     do i=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        uReconstructX_tl(:,iCell) = uReconstructX_tl(:,iCell) &
             + coeffs_reconstruct(1,i,iCell) * u_tl(:,iEdge)
        uReconstructX(:,iCell) = uReconstructX(:,iCell) &
             + coeffs_reconstruct(1,i,iCell) * u(:,iEdge)

        uReconstructY_tl(:,iCell) = uReconstructY_tl(:,iCell) &
             + coeffs_reconstruct(2,i,iCell) * u_tl(:,iEdge)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) &
             + coeffs_reconstruct(2,i,iCell) * u(:,iEdge)

        uReconstructZ_tl(:,iCell) = uReconstructZ_tl(:,iCell) &
             + coeffs_reconstruct(3,i,iCell) * u_tl(:,iEdge)
        uReconstructZ(:,iCell) = uReconstructZ(:,iCell) &
             + coeffs_reconstruct(3,i,iCell) * u(:,iEdge)
     enddo
  enddo

  if (on_a_sphere) then
     do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))
        uReconstructZonal_tl(:,iCell) = -uReconstructX_tl(:,iCell)*slon + &
             uReconstructY_tl(:,iCell)*clon
        uReconstructZonal(:,iCell) = -uReconstructX(:,iCell)*slon + &
             uReconstructY(:,iCell)*clon

        uReconstructMeridional_tl(:,iCell) = -(uReconstructX_tl(:,iCell)*clon &
             + uReconstructY_tl(:,iCell)*slon)*slat &
             + uReconstructZ_tl(:,iCell)*clat
        uReconstructMeridional(:,iCell) = -(uReconstructX(:,iCell)*clon       &
             + uReconstructY(:,iCell)*slon)*slat &
             + uReconstructZ(:,iCell)*clat
     end do
  else
     do iCell = 1, nCells
        uReconstructZonal_tl     (:,iCell) = uReconstructX_tl(:,iCell)
        uReconstructZonal     (:,iCell) = uReconstructX(:,iCell)
        uReconstructMeridional_tl(:,iCell) = uReconstructY_tl(:,iCell)
        uReconstructMeridional(:,iCell) = uReconstructY(:,iCell)
     end do
  end if

end subroutine mpas_reconstruct_2d_tlm
