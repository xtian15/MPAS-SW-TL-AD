subroutine sw_compute_solve_diagnostics_adj(dt, nCells, nEdges, nVertices, &
     maxEdges, maxEdges2, nVertLevels, vertexDegree, &
     config_monotonic, config_thickness_adv_order, config_apvm_upwinding, &
     h, h_ad, u, u_ad, v_in, v_ad, vh, vh_ad, h_edge, h_edge_ad, &
     h_vertex_in, h_vertex_ad, circulation_in, circulation_ad, &
     vorticity_in, vorticity_ad, divergence, divergence_ad, ke, ke_ad, &
     pv_edge, pv_edge_ad, pv_vertex_in, pv_vertex_ad, pv_cell_in, pv_cell_ad, &
     vorticity_cell, vorticity_cell_ad,  gradPVn_in, gradPVn_ad, &
     gradPVt_in, gradPVt_ad, weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, &
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


  ! ----- forward calculation requirement -----
  !needs updated values of h and u
  ! ----- end -----
  
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
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: vh_ad
  real(kind=RKIND), dimension(maxEdges2,nEdges) :: weightsOnEdge
  real(kind=RKIND), dimension(vertexDegree,nVertices) :: kiteAreasOnVertex
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge_ad
  real(kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: h
  real(kind=RKIND), dimension(nVertLevels,nCells) :: h_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: u
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: u_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v_in
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation_in
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity_in
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: pv_vertex_in
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: pv_vertex
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: pv_vertex_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: pv_cell_in
  real(kind=RKIND), dimension(nVertLevels,nCells) :: pv_cell
  real(kind=RKIND), dimension(nVertLevels,nCells) :: pv_cell_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVn_in
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVn
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVn_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVt_in
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVt
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: gradPVt_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex_in
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: vorticity_cell
  real(kind=RKIND), dimension(nVertLevels,nCells) :: vorticity_cell_ad

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
  real (kind=RKIND) :: d2fdx2_cell1_ad, d2fdx2_cell2_ad
  real (kind=RKIND) :: coef_3rd_order

  integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, cov
  real (kind=RKIND) :: flux, vorticity_abs, workpv
  real (kind=RKIND) :: flux_ad, vorticity_abs_ad, workpv_ad
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

  coef_3rd_order = 0.
  if (config_thickness_adv_order == 3) coef_3rd_order = 1.0
  if (config_thickness_adv_order == 3 .and. config_monotonic) coef_3rd_order = 0.25

  ! ----- forward calculation -----
  circulation(:,:) = 0.0
  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
        circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
     end do
  end do
  do iVertex = 1, nVertices
     do k = 1, nVertLevels
        vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
     end do
  end do

  v(:,:) = 0.0
  do iEdge = 1,nEdges
     do i = 1, nEdgesOnEdge(iEdge)
        eoe = edgesOnEdge(i,iEdge)
        do k = 1,nVertLevels
           v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
        end do
     end do
  end do

  do iVertex = 1,nVertices
     do k = 1, nVertLevels
        h_vertex(k,iVertex) = 0.0
        do i = 1, vertexDegree
           h_vertex(k,iVertex) = h_vertex(k,iVertex) + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
        end do
        h_vertex(k,iVertex) = h_vertex(k,iVertex) / areaTriangle(iVertex)
        pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex(k,iVertex)
     end do
  end do

  do iEdge = 1, nEdges
     do k = 1, nVertLevels
        gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / dvEdge(iEdge)
     enddo
  enddo

  pv_cell(:,:) = 0.0
  do iVertex = 1, nVertices
     do i = 1, vertexDegree
        iCell = cellsOnVertex(i,iVertex)
        if (iCell <= nCells) then
           do k = 1, nVertLevels
              pv_cell(k,iCell) = pv_cell(k,iCell) &
                   + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) &
                   / areaCell(iCell)
           end do
        end if
     end do
  end do

  gradPVn(:,:) = 0.0
  do iEdge = 1, nEdges
     if( cellsOnEdge(1,iEdge) <= nCells .and. cellsOnEdge(2,iEdge) <= nCells) then
        do k = 1, nVertLevels
           gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / dcEdge(iEdge)
        enddo
     endif
  enddo
  
  ! ----- end forward calculation -----

  ! Modify PV edge with upstream bias.
  do iEdge = nEdges, 1, -1
     do k = nVertLevels, 1, -1
        gradPVn_ad(k,iEdge)=gradPVn_ad(k,iEdge) &
             -config_apvm_upwinding * u(k,iEdge) * dt * pv_edge_ad(k,iEdge)
        u_ad(k,iEdge)=u_ad(k,iEdge) &
             -config_apvm_upwinding * dt * gradPVn(k,iEdge) * pv_edge_ad(k,iEdge)
     enddo
  enddo

  do iEdge = nEdges, 1, -1
     if( cellsOnEdge(1,iEdge) <= nCells .and. cellsOnEdge(2,iEdge) <= nCells) then
        do k = nVertLevels, 1, -1
           pv_cell_ad(k,cellsOnEdge(1,iEdge))=pv_cell_ad(k,cellsOnEdge(1,iEdge)) &
                -gradPVn_ad(k,iEdge)/dcEdge(iEdge)
           pv_cell_ad(k,cellsOnEdge(2,iEdge))=pv_cell_ad(k,cellsOnEdge(2,iEdge)) &
                +gradPVn_ad(k,iEdge)/dcEdge(iEdge)
           gradPVn_ad(k,iEdge)=0.
        enddo
     endif
  enddo
  gradPVn_ad(:,:) = 0.0

  do iVertex = nVertices, 1, -1
     do i = vertexDegree, 1, -1
        iCell = cellsOnVertex(i,iVertex)
        if (iCell <= nCells) then
           do k = nVertLevels, 1, -1
              vorticity_ad(k, iVertex)=vorticity_ad(k, iVertex) &
                   +kiteAreasOnVertex(i, iVertex)/areaCell(iCell) * vorticity_cell_ad(k,iCell)
              pv_vertex_ad(k, iVertex)=pv_vertex_ad(k, iVertex) &
                   +kiteAreasOnVertex(i, iVertex)/areaCell(iCell) * pv_cell_ad(k,iCell)
           enddo
        endif
     enddo
  enddo
  vorticity_cell_ad(:,:) = 0.0
  pv_cell_ad(:,:) = 0.0

  ! ----- forward calculation -----
  !gradPVt=???
  !v=???
  ! ----- end forward calculation -----
  do iEdge = nEdges, 1, -1
     do k = nVertLevels, 1, -1
        gradPVt_ad(k,iEdge)=gradPVt_ad(k,iEdge) &
             -config_apvm_upwinding * v(k,iEdge) * dt * pv_edge_ad(k,iEdge)
        v_ad(k,iEdge)=v_ad(k,iEdge) &
             -config_apvm_upwinding * dt * gradPVt(k,iEdge) * pv_edge_ad(k,iEdge)
     end do
  end do
  
  do iVertex = nVertices, 1, -1
     do i = vertexDegree, 1, -1
        iEdge = edgesOnVertex(i,iVertex)
        do k = nVertLevels, 1, -1
           pv_vertex_ad(k,iVertex)=pv_vertex_ad(k,iVertex) + 0.5*pv_edge_ad(k,iEdge)
        end do
     end do
  end do
  pv_edge_ad(:,:) = 0.0

  do iEdge = nEdges, 1, -1
     do k = nVertLevels, 1, -1
        pv_vertex_ad(k,verticesOnEdge(1,iEdge))=pv_vertex_ad(k,verticesOnEdge(1,iEdge)) &
             -gradPVt_ad(k,iEdge)/dvEdge(iEdge)
        pv_vertex_ad(k,verticesOnEdge(2,iEdge))=pv_vertex_ad(k,verticesOnEdge(2,iEdge)) &
             +gradPVt_ad(k,iEdge)/dvEdge(iEdge)
        gradPVt_ad(k,iEdge)=0.
     end do
  end do

  ! ----- forward calculation -----
  !vorticity=???
  !h_vertex=???
  ! ----- end forward calculation -----
  do iVertex = nVertices, 1, -1
     do k = nVertLevels, 1, -1
        ! ----- forward calculation -----
        !h_vertex(k,iVertex) = 0.0
        !do i = 1, vertexDegree
        !   h_vertex(k,iVertex) = h_vertex(k,iVertex) + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
        !end do
        !h_vertex(k,iVertex) = h_vertex(k,iVertex) / areaTriangle(iVertex)
        ! ----- end forward calculation -----
        h_vertex_ad(k,iVertex)=h_vertex_ad(k,iVertex) &
             -(fVertex(iVertex) + vorticity(k,iVertex))/(h_vertex(k,iVertex)**2)&
             *pv_vertex_ad(k,iVertex)
        vorticity_ad(k,iVertex)=vorticity_ad(k,iVertex) &
             +pv_vertex_ad(k,iVertex)/h_vertex(k,iVertex)
        pv_vertex_ad(k,iVertex)=0.

        h_vertex_ad(k,iVertex)=h_vertex_ad(k,iVertex)/areaTriangle(iVertex)
        
        do i = vertexDegree, 1, -1
           h_ad(k,cellsOnVertex(i,iVertex))=h_ad(k,cellsOnVertex(i,iVertex)) &
                +kiteAreasOnVertex(i,iVertex)*h_vertex_ad(k,iVertex)
        end do
        h_vertex_ad(k,iVertex) = 0.0
     end do
  end do

  ! vh omitted

  do iEdge = nEdges, 1, -1
     do i = nEdgesOnEdge(iEdge), 1, -1
        eoe = edgesOnEdge(i,iEdge)
        do k = nVertLevels, 1, -1
           u_ad(k, eoe)=u_ad(k, eoe) + weightsOnEdge(i,iEdge)*v_ad(k,iEdge)
        end do
     end do
  end do
  v_ad(:,:)=0.

  do iCell = nCells, 1, -1
     do k=nVertLevels, 1, -1
        ke_ad(k,iCell)=ke_ad(k,iCell)/areaCell(iCell)
     end do
     do i = nEdgesOnCell(iCell), 1, -1
        iEdge = edgesOnCell(i,iCell)
        do k = nVertLevels, 1, -1
           u_ad(k,iEdge)=u_ad(k,iEdge) &
                +0.25*dcEdge(iEdge)*dvEdge(iEdge)*2.*u(k,iEdge)*ke_ad(k,iCell)
        end do
     end do
  end do
  ke_ad(:,:)=0.


  do iCell = nCells, 1, -1
     r = 1.0 / areaCell(iCell)
     do k = nVertLevels, 1, -1
        divergence_ad(k,iCell)=divergence_ad(k,iCell)*r
     end do
  end do
  do iEdge = nEdges, 1, -1
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     if(cell2 <= nCells) then
        do k = nVertLevels, 1, -1
           u_ad(k,iEdge)=u_ad(k,iEdge) - dvEdge(iEdge)*divergence_ad(k,cell2)
        end do
     end if
     if(cell1 <= nCells) then
        do k = nVertLevels, 1, -1
           u_ad(k,iEdge)=u_ad(k,iEdge) + dvEdge(iEdge)*divergence_ad(k,cell1)
        end do
     end if
  end do
  divergence_ad(:,:) = 0.0

  do iVertex = nVertices, 1, -1
     do k = nVertLevels, 1, -1
        circulation_ad(k,iVertex)=circulation_ad(k,iVertex) &
          + vorticity_ad(k,iVertex)/areaTriangle(iVertex)
        vorticity_ad(k,iVertex)=0.
     end do
  end do
  do iEdge = nEdges, 1, -1
     do k = nVertLevels, 1, -1
        u_ad(k,iEdge)=u_ad(k,iEdge) + dcEdge(iEdge)*circulation_ad(k,verticesOnEdge(2,iEdge))
        u_ad(k,iEdge)=u_ad(k,iEdge) - dcEdge(iEdge)*circulation_ad(k,verticesOnEdge(1,iEdge))
     end do
  end do
  circulation_ad(:,:) = 0.0

  if (config_thickness_adv_order == 2) then
     do iEdge = nEdges, 1, -1
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        if (cell1 <= nCells .and. cell2 <= nCells) then
           do k = nVertLevels, 1, -1
              h_ad(k,cell2)=h_ad(k,cell2) + 0.5*h_edge_ad(k,iEdge)
              h_ad(k,cell1)=h_ad(k,cell1) + 0.5*h_edge_ad(k,iEdge)
              h_edge_ad(k,iEdge)=0.
           end do
        end if
     end do
  else if (config_thickness_adv_order == 3) then
     do iEdge = nEdges, 1, -1
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        if (cell1 <= nCells .and. cell2 <= nCells) then
           do k = 1, nVertLevels
              if (u(k,iEdge) > 0) then
                 d2fdx2_cell2_ad= (dcEdge(iEdge)**2)*coef_3rd_order/12.*h_edge_ad(k,iEdge)
                 d2fdx2_cell1_ad=-(dcEdge(iEdge)**2)*coef_3rd_order/12.*h_edge_ad(k,iEdge)
                 d2fdx2_cell2_ad= d2fdx2_cell2_ad &
                      -(dcEdge(iEdge)**2)/12.*h_edge_ad(k,iEdge)
                 d2fdx2_cell1_ad= d2fdx2_cell1_ad &
                      -(dcEdge(iEdge)**2)/12.*h_edge_ad(k,iEdge)
                 h_ad(k,cell2)=h_ad(k,cell2) + 0.5*h_edge_ad(k,iEdge)
                 h_ad(k,cell1)=h_ad(k,cell1) + 0.5*h_edge_ad(k,iEdge)
                 h_edge_ad(k,iEdge)=0.
              else
                 d2fdx2_cell2_ad=-(dcEdge(iEdge)**2)*coef_3rd_order/12.*h_edge_ad(k,iEdge)
                 d2fdx2_cell1_ad= (dcEdge(iEdge)**2)*coef_3rd_order/12.*h_edge_ad(k,iEdge)
                 d2fdx2_cell2_ad=d2fdx2_cell2_ad &
                      -(dcEdge(iEdge)**2)/12.*h_edge_ad(k,iEdge)
                 d2fdx2_cell1_ad=d2fdx2_cell1_ad &
                      -(dcEdge(iEdge)**2)/12.*h_edge_ad(k,iEdge)
                 h_ad(k,cell2)=h_ad(k,cell2) + 0.5*h_edge_ad(k,iEdge)
                 h_ad(k,cell1)=h_ad(k,cell1) + 0.5*h_edge_ad(k,iEdge)
                 h_edge_ad(k,iEdge)=0.
              end if

              if(boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then
                 do i=nEdgesOnCell(cell2), 1, -1
                    h_ad(k, cellsOnCell(i,cell2))=h_ad(k, cellsOnCell(i,cell2)) &
                         +deriv_two(i+1,2,iEdge)*d2fdx2_cell2_ad
                 end do
                 do i = nEdgesOnCell(cell1), 1, -1
                    h_ad(k, cellsOnCell(i,cell1))=h_ad(k, cellsOnCell(i,cell1)) &
                         +deriv_two(i+1,1,iEdge)*d2fdx2_cell1_ad
                 end do
                 h_ad(k,cell2)=h_ad(k,cell2) + deriv_two(1,2,iEdge)*d2fdx2_cell2_ad
                 d2fdx2_cell2_ad=0.
                 h_ad(k,cell1)=h_ad(k,cell1) + deriv_two(1,1,iEdge)*d2fdx2_cell1_ad
                 d2fdx2_cell1_ad=0.
              end if ! (boundaryCell(k,cell1) .eq. 0 .and.
              d2fdx2_cell2_ad = 0.0
              d2fdx2_cell1_ad = 0.0
           end do ! k = 1, nVertLevels
        end if ! (cell1 <= nCells .and. cell2 <= nCells) then
     end do ! iEdge = nEdges, 1, -1

  else if (config_thickness_adv_order == 4) then
     do iEdge = nEdges, 1, -1
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)

        if (cell1 <= nCells .and. cell2 <= nCells) then
           do k = nVertLevels, 1, -1
              d2fdx2_cell2_ad=-(dcEdge(iEdge) **2)/12. * h_edge_ad(k,iEdge)
              d2fdx2_cell1_ad=-(dcEdge(iEdge) **2)/12. * h_edge_ad(k,iEdge)
              h_ad(k,cell2)=h_ad(k,cell2) + 0.5*h_edge_ad(k,iEdge)
              h_ad(k,cell1)=h_ad(k,cell1) + 0.5*h_edge_ad(k,iEdge)
              h_edge_ad(k,iEdge)=0.

              if (boundaryCell(k,cell1) .eq. 0 .and. boundaryCell(k,cell2) .eq. 0) then
                 do i = nEdgesOnCell(cell2), 1, -1
                    h_ad(k, cellsOnCell(i,cell2))=h_ad(k, cellsOnCell(i,cell2)) &
                         +deriv_two(i+1,2,iEdge)*d2fdx2_cell2_ad
                 end do
                 do i = nEdgesOnCell(cell1), 1, -1
                    h_ad(k, cellsOnCell(i,cell1))=h_ad(k, cellsOnCell(i,cell1)) &
                         +deriv_two(i+1,1,iEdge)*d2fdx2_cell1_ad
                 end do
                 h_ad(k,cell2)=h_ad(k,cell2) + deriv_two(1,2,iEdge)*d2fdx2_cell2_ad
                 d2fdx2_cell2_ad=0.
                 h_ad(k,cell1)=h_ad(k,cell1) + deriv_two(1,1,iEdge)*d2fdx2_cell1_ad
                 d2fdx2_cell1_ad=0.
              end if ! (boundaryCell(k,cell1) .eq. 0 .and.
              d2fdx2_cell2_ad = 0.0
              d2fdx2_cell1_ad = 0.0
           end do ! k=nVertLevels, 1, -1
        end if ! (cell1 <= nCells .and. cell2 <= nCells) then
     end do ! iEdge = 1, nEdges
     
  end if ! (config_thickness_adv_order == 2) then
  
end subroutine sw_compute_solve_diagnostics_adj

subroutine sw_compute_tend_adj(nCells, nEdges, nVertices, nVertLevels, &
     maxEdges, maxEdges2, vertexDegree, &
     config_bottom_drag, config_wind_stress, config_h_mom_eddy_visc2, &
     config_h_mom_eddy_visc4, h, h_ad, u, u_ad, v, v_ad, h_edge, h_edge_ad, &
     circulation, circulation_ad, vorticity, vorticity_ad, &
     divergence, divergence_ad, ke, ke_ad, pv_edge, pv_edge_ad, vh, vh_ad, &
     weightsOnEdge, kiteAreasOnVertex, cellsOnEdge, cellsOnVertex, &
     verticesOnEdge, nEdgesOnCell, edgesOnCell, nEdgesOnEdge, &
     edgesOnEdge, edgesOnVertex, dcEdge, dvEdge, areaCell, areaTriangle, &
     h_s, fVertex, fEdge, u_src, meshScalingDel2, meshScalingDel4, &
     tend_h, tend_h_ad, tend_u, tend_u_ad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute height and normal wind tendencies, as well as diagnostic variables
  ! Input: s - current model state
  !        grid - grid metadata
  ! Output: tend - computed tendencies for prognostic variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ----- forward calculation requirement -----
  ! needs updated values of ke, h_edge, and pv_edge
  ! ----- end -----
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
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: vh_ad
  real(kind=RKIND), dimension(maxEdges2,nEdges) :: weightsOnEdge
  real(kind=RKIND), dimension(vertexDegree,nVertices) :: kiteAreasOnVertex
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: h_edge_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: h
  real(kind=RKIND), dimension(nVertLevels,nCells) :: h_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: u
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: u_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: v_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: tend_h
  real(kind=RKIND), dimension(nVertLevels,nCells) :: tend_h_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: tend_u
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: tend_u_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: circulation_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: vorticity_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke
  real(kind=RKIND), dimension(nVertLevels,nCells) :: ke_ad
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge
  real(kind=RKIND), dimension(nVertLevels,nEdges) :: pv_edge_ad
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence
  real(kind=RKIND), dimension(nVertLevels,nCells) :: divergence_ad
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex
  real(kind=RKIND), dimension(nVertLevels,nVertices) :: h_vertex_ad

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
  real (kind=RKIND) :: flux_ad, vorticity_abs_ad, workpv_ad, q_ad, upstream_bias_ad
  real (kind=RKIND) :: ke_edge, ke_edge_ad
  real (kind=RKIND), parameter :: rho_ref = 1000.0
  real (kind=RKIND) :: r, u_diffusion, u_diffusion_ad
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_divergence_ad
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u_ad
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity
  real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation_ad, delsq_vorticity_ad
  ! ----- end local variables

  if (config_bottom_drag) then
     do iEdge = nEdges, 1, -1
        ! ----- forward calculation -----
        ke_edge = 0.5 * ( ke(1,cellsOnEdge(1,iEdge)) + ke(1,cellsOnEdge(2,iEdge)))
        ! ----- end forward calculation -----
        
        h_edge_ad(1,iEdge)=h_edge_ad(1,iEdge) &
             +1.0e-3*u(1,iEdge)*sqrt(2.0*ke_edge)/(h_edge(1,iEdge)**2)*tend_u_ad(1,iEdge)
        ke_edge_ad=-0.5*1.0e-3*u(1,iEdge)*2.0/sqrt(2.0*ke_edge)/h_edge(1,iEdge)*tend_u_ad(1,iEdge)
        u_ad(1,iEdge)=u_ad(1,iEdge) &
             -1.0e-3*sqrt(2.0*ke_edge)/h_edge(1,iEdge)*tend_u_ad(1,iEdge)

        ke_ad(1,cellsOnEdge(2,iEdge))=ke_ad(1,cellsOnEdge(2,iEdge)) +0.5*ke_edge_ad
        ke_ad(1,cellsOnEdge(1,iEdge))=ke_ad(1,cellsOnEdge(1,iEdge)) +0.5*ke_edge_ad
        ke_edge_ad=0.
     end do
  end if

  if(config_wind_stress) then
     do iEdge = nEdges, 1, -1
        h_edge_ad(1,iEdge)=h_edge_ad(1,iEdge) &
             -u_src(1,iEdge) / rho_ref / (h_edge(1,iEdge)**2) * tend_u_ad(1,iEdge)
     end do
  end if

  if (config_h_mom_eddy_visc4 > 0.0) then
     allocate(delsq_divergence_ad(nVertLevels, nCells+1))
     allocate(delsq_divergence(nVertLevels, nCells+1))
     allocate(delsq_u_ad(nVertLevels, nEdges+1))
     allocate(delsq_u(nVertLevels, nEdges+1))
     allocate(delsq_circulation_ad(nVertLevels, nVertices+1))
     allocate(delsq_circulation(nVertLevels, nVertices+1))
     allocate(delsq_vorticity_ad(nVertLevels, nVertices+1))
     allocate(delsq_vorticity(nVertLevels, nVertices+1))
     delsq_divergence_ad =0.
     delsq_u_ad          =0.
     delsq_circulation_ad=0.
     delsq_vorticity_ad  =0.

     do iEdge = nEdges, 1, -1
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)

        do k=nVertLevels, 1, -1
           u_diffusion_ad=-tend_u_ad(k,iEdge)
           u_diffusion_ad=meshScalingDel4(iEdge)*config_h_mom_eddy_visc4 * u_diffusion_ad

           delsq_vorticity_ad(k,vertex1)=delsq_vorticity_ad(k,vertex1) &
                +u_diffusion_ad/dvEdge(iEdge)
           delsq_vorticity_ad(k,vertex2)=delsq_vorticity_ad(k,vertex2) &
                -u_diffusion_ad/dvEdge(iEdge)
           delsq_divergence_ad(k,cell1)=delsq_divergence_ad(k,cell1) &
                -u_diffusion_ad/dcEdge(iEdge)
           delsq_divergence_ad(k,cell2)=delsq_divergence_ad(k,cell2) &
                +u_diffusion_ad/dcEdge(iEdge)
           u_diffusion_ad=0.
        end do
     end do ! iEdge = nEdges, 1, -1

     do iCell = nCells, 1, -1
        r = 1.0 / areaCell(iCell)
        do k = nVertLevels, 1, -1
           delsq_divergence_ad(k,iCell)=delsq_divergence_ad(k,iCell)*r
        end do
     end do

     do iEdge = nEdges, 1, -1
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        do k=nVertLevels, 1, -1
           delsq_u_ad(k,iEdge)=delsq_u_ad(k,iEdge) -dvEdge(iEdge)*delsq_divergence_ad(k,cell2)
           delsq_u_ad(k,iEdge)=delsq_u_ad(k,iEdge) +dvEdge(iEdge)*delsq_divergence_ad(k,cell1)
        end do
     end do
     delsq_divergence_ad(:,:) = 0.0

     do iVertex = nVertices, 1, -1
        r = 1.0 / areaTriangle(iVertex)
        do k = nVertLevels, 1, -1
           delsq_circulation_ad(k,iVertex)=delsq_circulation_ad(k,iVertex) &
                +delsq_vorticity_ad(k,iVertex) * r
           delsq_vorticity_ad(k,iVertex)=0.
        end do
     end do

     do iEdge = nEdges, 1, -1
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)
        do k=nVertLevels, 1, -1
           delsq_u_ad(k,iEdge)=delsq_u_ad(k,iEdge) + dcEdge(iEdge)*delsq_circulation_ad(k,vertex2)
           delsq_u_ad(k,iEdge)=delsq_u_ad(k,iEdge) - dcEdge(iEdge)*delsq_circulation_ad(k,vertex1)
        end do
     end do
     delsq_circulation_ad(:,:) = 0.0

     do iEdge = nEdges, 1, -1
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)
        do k = nVertLevels, 1, -1
           vorticity_ad(k,vertex1)=vorticity_ad(k,vertex1)+delsq_u_ad(k,iEdge)/dvEdge(iEdge)
           vorticity_ad(k,vertex2)=vorticity_ad(k,vertex2)-delsq_u_ad(k,iEdge)/dvEdge(iEdge)

           divergence_ad(k,cell1)=divergence_ad(k,cell1)-delsq_u_ad(k,iEdge)/dcEdge(iEdge)
           divergence_ad(k,cell2)=divergence_ad(k,cell2)+delsq_u_ad(k,iEdge)/dcEdge(iEdge)
        end do
     end do
     delsq_u_ad(:,:) = 0.0

     deallocate(delsq_divergence_ad)
     deallocate(delsq_divergence)
     deallocate(delsq_u_ad)
     deallocate(delsq_u)
     deallocate(delsq_circulation_ad)
     deallocate(delsq_circulation)
     deallocate(delsq_vorticity_ad)
     deallocate(delsq_vorticity)
  end if ! (config_h_mom_eddy_visc4 > 0.0)

  if (config_h_mom_eddy_visc2 > 0.0) then
     do iEdge = 1, nEdges
        cell1 = cellsOnEdge(1,iEdge)
        cell2 = cellsOnEdge(2,iEdge)
        vertex1 = verticesOnEdge(1,iEdge)
        vertex2 = verticesOnEdge(2,iEdge)
        do k = 1, nVertLevels
           u_diffusion_ad=tend_u_ad(k,iEdge)
           u_diffusion_ad=meshScalingDel2(iEdge)*config_h_mom_eddy_visc2*u_diffusion_ad
           vorticity_ad(k,vertex1)=vorticity_ad(k,vertex1) +u_diffusion_ad/dvEdge(iEdge)
           vorticity_ad(k,vertex2)=vorticity_ad(k,vertex2) -u_diffusion_ad/dvEdge(iEdge)
           divergence_ad(k,cell1)=divergence_ad(k,cell1) -u_diffusion_ad/dcEdge(iEdge)
           divergence_ad(k,cell2)=divergence_ad(k,cell2) +u_diffusion_ad/dcEdge(iEdge)
           u_diffusion_ad=0.
        end do
     end do ! iEdge = 1, nEdges
  end if ! (config_h_mom_eddy_visc2 > 0.0) then

  do iEdge = nEdges, 1, -1
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     vertex1 = verticesOnEdge(1,iEdge)
     vertex2 = verticesOnEdge(2,iEdge)

     do k = nVertLevels, 1, -1
        h_ad(k,cell1)=h_ad(k,cell1) +gravity/dcEdge(iEdge) * tend_u_ad(k,iEdge)
        h_ad(k,cell2)=h_ad(k,cell2) -gravity/dcEdge(iEdge) * tend_u_ad(k,iEdge)
        ke_ad(k,cell1)=ke_ad(k,cell1) + tend_u_ad(k,iEdge)/dcEdge(iEdge)
        ke_ad(k,cell2)=ke_ad(k,cell2) - tend_u_ad(k,iEdge)/dcEdge(iEdge)
        q_ad=tend_u_ad(k,iEdge)
        tend_u_ad(k,iEdge)=0.

        do j = nEdgesOnEdge(iEdge), 1, -1
           eoe = edgesOnEdge(j,iEdge)
           ! ----- forward calculation -----
           workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
           ! ----- end forward calculation -----
           h_edge_ad(k,eoe)=h_edge_ad(k,eoe) + weightsOnEdge(j,iEdge)*u(k,eoe)*workpv*q_ad
           workpv_ad=weightsOnEdge(j,iEdge)*u(k,eoe)*h_edge(k,eoe) * q_ad
           u_ad(k,eoe)=u_ad(k,eoe) + weightsOnEdge(j,iEdge)*workpv*h_edge(k,eoe) * q_ad
           pv_edge_ad(k,eoe)=pv_edge_ad(k,eoe) + 0.5*workpv_ad
           pv_edge_ad(k,iEdge)=pv_edge_ad(k,iEdge) + 0.5*workpv_ad
           workpv_ad=0.
        end do
        q_ad=0.
     end do
  end do
  tend_u_ad(:,:) = 0.0

  do iCell = nCells, 1, -1
     do k = nVertLevels, 1, -1
        tend_h_ad(k,iCell)=tend_h_ad(k,iCell)/areaCell(iCell)
     end do
  end do
  do iEdge = 1, nEdges
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)
     do k = 1, nVertLevels
        flux_ad=tend_h_ad(k,cell2)
        flux_ad=flux_ad-tend_h_ad(k,cell1)
        h_edge_ad(k,iEdge)=h_edge_ad(k,iEdge) +u(k,iEdge)*dvEdge(iEdge)*flux_ad
        u_ad(k,iEdge)=u_ad(k,iEdge) +dvEdge(iEdge)*h_edge(k,iEdge)*flux_ad
        flux_ad=0.
     end do
  end do
  tend_h_ad(:,:) = 0.0

end subroutine sw_compute_tend_adj

subroutine mpas_reconstruct_2d_adj(nCells, nEdges, maxEdges, nVertLevels, R3, &
     on_a_sphere, edgesOnCell, nEdgesOnCell, latCell, lonCell, &
     coeffs_reconstruct, u, u_ad, uReconstructX, uReconstructX_ad, &
     uReconstructY, uReconstructY_ad, uReconstructZ, uReconstructZ_ad, &
     uReconstructZonal, uReconstructZonal_ad, &
     uReconstructMeridional, uReconstructMeridional_ad)
  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)

  integer :: nCells, nEdges, nVertLevels, maxEdges, R3
  real (kind=RKIND), dimension(nVertLevels,nEdges), intent(in) :: u
  real (kind=RKIND), dimension(nVertLevels,nEdges) :: u_ad
  real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: uReconstructX
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructX_ad
  real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: uReconstructY
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructY_ad
  real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: uReconstructZ
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZ_ad
  real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: uReconstructZonal
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZonal_ad
  real (kind=RKIND), dimension(nVertLevels,nCells), intent(in) :: uReconstructMeridional
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructMeridional_ad
  !logical, optional, intent(in) :: includeHalos

  integer, dimension(maxEdges,nCells) :: edgesOnCell
  integer, dimension(nCells) :: nEdgesOnCell
  real(kind=RKIND), dimension(nCells) :: latCell, lonCell
  real (kind=RKIND), dimension(R3,maxEdges,nCells) :: coeffs_reconstruct
  logical :: on_a_sphere

  !f2py intent(in) nCells, nEdges, nVertLevels, maxEdges, R3, on_a_sphere,
  !f2py intent(in) edgesOnCell, nEdgesOnCell, latCell, lonCell,
  !f2py intent(in) coeffs_reconstruct, u, uReconstructX, uReconstructY,
  !f2py intent(in) uReconstructZ, uReconstructZonal, uReconstructMeridional

  !f2py intent(in,out) u_ad, uReconstructX_ad, uReconstructY_ad, uReconstructZ_ad,
  !f2py intent(in,out) uReconstructZonal_ad, uReconstructMeridional_ad

  ! local variable
  logical :: includeHalosLocal
  integer :: iCell,iEdge, i
  real (kind=RKIND) :: clat, slat, clon, slon

  if (on_a_sphere) then
     do iCell = nCells, 1, -1
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))

        uReconstructZ_ad(:,iCell)=uReconstructZ_ad(:,iCell) &
             +clat*uReconstructMeridional_ad(:,iCell)
        uReconstructY_ad(:,iCell)=uReconstructY_ad(:,iCell) &
             -slon*slat*uReconstructMeridional_ad(:,iCell)
        uReconstructX_ad(:,iCell)=uReconstructX_ad(:,iCell) &
             -clon*slat*uReconstructMeridional_ad(:,iCell)
        uReconstructMeridional_ad(:,iCell)=0._RKIND

        uReconstructY_ad(:,iCell)=uReconstructY_ad(:,iCell) &
             +clon*uReconstructZonal_ad(:,iCell)
        uReconstructX_ad(:,iCell)=uReconstructX_ad(:,iCell) &
             -slon*uReconstructZonal_ad(:,iCell)
        uReconstructZonal_ad(:,iCell)=0._RKIND
     end do
  else
     do iCell = nCells, 1, -1
        uReconstructY_ad(:,iCell)=uReconstructY_ad(:,iCell) &
             +uReconstructMeridional_ad(:,iCell)
        uReconstructMeridional_ad(:,iCell)=0._RKIND
        uReconstructX_ad(:,iCell)=uReconstructX_ad(:,iCell) &
             +uReconstructZonal_ad(:,iCell)
        uReconstructZonal_ad(:,iCell)=0._RKIND
     end do
  end if

  do iCell = nCells, 1, -1
     do i=nEdgesOnCell(iCell), 1, -1
        iEdge = edgesOnCell(i,iCell)
        u_ad(:,iEdge)=u_ad(:,iEdge) &
             + coeffs_reconstruct(3,i,iCell)*uReconstructZ_ad(:,iCell)
        u_ad(:,iEdge)=u_ad(:,iEdge) &
             + coeffs_reconstruct(2,i,iCell)*uReconstructY_ad(:,iCell)
        u_ad(:,iEdge)=u_ad(:,iEdge) &
             + coeffs_reconstruct(1,i,iCell)*uReconstructX_ad(:,iCell)
     end do
     uReconstructX_ad(:,iCell) = 0.0_RKIND
     uReconstructY_ad(:,iCell) = 0.0_RKIND
     uReconstructZ_ad(:,iCell) = 0.0_RKIND
  end do

end subroutine mpas_reconstruct_2d_adj
