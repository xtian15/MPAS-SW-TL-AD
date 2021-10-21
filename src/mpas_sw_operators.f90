!subroutine mpas_initialize_vectors(meshPool)!{{{
subroutine mpas_initialize_vectors(nCells, nEdges, maxEdges, R3, &
     verticesOnEdge, cellsOnEdge, &
     edgesOnCell, xCell, yCell, zCell, xEdge, yEdge, zEdge, &
     localVerticalUnitVectors, edgeNormalVectors, cellTangentPlane, &
     on_a_sphere, is_periodic, x_period, y_period)
  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)
  
  integer :: nCells, nEdges, maxEdges, R3
  integer, dimension(2, nEdges) :: verticesOnEdge, cellsOnEdge
  integer, dimension(maxEdges, nCells) :: edgesOnCell

  real(kind=RKIND), dimension(nCells) :: xCell, yCell, zCell
  real(kind=RKIND), dimension(nEdges) :: xEdge, yEdge, zEdge
  real(kind=RKIND), dimension(R3,nCells) :: localVerticalUnitVectors
  real(kind=RKIND), dimension(R3,nEdges) :: edgeNormalVectors
  real(kind=RKIND), dimension(R3,2,nCells) :: cellTangentPlane

  logical :: on_a_sphere, is_periodic
  real(kind=RKIND) :: x_period, y_period

  !f2py intent(in) nCells, nEdges, maxEdges, R3
  !f2py intent(in) verticesOnEdge, cellsOnEdge, edgesOnCell
  !f2py intent(in) xCell, yCell, zCell, xEdge, yEdge, zEdge
  !f2py intent(in) on_a_sphere, is_periodic, x_period, y_period
  !f2py intent(out) localVerticalUnitVectors, edgeNormalVectors, cellTangentPlane

  ! local variables
  integer :: iEdge, iCell, cell1, cell2
  real(kind=RKIND) :: mpas_fix_periodicity
  real(kind=RKIND), dimension(3) :: xHatPlane, yHatPlane, rHat
  real(kind=RKIND) :: normalDotRHat

  ! init arrays
  edgeNormalVectors = 0
  localVerticalUnitVectors = 0

  ! loop over all cells to be solved on this block
  do iCell = 1, nCells
     if (on_a_sphere) then
        localVerticalUnitVectors(1,iCell) = xCell(iCell)
        localVerticalUnitVectors(2,iCell) = yCell(iCell)
        localVerticalUnitVectors(3,iCell) = zCell(iCell)
        call mpas_unit_vec_in_r3(localVerticalUnitVectors(:,iCell))
     else ! on a plane
        localVerticalUnitVectors(:,iCell) = (/ 0., 0., 1. /)
     end if
  end do

  ! Initialize normal unit vectors at each edge
  ! These vectors point from cell to cell.
  ! At boundaries, one cell does not exist, so it points from cell to edge or from edge to cell.

  do iEdge = 1,nEdges
     cell1 = cellsOnEdge(1,iEdge)
     cell2 = cellsOnEdge(2,iEdge)

     if (cell1 == nCells+1) then ! this is a boundary edge
        ! the normal points from the edge location to the cell location
        if (is_periodic) then
           edgeNormalVectors(1,iEdge) = xCell(cell2) - mpas_fix_periodicity(xEdge(iEdge), xCell(cell2), x_period)
           edgeNormalVectors(2,iEdge) = yCell(cell2) - mpas_fix_periodicity(yEdge(iEdge), yCell(cell2), y_period)
           edgeNormalVectors(3,iEdge) = zCell(cell2) - zEdge(iEdge)
        else
           edgeNormalVectors(1,iEdge) = xCell(cell2) - xEdge(iEdge)
           edgeNormalVectors(2,iEdge) = yCell(cell2) - yEdge(iEdge)
           edgeNormalVectors(3,iEdge) = zCell(cell2) - zEdge(iEdge)
        end if

     else if (cell2 == nCells+1) then ! this is a boundary edge
        ! the normal points from the cell location to the edge location
        if (is_periodic) then
           edgeNormalVectors(1,iEdge) = mpas_fix_periodicity(xEdge(iEdge), xCell(cell1), x_period) - xCell(cell1)
           edgeNormalVectors(2,iEdge) = mpas_fix_periodicity(yEdge(iEdge), yCell(cell1), y_period) - yCell(cell1)
           edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)
        else
           edgeNormalVectors(1,iEdge) = xEdge(iEdge) - xCell(cell1)
           edgeNormalVectors(2,iEdge) = yEdge(iEdge) - yCell(cell1)
           edgeNormalVectors(3,iEdge) = zEdge(iEdge) - zCell(cell1)
        end if

     else ! this is not a boundary cell
        ! the normal points from the cell 1 to cell2
        ! mrp problem: on periodic domains, vectors on edges of domain point the wrong way.
        if (is_periodic) then
           edgeNormalVectors(1,iEdge) = mpas_fix_periodicity(xCell(cell2), xCell(cell1), x_period) - xCell(cell1)
           edgeNormalVectors(2,iEdge) = mpas_fix_periodicity(yCell(cell2), yCell(cell1), y_period) - yCell(cell1)
           edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)
        else
           edgeNormalVectors(1,iEdge) = xCell(cell2) - xCell(cell1)
           edgeNormalVectors(2,iEdge) = yCell(cell2) - yCell(cell1)
           edgeNormalVectors(3,iEdge) = zCell(cell2) - zCell(cell1)
        end if

     end if
     call mpas_unit_vec_in_r3(edgeNormalVectors(:,iEdge))
  end do

  do iCell=1,nCells
     iEdge = edgesOnCell(1,iCell)
     ! xHat and yHat are a local basis in the plane of the horizontal cell
     ! we arbitrarily choose xHat to point toward the first edge
     rHat = localVerticalUnitVectors(:,iCell)
     normalDotRHat = sum(edgeNormalVectors(:,iEdge)*rHat)
     xHatPlane = edgeNormalVectors(:,iEdge) - normalDotRHat*rHat
     call mpas_unit_vec_in_r3(xHatPlane)

     call mpas_cross_product_in_r3(rHat, xHatPlane, yHatPlane)
     call mpas_unit_vec_in_r3(yHatPlane) ! just to be sure...
     cellTangentPlane(:,1,iCell) = xHatPlane
     cellTangentPlane(:,2,iCell) = yHatPlane
  end do

end subroutine mpas_initialize_vectors!}}}

subroutine mpas_unit_vec_in_r3(xin)!{{{
  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)
  
  real (kind=RKIND), dimension(3), intent(inout) :: xin !< Input/Output: Vector and unit vector
  real (kind=RKIND) :: mag
  mag = sqrt(xin(1)**2+xin(2)**2+xin(3)**2)
  xin(:) = xin(:) / mag
end subroutine mpas_unit_vec_in_r3!}}}

subroutine mpas_cross_product_in_r3(p_1,p_2,p_out)!{{{
  integer, parameter :: RKIND  = selected_real_kind(12)
  real (kind=RKIND), intent(in) :: p_1 (3) !< Input: Vector 1
  real (kind=RKIND), intent(in) :: p_2 (3) !< Input: Vector 2
  real (kind=RKIND), intent(out) :: p_out (3) !< Output: Cross product of vector 1 and vector 2

  p_out(1) = p_1(2)*p_2(3)-p_1(3)*p_2(2)
  p_out(2) = p_1(3)*p_2(1)-p_1(1)*p_2(3)
  p_out(3) = p_1(1)*p_2(2)-p_1(2)*p_2(1)
end subroutine mpas_cross_product_in_r3!}}}

function mpas_fix_periodicity(pxi, xci, xiRef) !{{{
  implicit none
  
  integer, parameter :: RKIND  = selected_real_kind(12)
  real (kind=RKIND), intent(in) :: pxi, xci, xiRef
  real (kind=RKIND) :: mpas_fix_periodicity

  ! local variables
  real (kind=RKIND) :: dist

  dist = pxi - xci

  if (abs(dist) > xiRef * 0.5_RKIND) then
     mpas_fix_periodicity = pxi - (dist/abs(dist)) * xiRef
  else
     mpas_fix_periodicity = pxi
  end if

end function mpas_fix_periodicity !}}}


!subroutine mpas_init_reconstruct(meshPool) ! ,includeHalos)
subroutine mpas_init_reconstruct(nCells, nEdges, maxEdges, R3, &
     edgesOnCell, nEdgesOnCell, &
     xCell, yCell, zCell, xEdge, yEdge, zEdge, edgeNormalVectors, &
     cellTangentPlane, coeffs_reconstruct, is_periodic, x_period, y_period)
  
  implicit none

  !logical, optional, intent(in) :: includeHalos

  integer, parameter :: RKIND  = selected_real_kind(12)
  
  ! temporary arrays needed in the (to be constructed) init procedure
  integer :: nCells, nEdges, maxEdges, R3
  integer, dimension(maxEdges, nCells) :: edgesOnCell
  integer, dimension(nCells) :: nEdgesOnCell
  real(kind=RKIND), dimension(nCells) :: xCell, yCell, zCell
  real(kind=RKIND), dimension(nEdges) :: xEdge, yEdge, zEdge
  real(kind=RKIND), dimension(R3,nEdges) :: edgeNormalVectors
  real(kind=RKIND), dimension(R3,2,nCells) :: cellTangentPlane

  real (kind=RKIND), dimension(R3,maxEdges,nCells) :: coeffs_reconstruct
  logical :: is_periodic
  real(kind=RKIND) :: x_period, y_period

  !f2py intent(in) nCells, nEdges, maxEdges, R3, edgesOnCell, nEdgesOnCell
  !f2py intent(in) xCell, yCell, zCell, xEdge, yEdge, zEdge, edgeNormalVectors
  !f2py intent(in) cellTangentPlane, is_periodic, x_period, y_period
  !f2py intent(out) coeffs_reconstruct

  ! local variables
  integer :: i, iCell, iEdge, pointCount, maxEdgeCount
  real (kind=RKIND) :: mpas_fix_periodicity
  real (kind=RKIND) :: r, cellCenter(3), alpha, tangentPlane(2,3)
  real (kind=RKIND), allocatable, dimension(:,:) :: edgeOnCellLocations, &
       edgeOnCellNormals, coeffs, edgeOnCellLocationsWork, &
       edgeOnCellNormalsWork, coeffsWork
  logical :: includeHalosLocal=.False.


  !if ( present(includeHalos) ) then
  !   includeHalosLocal = includeHalos
  !else
  !   includeHalosLocal = .false.
  !end if

  !if ( includeHalosLocal ) then
  !   call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
  !else
  !   call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCells)
  !end if

  ! init arrays
  coeffs_reconstruct = 0.0

  maxEdgeCount = maxval(nEdgesOnCell)

  allocate(edgeOnCellLocations(maxEdgeCount,3))
  allocate(edgeOnCellNormals(maxEdgeCount,3))
  allocate(coeffs(maxEdgeCount,3))

  ! loop over all cells to be solved on this block
  do iCell=1,nCells
     pointCount = nEdgesOnCell(iCell)
     cellCenter(1) = xCell(iCell)
     cellCenter(2) = yCell(iCell)
     cellCenter(3) = zCell(iCell)

     do i=1,pointCount
        iEdge = edgesOnCell(i,iCell)
        if (is_periodic) then
           edgeOnCellLocations(i,1)  = mpas_fix_periodicity(xEdge(iEdge), cellCenter(1), x_period)
           edgeOnCellLocations(i,2)  = mpas_fix_periodicity(yEdge(iEdge), cellCenter(2), y_period)
           edgeOnCellLocations(i,3)  = zEdge(iEdge)
        else
           edgeOnCellLocations(i,1)  = xEdge(iEdge)
           edgeOnCellLocations(i,2)  = yEdge(iEdge)
           edgeOnCellLocations(i,3)  = zEdge(iEdge)
        end if
        edgeOnCellNormals(i,:)  = edgeNormalVectors(:, iEdge)
     end do

     alpha = 0.0
     do i=1,pointCount
        r = sqrt(sum((cellCenter - edgeOnCellLocations(i,:))**2))
        alpha = alpha + r
     enddo
     alpha = alpha/pointCount

     tangentPlane(1,:) = cellTangentPlane(:,1,iCell)
     tangentPlane(2,:) = cellTangentPlane(:,2,iCell)

     allocate(edgeOnCellLocationsWork(pointCount,3))
     allocate(edgeOnCellNormalsWork(pointCount,3))
     allocate(coeffsWork(pointCount,3))

     edgeOnCellLocationsWork = edgeOnCellLocations(1:pointCount,:)
     edgeOnCellNormalsWork = edgeOnCellNormals(1:pointCount,:)

     call mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pointCount, &
          edgeOnCellLocationsWork, edgeOnCellNormalsWork, &
          cellCenter, alpha, tangentPlane, coeffsWork)

     coeffs(1:pointCount,:) = coeffsWork

     deallocate(edgeOnCellLocationsWork)
     deallocate(edgeOnCellNormalsWork)
     deallocate(coeffsWork)

     do i=1,pointCount
        coeffs_reconstruct(:,i,iCell) = coeffs(i,:)
     end do
  enddo   ! iCell

  deallocate(edgeOnCellLocations)
  deallocate(edgeOnCellNormals)
  deallocate(coeffs)

end subroutine mpas_init_reconstruct!}}}

!subroutine mpas_reconstruct_2d(meshPool, u, uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, uReconstructMeridional, includeHalos)!{{{

subroutine mpas_reconstruct_2d(nCells, nEdges, maxEdges, nVertLevels, R3, &
     on_a_sphere, &
     edgesOnCell, nEdgesOnCell, latCell, lonCell, coeffs_reconstruct, u, &
     uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal, &
     uReconstructMeridional)
  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)

  integer :: nCells, nEdges, nVertLevels, maxEdges, R3
  real (kind=RKIND), dimension(nVertLevels,nEdges) :: u
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructX
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructY
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZ
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructZonal
  real (kind=RKIND), dimension(nVertLevels,nCells) :: uReconstructMeridional
  !logical, optional, intent(in) :: includeHalos

  integer, dimension(maxEdges,nCells) :: edgesOnCell
  integer, dimension(nCells) :: nEdgesOnCell
  real(kind=RKIND), dimension(nCells) :: latCell, lonCell
  real (kind=RKIND), dimension(R3,maxEdges,nCells) :: coeffs_reconstruct
  logical :: on_a_sphere

  !f2py intent(in) nCells, nEdges, nVertLevels, maxEdges, R3, on_a_sphere
  !f2py intent(in) edgesOnCell, nEdgesOnCell, latCell, lonCell,
  !f2py intent(in) coeffs_reconstruct, u
  !f2py intent(out) uReconstructX, uReconstructY, uReconstructZ, uReconstructZonal
  !f2py intent(out) uReconstructMeridional

  ! local variable
  logical :: includeHalosLocal
  integer :: iCell,iEdge, i
  real (kind=RKIND) :: clat, slat, clon, slon

  !if ( present(includeHalos) ) then
  !   includeHalosLocal = includeHalos
  !else
  !   includeHalosLocal = .false.
  !end if

  ! loop over cell centers
  !$omp do schedule(runtime)

  do iCell = 1, nCells
     ! initialize the reconstructed vectors
     uReconstructX(:,iCell) = 0.0
     uReconstructY(:,iCell) = 0.0
     uReconstructZ(:,iCell) = 0.0

     ! a more efficient reconstruction where rbf_values*matrix_reconstruct
     ! has been precomputed in coeffs_reconstruct
     do i=1,nEdgesOnCell(iCell)
        iEdge = edgesOnCell(i,iCell)
        uReconstructX(:,iCell) = uReconstructX(:,iCell) &
             + coeffs_reconstruct(1,i,iCell) * u(:,iEdge)
        uReconstructY(:,iCell) = uReconstructY(:,iCell) &
             + coeffs_reconstruct(2,i,iCell) * u(:,iEdge)
        uReconstructZ(:,iCell) = uReconstructZ(:,iCell) &
             + coeffs_reconstruct(3,i,iCell) * u(:,iEdge)
     enddo
  enddo   ! iCell
  !$omp end do

  if (on_a_sphere) then
     !$omp do schedule(runtime)
     do iCell = 1, nCells
        clat = cos(latCell(iCell))
        slat = sin(latCell(iCell))
        clon = cos(lonCell(iCell))
        slon = sin(lonCell(iCell))
        uReconstructZonal(:,iCell) = -uReconstructX(:,iCell)*slon + &
             uReconstructY(:,iCell)*clon
        uReconstructMeridional(:,iCell) = -(uReconstructX(:,iCell)*clon       &
             + uReconstructY(:,iCell)*slon)*slat &
             + uReconstructZ(:,iCell)*clat
     end do
     !$omp end do
  else
     !$omp do schedule(runtime)
     do iCell = 1, nCells
        uReconstructZonal     (:,iCell) = uReconstructX(:,iCell)
        uReconstructMeridional(:,iCell) = uReconstructY(:,iCell)
     end do
     !$omp end do
  end if

end subroutine mpas_reconstruct_2d!}}}


subroutine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs(pointCount, &!{{{
     sourcePoints, unitVectors, destinationPoint, &
     alpha, planeBasisVectors, coefficients)

  integer, parameter :: RKIND  = selected_real_kind(12)

  integer, intent(in) :: pointCount !< Input: Number of points
  real(kind=RKIND), dimension(pointCount,3), intent(in) :: sourcePoints !< Input: List of points
  real(kind=RKIND), dimension(pointCount,3), intent(in) :: unitVectors !< Input: List of unit vectors
  real(kind=RKIND), dimension(3), intent(in) :: destinationPoint !< Input: Destination point
  real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
  real(kind=RKIND), dimension(2,3) :: planeBasisVectors !< Input: Basis vectors for interpolation plane
  real(kind=RKIND), dimension(pointCount, 3), intent(out) :: coefficients !< Output: List of coefficients

  integer :: i
  integer :: matrixSize

  real(kind=RKIND), dimension(pointCount,2) :: planarSourcePoints
  real(kind=RKIND), dimension(pointCount,2) :: planarUnitVectors
  real(kind=RKIND), dimension(2) :: planarDestinationPoint

  real(kind=RKIND), dimension(pointCount+2, pointCount+2) :: matrix, matrixCopy
  real(kind=RKIND), dimension(pointCount, pointCount) :: matrixWork
  real(kind=RKIND), dimension(pointCount+2, 2) :: rhs, coeffs
  real(kind=RKIND), dimension(pointCount,2) :: rhsWork
  integer, dimension(pointCount+2) :: pivotIndices

  matrixSize = pointCount+2 ! space for constant vector in plane

  matrix = 0.0
  rhs = 0.0
  coeffs = 0.0

  do i = 1, pointCount
     planarSourcePoints(i,1) = sum(sourcePoints(i,:)*planeBasisVectors(1,:))
     planarSourcePoints(i,2) = sum(sourcePoints(i,:)*planeBasisVectors(2,:))
     planarUnitVectors(i,1) = sum(unitVectors(i,:)*planeBasisVectors(1,:))
     planarUnitVectors(i,2) = sum(unitVectors(i,:)*planeBasisVectors(2,:))
  end do
  planarDestinationPoint(1) = sum(destinationPoint*planeBasisVectors(1,:))
  planarDestinationPoint(2) = sum(destinationPoint*planeBasisVectors(2,:))

  call mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, 2, &
       planarSourcePoints, planarUnitVectors, planarDestinationPoint, &
       alpha, matrixWork, rhsWork)

  matrix(1:pointCount,1:pointCount) = matrixWork
  rhs(1:pointCount,:) = rhsWork

  do i = 1, pointCount
     matrix(i,pointCount+1:pointCount+2) = planarUnitVectors(i,:)
     matrix(pointCount+1:pointCount+2,i) = matrix(i,pointCount+1:pointCount+2)
  end do
  do i = 1,2
     rhs(pointCount+i,i) = 1.0 ! the unit vector in the ith direction
  end do

  ! solve each linear system
  matrixCopy = matrix
  call mpas_legs(matrix, matrixSize, rhs(:,1), coeffs(:,1), pivotIndices)
  call mpas_legs(matrixCopy, matrixSize, rhs(:,2), coeffs(:,2), pivotIndices)


  do i = 1,3
     coefficients(:,i) = planeBasisVectors(1,i)*coeffs(1:pointCount,1) &
          + planeBasisVectors(2,i)*coeffs(1:pointCount,2)
  end do

end subroutine mpas_rbf_interp_func_3D_plane_vec_const_dir_comp_coeffs !}}}

subroutine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs(pointCount, dimensions, &!{{{
     sourcePoints, unitVectors, destinationPoint, &
     alpha, matrix, rhs)

  integer, parameter :: RKIND  = selected_real_kind(12)
  integer, intent(in) :: pointCount !< Input: Number of points
  integer, intent(in) :: dimensions !< Input: Number of dimensions
  real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: sourcePoints !< Input: List of points
  real(kind=RKIND), dimension(pointCount,dimensions), intent(in) :: unitVectors !< Input: List of unit vectors
  real(kind=RKIND), dimension(dimensions), intent(in) :: destinationPoint !< Input: Destination point
  real(kind=RKIND), intent(in) :: alpha !< Input: Characteristic length scale of RBFs
  real(kind=RKIND), dimension(pointCount,pointCount), intent(out) :: matrix !< Output: Matrix
  real(kind=RKIND), dimension(pointCount,dimensions), intent(out) :: rhs !< Output: Right hand side

  integer :: i, j

  real(kind=RKIND) :: evaluate_rbf
  real(kind=RKIND) :: rSquared, rbfValue, unitVectorDotProduct

  do j = 1, pointCount
     do i = j, pointCount
        rSquared = sum((sourcePoints(i,:)-sourcePoints(j,:))**2)/alpha**2
        rbfValue = evaluate_rbf(rSquared)
        unitVectorDotProduct = sum(unitVectors(i,:)*unitVectors(j,:))
        matrix(i,j) = rbfValue*unitVectorDotProduct
        matrix(j,i) = matrix(i,j)
     end do
  end do

  do j = 1, pointCount
     rSquared = sum((destinationPoint-sourcePoints(j,:))**2)/alpha**2
     rhs(j,:) = evaluate_rbf(rSquared)*unitVectors(j,:)
  end do

end subroutine mpas_set_up_vector_dirichlet_rbf_matrix_and_rhs!}}}

function evaluate_rbf(rSquared) result(rbfValue)!{{{
  integer, parameter :: RKIND  = selected_real_kind(12)
  real(kind=RKIND), intent(in) :: rSquared !< Input: Squared value of r
  real(kind=RKIND) :: rbfValue

  ! inverse multiquadratic
  rbfValue = 1/sqrt(1 + rSquared)

    end function evaluate_rbf!}}}

subroutine mpas_legs (A,N,B,X,INDX)!{{{
  IMPLICIT NONE

  integer, parameter :: RKIND  = selected_real_kind(12)
  integer, INTENT (IN) :: N !< Input: Size of matrix and vectors
  integer, INTENT (OUT), DIMENSION (N) :: INDX !< Output: Pivot vector
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A !< Input/Output: Matrix
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N) :: B !< Input/Output: Right hand side vector
  real(kind=RKIND), INTENT (OUT), DIMENSION (N) :: X !< Output: Solution

  integer :: I,J
  !
  CALL elgs (A,N,INDX)
  !
  DO I = 1, N-1
     DO J = I+1, N
        B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
     END DO
  END DO
  !
  X(N) = B(INDX(N))/A(INDX(N),N)
  DO I = N-1, 1, -1
     X(I) = B(INDX(I))
     DO J = I+1, N
        X(I) = X(I)-A(INDX(I),J)*X(J)
     END DO
     X(I) =  X(I)/A(INDX(I),I)
  END DO
  !
END subroutine mpas_legs!}}}

subroutine elgs (A,N,INDX)!{{{
  !
  ! subroutine to perform the partial-pivoting Gaussian elimination.
  ! A(N,N) is the original matrix in the input and transformed matrix
  ! plus the pivoting element ratios below the diagonal in the output.
  ! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
  !
  IMPLICIT NONE

  integer, parameter :: RKIND  = selected_real_kind(12)
  integer, INTENT (IN) :: N !< Input: Size of matrix
  integer, INTENT (OUT), DIMENSION (N) :: INDX !< Output: Pivot vector
  real(kind=RKIND), INTENT (INOUT), DIMENSION (N,N) :: A !< Input/Output: Matrix and solution
  integer :: I,J,K,ITMP
  real(kind=RKIND) :: C1,PI,PI1,PJ
  real(kind=RKIND), DIMENSION (N) :: C
  !
  ! Initialize the index
  !
  DO I = 1, N
     INDX(I) = I
  END DO
  !
  ! Find the rescaling factors, one from each row
  !
  DO I = 1, N
     C1= 0.0
     DO J = 1, N
        !C1 = AMAX1(C1,ABS(A(I,J)))
        C1 = MAX(C1,ABS(A(I,J)))
     END DO
     C(I) = C1
  END DO
  !
  ! Search the pivoting (largest) element from each column
  !
  DO J = 1, N-1
     PI1 = 0.0
     DO I = J, N
        PI = ABS(A(INDX(I),J))/C(INDX(I))
        IF (PI>PI1) THEN
           PI1 = PI
           K   = I
        ENDIF
     END DO
     !
     ! Interchange the rows via INDX(N) to record pivoting order
     !
     ITMP    = INDX(J)
     INDX(J) = INDX(K)
     INDX(K) = ITMP
     DO I = J+1, N
        PJ  = A(INDX(I),J)/A(INDX(J),J)
        !
        ! Record pivoting ratios below the diagonal
        A(INDX(I),J) = PJ
        !
        ! Modify other elements accordingly
        !
        DO K = J+1, N
           A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
        END DO
     END DO
  END DO
  !
END subroutine elgs!}}}
