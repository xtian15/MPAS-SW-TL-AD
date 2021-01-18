subroutine sphere_distance(lat1, lon1, nlat, nlon, gridlat, gridlon, dist)

 implicit none

 integer, parameter :: RKIND  = selected_real_kind(12)
 integer :: nlat, nlon
 real (kind=RKIND), intent(in) :: lat1, lon1
 real (kind=RKIND), dimension(nlat,nlon) :: gridlat, gridlon
 real (kind=RKIND), dimension(nlat,nlon) :: dist

 ! local
 real (kind=RKIND), dimension(nlat,nlon) :: arg1
 real (kind=RKIND), parameter :: radius=1.0_RKIND

 arg1 = sqrt( sin(0.5*(gridlat-lat1))**2 +  &
              cos(lat1)*cos(gridlat)*sin(0.5*(gridlon-lon1))**2 )
 dist = 2.*radius*asin(arg1)

end subroutine sphere_distance

subroutine grid2mpas(nlat, nlon, gridlat, gridlon, gridvar, &
     ncell, latcell, loncell, var)

  implicit none

  integer, parameter :: RKIND  = selected_real_kind(12)
  integer :: nlat, nlon, ncell
  real(kind=RKIND), dimension(nlat, nlon) :: gridlat   ! [90, -90]
  real(kind=RKIND), dimension(nlat, nlon) :: gridlon   ! [0, 360]
  real(kind=RKIND), dimension(nlat, nlon) :: gridvar
  real(kind=RKIND), dimension(ncell) :: latcell, loncell, var

  !f2py intent(in) nlat, nlon, gridlat, gridlon, gridvar
  !f2py intent(in) ncell, latcell, loncell
  !f2py intent(out) var

  ! ----- local vars -----
  integer :: icell
  integer :: latloc1, latloc2, lonloc1, lonloc2
  integer, dimension(1) :: loc1d
  integer, dimension(2) :: shp
  real(kind=RKIND), dimension(nlat) :: lat1d
  real(kind=RKIND), dimension(nlon) :: lon1d

  where(loncell<0) latcell=loncell+360
  where(gridlon<0) gridlat=gridlat+360

  lat1d=gridlat(:,1)
  lon1d=gridlon(1,:)
  do icell=1, ncell

     loc1d=minloc(abs(latcell(icell)-lat1d))
     latloc1=loc1d(1)-1
     latloc2=loc1d(1)+1
     if(latloc1<1   ) latloc1=1
     if(latloc2>nlat) latloc2=nlat

     loc1d=minloc(abs(loncell(icell)-lon1d))
     lonloc1=loc1d(1)-1
     lonloc2=loc1d(1)+1
     if(lonloc1<1   ) lonloc1=1
     if(lonloc2>nlon) lonloc2=nlon

     shp=shape(gridvar(latloc1:latloc2,lonloc1:lonloc2))     
     var(icell)=SUM(gridvar(latloc1:latloc2,lonloc1:lonloc2))/(shp(1)*shp(2))

  end do

  return
end subroutine grid2mpas
     
        
