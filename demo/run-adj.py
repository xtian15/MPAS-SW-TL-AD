import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta
from scipy import interpolate

import sys
sys.path.append('../src/')

from module_sw_mpas import mpas_sw_module as mpsw
from mpas_namelist import namelist
from mpas_sw_driver import read_configs, read_dims, read_vars, \
    initial_conditions, clock_namelist

def run_sw_adj():

    # ----- nonlinear state trajectory -----
    fwdname='x1.10242.state.nc'
    fwdfile=nc.Dataset(fwdname, 'r')
    # ----- end nonlinear state trajectory -----
    
    nml=namelist(nmlname='namelist.sw.x1.10242')
    
    read_configs(nml)
    read_dims(nml)
    read_vars(nml)

    initial_conditions(nml)
    mpsw.var_allocation_adj()
    mpsw.sw_mpas_init_block()
    clock_start, clock_end=clock_namelist(nml)

    # ----- adj initial conditions -----
    mpsw.u_ad[0,0,0]=1.
    mpsw.h_ad[0,0,0]=1.
    # ----- end adj initial conditions -----

    today=clock_end
    ulist, vlist, hlist=[], [], []
    unorm=[]
    itimestep=int( (clock_end-clock_start).total_seconds()/mpsw.config_dt )
    while today>=clock_start:
        print('{:%Y-%m-%d %H:%M}'.format(today), mpsw.u_ad.sum(), mpsw.h_ad.sum())
        
        mpsw.u[0,0]=fwdfile['u'][itimestep]
        mpsw.h[0,0]=fwdfile['h'][itimestep]
        if (today-clock_start).total_seconds()%nml.output_interval==0:
            ulist.append(mpsw.ureconstructzonal_ad[0].copy())
            vlist.append(mpsw.ureconstructmeridional_ad[0].copy())
            hlist.append(mpsw.h_ad[0,0].copy())
            unorm.append(mpsw.u_ad[0,0].copy())
        
        mpsw.sw_rk4_adj()

        itimestep-=1
        today-=timedelta(seconds=int(mpsw.config_dt))

    ulist, vlist=np.array(ulist), np.array(vlist)
    hlist=np.array(hlist)
    unorm=np.array(unorm)
        
    r2d=180./np.pi
    outname=nml.file_output.replace('.nc', '.adj.nc')
    with nc.Dataset(outname, 'w') as of:
        of.createDimension('nTime', hlist.shape[0])
        of.createDimension('nCell', mpsw.latcell.shape[0])
        of.createDimension('nEdge', mpsw.latedge.shape[0])

        of.createVariable('latCell', 'f4', ('nCell'))[:]=mpsw.latcell*r2d
        of.createVariable('lonCell', 'f4', ('nCell'))[:]=mpsw.loncell*r2d
        of.createVariable('latEdge', 'f4', ('nEdge'))[:]=mpsw.latedge*r2d
        of.createVariable('lonEdge', 'f4', ('nEdge'))[:]=mpsw.lonedge*r2d

        of.createVariable('ux_ad', 'f4', ('nTime', 'nCell'))[:]=ulist
        of.createVariable('uy_ad', 'f4', ('nTime', 'nCell'))[:]=vlist
        of.createVariable('h_ad',  'f4', ('nTime', 'nCell'))[:]=hlist
        of.createVariable('u_ad',  'f4', ('nTime', 'nEdge'))[:]=unorm
        
if __name__=='__main__':
    run_sw_adj()
