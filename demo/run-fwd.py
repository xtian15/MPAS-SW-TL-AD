import numpy as np
import netCDF4 as nc
from datetime import datetime, timedelta

# modify the source directory accordingly
import sys
sys.path.append('../src/')

from module_sw_mpas import mpas_sw_module as mpsw
from mpas_namelist import namelist
from mpas_sw_driver import read_configs, read_dims, read_vars, \
    initial_conditions, clock_namelist, insert_huv

def run_sw_fwd():
    nml=namelist(nmlname='namelist.sw.x1.10242')
    read_configs(nml)
    read_dims(nml)
    read_vars(nml)

    initial_conditions(nml)
    """
    # low res (insert results from sources other than ERA5)
    uedge, hcell=insert_huv(lat_hr, lon_hr, h_hr, ux_hr, uy_hr)
    mpsw.u[0,0]=uedge.copy()
    mpsw.h[0,0]=hcell.copy()
    """
        
    mpsw.sw_mpas_init_block()
    clock_start, clock_end=clock_namelist(nml)

    today=clock_start
    ulist, vlist, hlist=[], [], []
    unorm=[]
    while today<=clock_end:
        if (today-clock_start).total_seconds()%nml.output_interval==0:
            print(today, mpsw.u[0].min(), mpsw.u[0].max())
            ulist.append(mpsw.ureconstructzonal[0].copy())
            vlist.append(mpsw.ureconstructmeridional[0].copy())
            hlist.append(mpsw.h[0,0].copy())
            unorm.append(mpsw.u[0,0].copy())
        
        mpsw.sw_rk4()

        today+=timedelta(seconds=int(mpsw.config_dt))
    ulist, vlist, hlist=np.array(ulist), np.array(vlist), np.array(hlist)
    unorm=np.array(unorm)

    r2d=180./np.pi
    with nc.Dataset(nml.file_output, 'w') as of:
        of.createDimension('nCell', mpsw.latcell.shape[0])
        of.createDimension('nTime', ulist.shape[0])
        of.createVariable('latCell', 'f4', ('nCell'))[:]=mpsw.latcell*r2d
        of.createVariable('lonCell', 'f4', ('nCell'))[:]=mpsw.loncell*r2d
        of.createVariable('ux', 'f8', ('nTime', 'nCell'))[:]=ulist
        of.createVariable('uy', 'f8', ('nTime', 'nCell'))[:]=vlist
        of.createVariable('h',  'f8', ('nTime', 'nCell'))[:]=hlist

        of.createDimension('nEdge', mpsw.latedge.shape[0])
        of.createVariable('latEdge', 'f4', ('nEdge'))[:]=mpsw.latedge*r2d
        of.createVariable('lonEdge', 'f4', ('nEdge'))[:]=mpsw.lonedge*r2d
        of.createVariable('u', 'f8', ('nTime', 'nEdge'))[:]=unorm
        
if __name__=='__main__':
    run_sw_fwd()
