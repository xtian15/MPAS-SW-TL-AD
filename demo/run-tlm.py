# modify the source directory accordingly
import sys
from datetime import timedelta

import netCDF4 as nc
import numpy as np

sys.path.append('../src/')

from module_sw_mpas import mpas_sw_module as mpsw
from mpas_namelist import namelist
from mpas_sw_driver import read_configs, read_dims, read_vars, \
    initial_conditions, clock_namelist


def run_sw_tlm():
    nml = namelist(nmlname='namelist.sw.x1.10242')
    read_configs(nml)
    read_dims(nml)
    read_vars(nml)

    initial_conditions(nml)
    mpsw.var_allocation_tlm()
    mpsw.sw_mpas_init_block()

    # initializing TLM variables
    alpha = 1.E-3
    mpsw.u_tl[0, 0] = mpsw.u[0, 0] * alpha
    mpsw.h_tl[0, 0] = mpsw.h[0, 0] * alpha
    # end initializing TLM variables

    clock_start, clock_end = clock_namelist(nml)
    today = clock_start
    unorm, hlist, ulist, vlist = [], [], [], []
    while today <= clock_end:
        if (today - clock_start).total_seconds() % nml.output_interval == 0:
            print(today, mpsw.u[0].min(), mpsw.u[0].max())
            ulist.append(mpsw.ureconstructzonal_tl[0].copy())
            vlist.append(mpsw.ureconstructmeridional_tl[0].copy())
            hlist.append(mpsw.h_tl[0, 0].copy())
            unorm.append(mpsw.u_tl[0, 0].copy())
        mpsw.sw_rk4_tlm()
        today += timedelta(seconds=int(mpsw.config_dt))

    # outputing results
    ulist, vlist, hlist = np.array(ulist), np.array(vlist), np.array(hlist)
    unorm = np.array(unorm)

    r2d = 180. / np.pi
    outname = nml.file_output.replace('.nc', '.tlm.nc')
    with nc.Dataset(outname, 'w') as of:
        of.createDimension('nCell', mpsw.latcell.shape[0])
        of.createDimension('nTime', ulist.shape[0])
        of.createVariable('latCell', 'f4', ('nCell'))[:] = mpsw.latcell * r2d
        of.createVariable('lonCell', 'f4', ('nCell'))[:] = mpsw.loncell * r2d
        of.createVariable('ux_tl', 'f8', ('nTime', 'nCell'))[:] = ulist
        of.createVariable('uy_tl', 'f8', ('nTime', 'nCell'))[:] = vlist
        of.createVariable('h_tl', 'f8', ('nTime', 'nCell'))[:] = hlist

        of.createDimension('nEdge', mpsw.latedge.shape[0])
        of.createVariable('latEdge', 'f4', ('nEdge'))[:] = mpsw.latedge * r2d
        of.createVariable('lonEdge', 'f4', ('nEdge'))[:] = mpsw.lonedge * r2d
        of.createVariable('u_tl', 'f8', ('nTime', 'nEdge'))[:] = unorm


if __name__ == '__main__':
    run_sw_tlm()
