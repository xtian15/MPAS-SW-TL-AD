import numpy as np
import netCDF4 as nc

from scipy import interpolate

from pandas import to_timedelta
from datetime import datetime, timedelta
from module_sw_mpas import mpas_sw_module as mpsw
from mpas_namelist import namelist


def read_configs(nml):
    nlatm = nml.nlatm
    mpsw.config_test_case = nlatm['sw_model']['config_test_case']
    mpsw.config_dt = nlatm['sw_model']['config_dt']
    mpsw.config_calendar_type = nlatm['sw_model']['config_calendar_type']
    mpsw.config_stats_interval = nlatm['sw_model']['config_stats_interval']
    mpsw.config_h_ScaleWithMesh = nlatm['sw_model']['config_h_ScaleWithMesh']
    mpsw.config_h_mom_eddy_visc2 = nlatm['sw_model']['config_h_mom_eddy_visc2']
    mpsw.config_h_mom_eddy_visc4 = nlatm['sw_model']['config_h_mom_eddy_visc4']
    mpsw.config_h_tracer_eddy_diff2 = nlatm['sw_model']['config_h_tracer_eddy_diff2']
    mpsw.config_h_tracer_eddy_diff4 = nlatm['sw_model']['config_h_tracer_eddy_diff4']
    mpsw.config_thickness_adv_order = nlatm['sw_model']['config_thickness_adv_order']
    mpsw.config_tracer_adv_order = nlatm['sw_model']['config_tracer_adv_order']
    mpsw.config_positive_definite = nlatm['sw_model']['config_positive_definite']
    mpsw.config_monotonic = nlatm['sw_model']['config_monotonic']
    mpsw.config_wind_stress = nlatm['sw_model']['config_wind_stress']
    mpsw.config_bottom_drag = nlatm['sw_model']['config_bottom_drag']
    mpsw.config_apvm_upwinding = nlatm['sw_model']['config_apvm_upwinding']
    mpsw.config_num_halos = nlatm['sw_model']['config_num_halos']


def read_dims(nml):
    with nc.Dataset(nml.file_t0, 'r') as infile:
        try:
            mpsw.nvertlevels = infile.dimensions['nVertLevels'].size
        except:
            mpsw.nvertlevels = 1
        mpsw.ncells = infile.dimensions['nCells'].size
        mpsw.nedges = infile.dimensions['nEdges'].size
        mpsw.nvertices = infile.dimensions['nVertices'].size
        mpsw.maxedges = infile.dimensions['maxEdges'].size
        mpsw.maxedges2 = infile.dimensions['maxEdges2'].size
        mpsw.vertexdegree = infile.dimensions['vertexDegree'].size
        mpsw.r3 = 3  # ? under question
    mpsw.ntracers = 1
    mpsw.ncellssolve = mpsw.ncells
    mpsw.nedgessolve = mpsw.nedges
    mpsw.nverticessolve = mpsw.nvertices
    mpsw.cellstart = 1
    mpsw.cellend = mpsw.ncells
    mpsw.vertexstart = 1
    mpsw.vertexend = mpsw.nvertices
    mpsw.edgestart = 1
    mpsw.edgeend = mpsw.nedges
    mpsw.cellsolvestart = 1
    mpsw.cellsolveend = mpsw.ncells
    mpsw.vertexsolvestart = 1
    mpsw.vertexsolveend = mpsw.nvertices
    mpsw.edgesolvestart = 1
    mpsw.edgesolveend = mpsw.nedges


def read_vars(nml):
    mpsw.var_allocation()  # shapes in place
    with nc.Dataset(nml.file_t0, 'r') as infile:
        """
        if infile.on_a_sphere.strip()=='YES':
            mpsw.on_a_sphere=True
        else:
            mpsw.on_a_sphere=False
        """
        mpsw.on_a_sphere = True
        mpsw.is_periodic = False

        try:
            mpsw.sphere_radius = infile.sphere_radius
        except:
            mpsw.sphere_radius = 1.

        # ===== state =====
        # mpsw.u[0,:,:]=np.squeeze(infile['u'][()].T)
        # mpsw.h[0,:,:]=np.squeeze(infile['h'][()].T)
        # mpsw.tracers[0]=0. # no tracers enabled

        # ===== mesh =====
        mpsw.latcell = infile['latCell'][()].T
        mpsw.loncell = infile['lonCell'][()].T
        mpsw.xcell = infile['xCell'][()].T
        mpsw.ycell = infile['yCell'][()].T
        mpsw.zcell = infile['zCell'][()].T
        mpsw.latedge = infile['latEdge'][()].T
        mpsw.lonedge = infile['lonEdge'][()].T
        mpsw.xedge = infile['xEdge'][()].T
        mpsw.yedge = infile['yEdge'][()].T
        mpsw.zedge = infile['zEdge'][()].T
        mpsw.latvertex = infile['latVertex'][()].T
        mpsw.lonvertex = infile['lonVertex'][()].T
        mpsw.xvertex = infile['xVertex'][()].T
        mpsw.yvertex = infile['yVertex'][()].T
        mpsw.zvertex = infile['zVertex'][()].T
        mpsw.cellsonedge = infile['cellsOnEdge'][()].T
        mpsw.nedgesoncell = infile['nEdgesOnCell'][()].T
        mpsw.nedgesonedge = infile['nEdgesOnEdge'][()].T
        mpsw.edgesoncell = infile['edgesOnCell'][()].T
        mpsw.edgesonedge = infile['edgesOnEdge'][()].T
        mpsw.weightsonedge = infile['weightsOnEdge'][()].T
        mpsw.dvedge = infile['dvEdge'][()].T
        mpsw.dcedge = infile['dcEdge'][()].T
        mpsw.angleedge = infile['angleEdge'][()].T
        mpsw.areacell = infile['areaCell'][()].T
        mpsw.areatriangle = infile['areaTriangle'][()].T
        mpsw.cellsoncell = infile['cellsOnCell'][()].T
        mpsw.verticesoncell = infile['verticesOnCell'][()].T
        mpsw.verticesonedge = infile['verticesOnEdge'][()].T
        mpsw.edgesonvertex = infile['edgesOnVertex'][()].T
        mpsw.cellsonvertex = infile['cellsOnVertex'][()].T
        mpsw.kiteareasonvertex = infile['kiteAreasOnVertex'][()].T
        mpsw.meshdensity = infile['meshDensity'][()].T


def clock_namelist(nml):
    nlatm = nml.nlatm

    clock_start = datetime.strptime(nlatm['sw_model']['config_start_time'],
                                    '%Y-%m-%d_%H:%M:%S')
    clock_duration_str = nlatm['sw_model']['config_run_duration'].replace(
        '_', ' days ')
    clock_duration = to_timedelta(clock_duration_str)
    clock_end = clock_start + clock_duration

    return [clock_start, clock_end]


def era2mpas_old(latGrid, lonGrid, varGrid, lat, lon_in):
    latGrid = latGrid.flatten()
    lonGrid = lonGrid.flatten()
    lonGrid[lonGrid < 0] += 360.

    lon = lon_in.copy()
    lon[lon < 0] += 360

    var = interpolate.griddata(np.array([lonGrid, latGrid]).T, varGrid.flatten(),
                               np.array([lon, lat]).T, method='linear')

    ave = var[np.logical_not(np.isnan(var))].mean()
    var[np.isnan(var)] = ave
    return var


def era2mpas(latGrid, lonGrid, varGrid, lat, lon_in):
    from grid2mpas import grid2mpas

    lon = lon_in.copy()
    lon[lon < 0] += 360

    var = grid2mpas(latGrid, lonGrid, varGrid, lat, lon)

    ave = var[np.logical_not(np.isnan(var))].mean()
    var[np.isnan(var)] = ave
    return var


def mpsw_init_real(nml):
    # Scale all distances and areas from a unit sphere to one with radius a
    mpsw.xcell[()] *= mpsw.a
    mpsw.ycell[()] *= mpsw.a
    mpsw.zcell[()] *= mpsw.a
    mpsw.xvertex[()] *= mpsw.a
    mpsw.yvertex[()] *= mpsw.a
    mpsw.zvertex[()] *= mpsw.a
    mpsw.xedge[()] *= mpsw.a
    mpsw.yedge[()] *= mpsw.a
    mpsw.zedge[()] *= mpsw.a
    mpsw.dvedge[()] *= mpsw.a
    mpsw.dcedge[()] *= mpsw.a
    mpsw.areacell[()] *= mpsw.a ** 2.0
    mpsw.areatriangle[()] *= mpsw.a ** 2.0
    mpsw.kiteareasonvertex[()] *= mpsw.a ** 2.0

    alpha = 0.0
    # ----- initialize the Coriolis parameters -----
    mpsw.fedge[()] = 2. * mpsw.omega * (
            -np.cos(mpsw.lonedge) * np.cos(mpsw.latedge) * np.sin(alpha)
            + np.sin(mpsw.latedge) * np.cos(alpha))
    mpsw.fvertex[()] = 2. * mpsw.omega * (
            -np.cos(mpsw.lonvertex) * np.cos(mpsw.latvertex) * np.sin(alpha)
            + np.sin(mpsw.latvertex) * np.cos(alpha))
    # ----- end initialize the Coriolis parameters -----

    # initializing with ERA5
    # eraname: ERA5 reanalysis dataset in netCDF
    # anatime: datetime of interest within the file eraname

    day1 = datetime(1900, 1, 1)
    clock_start, clock_end = clock_namelist(nml)
    dhour = int((clock_start - day1).total_seconds() / 3600)

    eraname = nml.nlatm['io']['eraname']

    with nc.Dataset(eraname, 'r') as infile:
        time = infile['time'][()]
        itime = np.argmin(np.absolute(dhour - time))
        latGrid = infile['latitude'][()]
        lonGrid = infile['longitude'][()]
        lonGrid, latGrid = np.meshgrid(lonGrid, latGrid)
        zGrid = infile['z'][itime] / mpsw.gravity  # (nlat,nlon)
        uGrid = infile['u'][itime]
        vGrid = infile['v'][itime]

    r2d = 180 / np.pi
    zCell = era2mpas(latGrid, lonGrid, zGrid, mpsw.latcell * r2d, mpsw.loncell * r2d)
    mpsw.h[0, 0] = zCell
    uEdge = era2mpas(latGrid, lonGrid, uGrid, mpsw.latedge * r2d, mpsw.lonedge * r2d)
    vEdge = era2mpas(latGrid, lonGrid, vGrid, mpsw.latedge * r2d, mpsw.lonedge * r2d)
    uERA = np.cos(mpsw.angleedge) * uEdge + np.sin(mpsw.angleedge) * vEdge
    mpsw.u[0, 0] = uERA


def initial_conditions(nml):
    nlatm = nml.nlatm
    case_no = nlatm['sw_model']['config_test_case']
    if case_no == 1:
        mpsw.sw_test_case_1()
    elif case_no == 2:
        mpsw.sw_test_case_2()
    elif case_no == 5:
        mpsw.sw_test_case_5()
    elif case_no == 6:
        mpsw.sw_test_case_6()
    elif case_no == 7:
        mpsw_init_real(nml)
    else:
        exit('other cases not supported yet')


def insert_huv(lat, lon, h, u, v):
    r2d = 180. / np.pi
    latCell, lonCell = mpsw.latcell * r2d, mpsw.loncell * r2d
    latEdge, lonEdge = mpsw.latedge * r2d, mpsw.lonedge * r2d

    hcell = interpolate.griddata(np.array([lon, lat]).T, h,
                                 np.array([lonCell, latCell]).T, method='linear')
    hfill = interpolate.griddata(np.array([lon, lat]).T, h,
                                 np.array([lonCell, latCell]).T, method='nearest')
    hcell[np.isnan(hcell)] = hfill[np.isnan(hcell)]

    uxedge = interpolate.griddata(np.array([lon, lat]).T, u,
                                  np.array([lonEdge, latEdge]).T, method='linear')
    uxfill = interpolate.griddata(np.array([lon, lat]).T, u,
                                  np.array([lonEdge, latEdge]).T, method='nearest')
    uxedge[np.isnan(uxedge)] = uxfill[np.isnan(uxedge)]

    uyedge = interpolate.griddata(np.array([lon, lat]).T, v,
                                  np.array([lonEdge, latEdge]).T, method='linear')
    uyfill = interpolate.griddata(np.array([lon, lat]).T, v,
                                  np.array([lonEdge, latEdge]).T, method='nearest')
    uyedge[np.isnan(uyedge)] = uyfill[np.isnan(uyedge)]

    uedge = np.cos(mpsw.angleedge) * uxedge + np.sin(mpsw.angleedge) * uyedge

    return [uedge, hcell]


def run_sw_fwd():
    nml = namelist(nmlname='namelist.sw.x1.10242')
    read_configs(nml)
    read_dims(nml)
    read_vars(nml)

    initial_conditions(nml)
    mpsw.sw_mpas_init_block()

    clock_start, clock_end = clock_namelist(nml)

    today = clock_start
    ulist, vlist, hlist = [], [], []
    unorm = []
    while today <= clock_end:
        if (today - clock_start).total_seconds() % nml.output_interval == 0:
            print(today, mpsw.u[0].min(), mpsw.u[0].max())
            ulist.append(mpsw.ureconstructzonal[0].copy())
            vlist.append(mpsw.ureconstructmeridional[0].copy())
            hlist.append(mpsw.h[0, 0].copy())
            unorm.append(mpsw.u[0, 0].copy())

        mpsw.sw_rk4()

        today += timedelta(seconds=int(mpsw.config_dt))

    # Output your results
    ulist, vlist, hlist = np.array(ulist), np.array(vlist), np.array(hlist)
    unorm = np.array(unorm)

    r2d = 180. / np.pi
    with nc.Dataset(nml.file_output, 'w') as of:
        of.createDimension('nCell', mpsw.latcell.shape[0])
        of.createDimension('nTime', ulist.shape[0])
        of.createVariable('latCell', 'f4', ('nCell'))[:] = mpsw.latcell * r2d
        of.createVariable('lonCell', 'f4', ('nCell'))[:] = mpsw.loncell * r2d
        of.createVariable('ux', 'f8', ('nTime', 'nCell'))[:] = ulist
        of.createVariable('uy', 'f8', ('nTime', 'nCell'))[:] = vlist
        of.createVariable('h', 'f8', ('nTime', 'nCell'))[:] = hlist

        of.createDimension('nEdge', mpsw.latedge.shape[0])
        of.createVariable('latEdge', 'f4', ('nEdge'))[:] = mpsw.latedge * r2d
        of.createVariable('lonEdge', 'f4', ('nEdge'))[:] = mpsw.lonedge * r2d
        of.createVariable('u', 'f8', ('nTime', 'nEdge'))[:] = unorm


def run_sw_tlm(pt=1.E-3):
    nml = namelist(nmlname='namelist.sw.x1.10242')
    read_configs(nml)
    read_dims(nml)
    read_vars(nml)

    initial_conditions(nml)
    mpsw.var_allocation_tlm()
    mpsw.sw_mpas_init_block()

    # initializing TLM variables
    mpsw.u_tl[0] = mpsw.u[0] * pt
    mpsw.h_tl[0] = mpsw.h[0] * pt
    # end initializing TLM variables

    clock_start, clock_end = clock_namelist(nml)

    today = clock_start

    while today <= clock_end:
        mpsw.sw_rk4_tlm()
        today += timedelta(seconds=int(mpsw.config_dt))


def check_tlm():
    nml = namelist(nmlname='namelist.sw.tlmcheck')
    read_configs(nml)
    read_dims(nml)
    read_vars(nml)

    initial_conditions(nml)
    mpsw.var_allocation_tlm()
    mpsw.sw_mpas_init_block()

    idx1, idx2 = 0, 500

    clock_start, clock_end = clock_namelist(nml)

    u_back, h_back = mpsw.u[0].copy(), mpsw.h[0].copy()
    for ipt in range(15):
        pt = 10. ** (-ipt)

        mpsw.u[0] = u_back * (1. + pt)
        mpsw.h[0] = h_back * (1. + pt)
        # ----- running fwd -----
        today = clock_start
        while today <= clock_end:
            mpsw.sw_rk4()
            today += timedelta(seconds=int(mpsw.config_dt))
        u1, h1 = mpsw.u[0, idx1, idx2].copy(), mpsw.h[0, idx1, idx2].copy()
        # ----- end running fwd -----

        # ----- running tlm -----
        mpsw.u[0] = u_back.copy()
        mpsw.h[0] = h_back.copy()
        mpsw.u_tl[0] = mpsw.u[0] * pt
        mpsw.h_tl[0] = mpsw.h[0] * pt

        today = clock_start
        while today <= clock_end:
            mpsw.sw_rk4_tlm()
            today += timedelta(seconds=int(mpsw.config_dt))
        u0, u_tl = mpsw.u[0, idx1, idx2].copy(), mpsw.u_tl[0, idx1, idx2].copy()
        h0, h_tl = mpsw.h[0, idx1, idx2].copy(), mpsw.h_tl[0, idx1, idx2].copy()
        # ----- end running tlm -----

        print(ipt, (h1 - h0) / h_tl, (u1 - u0) / u_tl)


if __name__ == '__main__':
    check_tlm()
    # run_sw_fwd()
    # run_sw_fwd(pt=1.E-2)
    # run_sw_tlm(pt=1.E-2)
