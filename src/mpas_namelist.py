import f90nml as fnl
import json

"""
Purpose:
Parsing MPAS namelist.atmosphere
Note:
For the convenience of package f90nml, output namelist is in the form of dictionaries

Input:
    nmlname: filename of the ascii namelist.atmosphere file
Output:
    namelist(object): in which nlatm is the namelist dictionary

"""

nl_default = \
    {
        "sw_model": {
            "config_test_case": 5,
            "config_time_integration": "RK4",
            "config_dt": 172.8,
            "config_calendar_type": "gregorian_noleap",
            "config_start_time": "0000-01-01_00:00:00",
            "config_stop_time": "none",
            "config_run_duration": "none",
            "config_stats_interval": 100,
            "config_h_ScaleWithMesh": False,
            "config_h_mom_eddy_visc2": 0.0,
            "config_h_mom_eddy_visc4": 0.0,
            "config_h_tracer_eddy_diff2": 0.0,
            "config_h_tracer_eddy_diff4": 0.0,
            "config_thickness_adv_order": 2,
            "config_tracer_adv_order": 2,
            "config_positive_definite": False,
            "config_monotonic": False,
            "config_wind_stress": False,
            "config_bottom_drag": False,
            "config_apvm_upwinding": 0.5,
            "config_num_halos": 2
        },
        "io": {
            "input": "x1.40962.grid.nc",
            "output": "x1.40962.output.nc",
            "output_interval": "06:00:00"
        }
    }


class namelist(object):
    def __init__(self, nmlname='namelist.sw'):
        # ----- reading namelist.atmosphere -----
        # initializing with default values
        self.nlatm = nl_default  # fnl.read(self.nmlname_default)

        # replace default values with the customized
        nl = fnl.read(nmlname)
        for group in nl:
            for config in nl[group]:
                self.nlatm[group][config] = nl[group][config]
        # ----- end reading namelist.atmosphere -----

        self.file_t0 = self.nlatm['io']['input']
        self.file_output = self.nlatm['io']['output']
        output_interval = self.nlatm['io']['output_interval']

        # ----- parsing time format -----
        from pandas import to_timedelta
        self.output_interval = int((to_timedelta(output_interval)).total_seconds())
        # ----- end reading streams.atmosphere -----

    def __str__(self):
        return json.dumps(self.nlatm, indent=4)
