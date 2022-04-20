import datetime, numpy as np
import zutils.time as ztime
import xarray as xr
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import warnings


def write_temp_prof_file(filename, date, z, T):
    """
    Writes a GOTM temperature profile file

    ** NOTE: format the same as other files also. 
    """
    
    start = date.strftime('%Y/%m/%d %H:%M:%S')
    
    with open(filename, 'w') as f:
        
        f.write('{} {} 2 '.format(start, len(z)) + '\n')
        for z_, T_ in zip(z, T):
            f.write('{} {}'.format(z_, T_) + '\n')
    
    print('Written ' + filename)


def run_gotm(gotm_root='.', gotm_exe='gotm.exe'):
    
    st = os.path.join(gotm_root, gotm_exe)
    ret = os.system(st)
    
    print('System exit status: {}'.format(ret))
    
    # I'm actually not too sure about any of this.  
    if ret==0:
        print('Normal run completion')
    else:
        raise(Exception('Abormal run completion [maybe]'))

    
# def load_gotm(gotm_root='.', gotm_results='GOTM_output.nc'):
    
#     st = os.path.join(gotm_root, gotm_results)
    
# #     ds = load_gotm(gotm_results='GOTM_output.nc' )

# #     read_vars = [v for v in ds.variables.keys() if not v in ['z', 'zi']]
# #     read_vars

#     ds = xr.open_dataset(st, drop_variables=['z', 'zi']) 

#     return ds

def load_gotm(rgnml):
    
    # Just assume NC for now. Could get this from the setup.
    st = os.path.join(rgnml.gotm_root, rgnml.out_fn+'.nc')

    ds = xr.open_dataset(st, drop_variables=['z', 'zi']) 
    ds.load()
    ds.close()
    
    z_coord = np.linspace(-rgnml.depth, 0, rgnml.nlev) 

    ### THIS IS VERY SLOPPY
#     ds['z'].values = z_coord
    warnings.warn('JUST ASSUMING UNIFORM SPACING HERE')
    ds = ds.assign_coords({'z': z_coord})
    
    return ds
    
def make_meteo_file(time, u10, v10, P, T, rh, cc, filename=r'meteo_file.dat', tau_ramp=None, ramp_fcn='exp'):
    """

    function for making the meteo_file.dat file

    file with meteo data (for calc_fluxes=.true.) with
    date: yyyy-mm-dd hh:mm:ss
    x-component of wind (10 m) in m s−1
    y-component of wind (10 m) in m s−1
    air pressure (2 m) in hectopascal
    dry air temperature (2 m) in Celsius
    rel. hum. in % or wet bulb temp. in C or dew point temp. in C
    cloud cover in 1/10
    Example:
    1998-01-01 00:00:00 6.87 10.95 1013.0 6.80 73.2 0.91


    tau_ramp is a ramping timescale
    ramp_fcn is the ramping function ['exp' for exponential or 'pwl' for piecewise linear]
    """
    
    print('u10 is {}'.format(u10))
    
    t_sec = ztime.seconds_since(time, time[0])
    if not tau_ramp is None:
        if ramp_fcn.lower() == 'exp':
            rf = 1-1/np.exp(3*t_sec/tau_ramp)
        elif ramp_fcn.lower() == 'pwl':
            rf = np.ones_like(time).astype(float)
            rf[t_sec<tau_ramp] = t_sec[t_sec<tau_ramp]/tau_ramp
        else:
            raise(Exception('Unrecognised ramp function'))
    else:
        rf = np.ones_like(time).astype(float)
    
    with open(filename, 'w') as f:
        
        for t, r in zip(time, rf):
            
            ts = ztime.datetime64_2_datetime_single(t).strftime('%Y-%m-%d %H:%M:%S')
            st = '{} {} {} {} {} {} {} \n'.format(ts, u10*r, v10*r, P, T, rh, cc) 
            f.write(st)
            
    print('Written ' + filename)
    

class airsea_nml():
    """
    Class for making the airsea.nml file

    swr_method  1: constant value for short wave radiation (const_swr)
                2: read short wave radiation from file
                3: Solar radiation is calculated from time, longitude, latitude,
                and cloud cover.
    """
    
    def __init__(self, gotm_root='.'):
        
        self.gotm_root = gotm_root
        
        self.swr_method = 3
        
    def make(self):
        
        filename = os.path.join(self.gotm_root, 'airsea.nml')
        with open(filename, 'w' ) as f:
            f.writelines("&airsea" + "\n")
            f.writelines("calc_fluxes=     .true."+ "\n")
            f.writelines("fluxes_method = 2"+ "\n")
            f.writelines("back_radiation_method= 4"+ "\n")
            f.writelines("meteo_file=      'meteo_file.dat'"+ "\n")
            f.writelines("hum_method= 1"+ "\n")
            f.writelines("heat_method= 0"+ "\n")
            f.writelines("rain_impact= .false."+ "\n")
            f.writelines("calc_evaporation= .false."+ "\n")
            f.writelines("swr_method = {}".format(self.swr_method)+ "\n")
            f.writelines("const_swr=  0."+ "\n")
            f.writelines("swr_file = 'swrad.dat'"+ "\n")
            f.writelines("swr_factor = 1"+ "\n")
            f.writelines("heatflux_file=   'heatflux.dat'"+ "\n")
            f.writelines("const_heat=      0."+ "\n")
            f.writelines("momentum_method= 0"+ "\n")
            f.writelines("const_tx=        0."+ "\n")
            f.writelines("const_ty=        0."+ "\n")
            f.writelines("momentumflux_file='momentumflux.dat'"+ "\n")
            f.writelines("precip_method=   0"+ "\n")
            f.writelines("const_precip=    2.0e-8"+ "\n")
            f.writelines("precip_file=   'precip.dat'"+ "\n")
            f.writelines("precip_factor = 0.001"+ "\n")
            f.writelines("sst_method=      0"+ "\n")
            f.writelines("sst_file=        'sst.dat'"+ "\n")
            f.writelines("sss_method=      0"+ "\n")
            f.writelines("sss_file=        'sss.dat'"+ "\n")
            f.writelines("/")

        print('Written ' + filename)

class rungotm_nml():
    """
    Class for making the rungotm.nml file


    """
#     !
#     !-------------------------------------------------------------------------------
#     !
#     !-------------------------------------------------------------------------------
#     ! general model setup
#     !
#     ! title            -> title of simulation 
#     ! nlev             -> number of levels
#     ! dt               -> time step in seconds
#     ! cnpar            -> parameter for "explicitness" of numerical scheme
#     !                     (between 0.0 and 1.0)
#     ! buoy_method      -> method to compute mean buoyancy
#     !                  1: from equation of state
#     !                     (i.e. from potential temperature and salinity)
#     !                  2: from prognostic equation
#     !
#     !-------------------------------------------------------------------------------
#      &model_setup
#      /
#     !-------------------------------------------------------------------------------
#     ! geographic location
#     !
#     ! name             -> name of the station
#     ! latitude         -> latitude  in degree (north is positive)
#     ! longitude        -> longitude in degree (east  is positive)
#     ! depth            -> water depth in meters
#     !
#     !-------------------------------------------------------------------------------
#      &station
#      /
#     !
#     !-------------------------------------------------------------------------------
#     ! duration of run
#     !
#     ! timefmt          -> method to specify start and duration of model run
#     !                  1: duration computed from number of time steps, MaxN 
#     !                     (bogus start date used)
#     !                  2: duration computed from given start and stop dates 
#     !                     (number of time steps MaxN computed)
#     !                  3: duration computed from number of time steps, MaxN 
#     !                     (start date as specified, stop date computed)
#     !
#     ! MaxN             -> nominal number of time steps             (see "timefmt")
#     ! start            -> nominal start date: 2018/01/DD HH:01:SS  (see "timefmt")
#     ! stop             -> nominal stop  date: 2018/01/DD HH:01:SS  (see "timefmt")
#     !
#     !-------------------------------------------------------------------------------
#      &time
#      /
#     !-------------------------------------------------------------------------------
#     ! format for output and filename(s).
#     !
#     ! out_fmt          -> format for GOTM output
#     !                  1: ASCII
#     !                  2: NetCDF
#     !                  3: GrADS
#     !
#     ! out_dir          -> path to output directory (set permissions)
#     ! out_fn           -> output string used to generate output file names
#     ! nsave            -> save results every 'nsave' timesteps
#     ! diagnostics      -> diagnostics are written to output (if .true.)
#     !
#     ! mld_method       -> how to diagnose mixed layer depth
#     !                  1: mixed layer depth computed from TKE threshold
#     !                  2: mixed layer depth from Ri threshold
#     ! diff_k           -> TKE threshold [m^2/s^2] for mixed layer depth
#     ! ri_crit          -> Ri threshold for mixed layer depth
#     !
#     ! rad_corr         -> correct surface buoyancy flux for solar radiation
#     !                     for output (if true)
#     !
#     !-------------------------------------------------------------------------------
#      &output
#      /
#     !-------------------------------------------------------------------------------
#     ! Specify variables related to the equation of state.
#     !
#     ! eq_state_mode    -> choice for empirical formula for equation of state 
#     !                  1: UNESCO equation of state by Fofonoff and Millard (1983)
#     !                  2: equation of state according Jackett et al. (2005)
#     !
#     ! eq_state_method  -> method to compute density and buoyancy from salinity,
#     !                     potential temperature and pressure
#     !                  1: full equation of state (i.e. with the LOCAL
#     !                     pressure). This implies that T is NOT treated as
#     !                     the potential temperature but rather as the in-situ
#     !                     temperature!
#     !                  2: equation of state with pressure evaluated at the surface.
#     !                     This implies that T is treated as the potential 
#     !                     temperature and thus rho as the potential density.
#     !                  3: linearized equation of state at T0,S0,p0 
#     !                     (again, use p0=p_surf to work with potential
#     !                     temperature and density.)
#     !                  4: linear equation of state with T0,S0,dtr0,dsr0
#     !
#     ! For the precise definition of the following quantities, see 
#     ! GOTM documentation:
#     !
#     ! T0               -> reference temperature (deg C) for linear equation of state
#     ! S0               -> reference salinity (psu) for linear equation of state
#     ! p0               -> reference pressure (bar) for linear equation of state
#     ! dtr0             -> thermal expansion coefficient for linear equation of state
#     ! dsr0             -> saline expansion coefficient for linear equation of state
#     !-------------------------------------------------------------------------------
#      &eqstate
#      /

    """
    swr_method  1: constant value for short wave radiation (const_swr)
                2: read short wave radiation from file
                3: Solar radiation is calculated from time, longitude, latitude,
                and cloud cover.
    """
    
    def __init__(self, gotm_root='.'):
        
        self.gotm_root = gotm_root
        
        self.title = "GOTM setup made in python zutils"
        self.nlev = 500
        self.dt = 600

        self.name = "Prelude station"
        self.lat = -13.78637
        self.lon = 123.31754
        self.depth = 250.0
        
        self.start = datetime.datetime(2017, 3, 1)
        self.stop = datetime.datetime(2017, 4, 30)
        self.MaxN = 900
        
        self.out_fn = "GOTM_output"
        self.nsave = 6
        
    def get_time_vector(self):
        """
        Just return a time vector spanning the sim to help write other input files
        """
        
        time = np.arange(self.start, self.stop+datetime.timedelta(days=1), datetime.timedelta(seconds=self.dt))
        
        return time
    
    def make(self):
        
        filename = os.path.join(self.gotm_root, 'gotmrun.nml')
        with open(filename, 'w' ) as f:

            start = self.start.strftime('%Y/%m/%d %H:%M:%S')
            stop  = self.stop.strftime('%Y/%m/%d %H:%M:%S')

            f.writelines('&model_setup' + "\n")
            f.writelines(' title=           "{}",'.format(self.title)  + "\n")
            f.writelines(' nlev=              {},'.format(self.nlev) + "\n")
            f.writelines(' dt=                {},'.format(self.dt) + "\n")
            f.writelines(' cnpar=            1.0,' + "\n")
            f.writelines(' buoy_method=         1' + "\n")
            f.writelines('/' + "\n")

            f.writelines('&station' + "\n")
            f.writelines(' name=            "{}",'.format(self.name)  + "\n")
            f.writelines(' latitude=          {},'.format(self.lat)  + "\n")
            f.writelines(' longitude=         {},'.format(self.lon)  + "\n")
            f.writelines(' depth=              {}'.format(self.depth) + "\n")
            f.writelines('/' + "\n")

            f.writelines('&time' + "\n")
            f.writelines(' timefmt=         2,' + "\n")
            f.writelines(' MaxN=            {},'.format(self.MaxN) + "\n")
            f.writelines(' start=           \'{}\','.format(start)  + "\n")
            f.writelines(' stop=            \'{}\''.format(stop)  + "\n")
            f.writelines('/' + "\n")

            f.writelines('&output' + "\n")
            f.writelines(' out_fmt=         2,' + "\n")
            f.writelines(' out_dir=         "{}",'.format(self.gotm_root)  + "\n")
            f.writelines(' out_fn=          "{}",'.format(self.out_fn)  + "\n")
            f.writelines(' nsave=           {},'.format(self.nsave) + "\n")
            f.writelines(' diagnostics=     .false.,' + "\n")
            f.writelines(' mld_method=      1,' + "\n")
            f.writelines(' diff_k=          1.e-5,' + "\n")
            f.writelines(' Ri_crit=         0.5,' + "\n")
            f.writelines(' rad_corr=        .true.' + "\n")
            f.writelines('/' + "\n")

            f.writelines('&eqstate' + "\n")
            f.writelines(' eq_state_mode  = 1,' + "\n")
            f.writelines(' eq_state_method= 1,' + "\n")
            f.writelines(' T0=              10.,' + "\n")
            f.writelines(' S0=              35.,' + "\n")
            f.writelines(' p0=              0.,' + "\n")
            f.writelines(' dtr0=            -0.17,' + "\n")
            f.writelines(' dsr0=            0.78' + "\n")
            f.writelines('/' + "\n")

        print('Written ' + filename)


























####### Test functions

def get_temp_prof_file_ivica():
    
    data = [[-0.921810, 30.073167],
            [-4.136670, 29.950540],
            [-7.427881, 29.897431],
            [-10.801273, 29.852075],
            [-14.265800, 29.811374],
            [-17.833935, 29.782627],
            [-21.521985, 29.751206],
            [-25.350688, 29.707051],
            [-29.345851, 29.633556],
            [-33.539046, 29.520188],
            [-37.968555, 29.353866],
            [-42.680220, 29.124598],
            [-47.728443, 28.823971],
            [-53.177050, 28.451698],
            [-59.100018, 28.003671],
            [-65.581745, 27.500069],
            [-72.716677, 26.953902],
            [-80.607781, 26.331786],
            [-89.363277, 25.531090],
            [-99.091100, 24.476664],
            [-109.889951, 23.119767],
            [-121.836543, 21.589516],
            [-134.968049, 20.077908],
            [-149.260059, 18.740174],
            [-164.601669, 17.671567],
            [-180.770789, 16.819960],
            [-197.416678, 16.050792],
            [-214.058885, 15.335626],
            [-230.114302, 14.835393],
            [-244.960745, 14.256238]]

    return np.array(data)

def get_temp_prof_file_ivica_clipped():
    
    data = [[-0.921810, 29.852075],
            [-4.136670, 29.852075],
            [-7.427881, 29.852075],
            [-10.801273, 29.852075],
            [-14.265800, 29.811374],
            [-17.833935, 29.782627],
            [-21.521985, 29.751206],
            [-25.350688, 29.707051],
            [-29.345851, 29.633556],
            [-33.539046, 29.520188],
            [-37.968555, 29.353866],
            [-42.680220, 29.124598],
            [-47.728443, 28.823971],
            [-53.177050, 28.451698],
            [-59.100018, 28.003671],
            [-65.581745, 27.500069],
            [-72.716677, 26.953902],
            [-80.607781, 26.331786],
            [-89.363277, 25.531090],
            [-99.091100, 24.476664],
            [-109.889951, 23.119767],
            [-121.836543, 21.589516],
            [-134.968049, 20.077908],
            [-149.260059, 18.740174],
            [-164.601669, 17.671567],
            [-180.770789, 16.819960],
            [-197.416678, 16.050792],
            [-214.058885, 15.335626],
            [-230.114302, 14.835393],
            [-244.960745, 14.256238]]

    return np.array(data)