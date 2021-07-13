import xarray as xr
import numpy as np
import warnings
import datetime

"""
Module for checking CF compliance of netcdf files or xarray Datasets. 

The conventions are not coded in entirity, and are not coded perfectly. This may lead to both false positives and negatives.
    e.g. 1: Only canonical units are coded, so a height of millimetres (mm) would fail.
    e.g. 2: Some additional Standard Names have been added which are not (yet) in the CF Conventions.

Currently it does not check variables called time for CF compliance as xarray kind of messes with this stuff. Need to think about this.

Variables with the attribute cf_compliant set to False will not be checked for cf_compliance. This is typically 
used for intermediate processing quantities, QAQC flags, and plotting aid tools. 

Data variables used to assist with QAQC, plotting, or processing can take a cf_compliance flag. If this 
is set to boolean false, the variable will not be checked for CF compliance. 

"""

class not_implemented_error(Exception):
    pass

class cf_version_error(Exception):
    pass

class attribute_error(Exception):
    pass

class unit_error(Exception):
    pass

def check_all(file, CF_version):
    """
    Pass a netcdf file and a CF version to check compliance.
    
    If the netcdf file is a path to a file it will be opened as an xarray Dataset. 
    
    If the netcdf can also be passed as an xarray Dataset directly.
    
    Passing xarray DataArrays and netCDF4 objects is not yet implemented. 
    """
    
    if isinstance(file, str):
        print('Filename given, opening dataset')
        ds = xr.open_dataset(file)
    elif isinstance(file, xr.Dataset):
        ds = file
    elif isinstance(file, xr.DataArray):
        raise(not_implemented_error('Passing xarray DataArrays has not yet been implemented'))
              
    check_var_names(ds, CF_version)
    check_global_attrs(ds, CF_version)
    check_var_attrs(ds, CF_version)
    
    print('All checks passed')
    
def check_global_attrs(ds, CF_version):

    CF_version = CF_version.lower()
    
    required_attrs = get_global_attrs(CF_version)
    
    for required_attr in required_attrs:
        if not required_attr in ds.attrs.keys():
            raise(attribute_error('No global attribute "{}"'.format(required_attr)))
        elif ds.attrs[required_attr] in ['?', '', 'blank']:
            warnings.warn('Global attribute "{}" is "{}"'.format(required_attr, ds.attrs[required_attr]))
    
    print('Global attribute checks passed')
    
def check_var_names(ds, CF_version):
    
    CF_version = CF_version.lower()
    
    if CF_version in ['1.7']:
        "There are no particular requirements for variable names in this version."
        pass
    else:
        raise(cf_version_error("CF checker module has not been developed for version {}".format(CF_version)))
        
    print('Variable name checks passed')
    
def check_var_attrs(ds, CF_version):
    """
    Checks that all data variables and coordinates are CF compliant.
    
    NOTE: if a variable has an attribute cf_compliant, and this is False, the variable will not be checked.
    """
    
    stardard_name_to_long_name(ds)
    
    CF_version = CF_version.lower()
    
    required_attrs = get_var_attrs(CF_version)
    
    keys = list(ds.data_vars.keys()) + list(ds.coords.keys())
    
    if CF_version in ['1.7']:
        for key in keys:
            
            if key == 'time':
                print('Variable {} is not being checked for CF compliance'.format(key))
                continue            
            
            if 'cf_compliant' in ds[key].attrs.keys() and not ds[key].attrs['cf_compliant']:
                print('Variable {} is not being checked for CF compliance'.format(key))
                continue
                
            for required_attr in required_attrs:
                if not required_attr in ds[key].attrs.keys():
                    raise(attribute_error('No attribute "{}" in variable {}'.format(required_attr, key))) 
                
                
                if ds[key].attrs[required_attr] in ['?', '', 'blank']:
                    warnings.warn('Attribute {} in variable {} is {}'.format(required_attr, key, ds[key].attrs[required_attr]))
         
            check_standard_names_and_canonical_units(ds[key])
    else:
        raise(cf_version_error("CF checker module has not been developed for version {}".format(CF_version)))
        
    print('Variable name checks passed')

    ds.coords
    
def check_standard_names_and_canonical_units(da):
    """
    For a given DataArray we check that the standard name exists within a subset standard names defined in this module. 
    We will also check whether the units are equal to or at least equivalent to canonical units.
        
    """
    
    standard_name = da.attrs['standard_name']
    units = da.attrs['units']
    
    canonical_units = get_canonical_units()

    if standard_name in canonical_units.keys():
        if not units in canonical_units[standard_name]:
            raise(Exception('{} is not a recognised unit for {} [Check use of Capitals]'.format(units, standard_name)))
        else:
            return
    else:
        raise(Exception('{} is not a recognised standard_name'.format(standard_name)))    

    raise(Exception('{} is not a recognised standard_name'.format(standard_name)))
    
"""
GET_* methods. Generally return lists.
"""
def get_global_attrs(CF_version):
    
    if CF_version == '1.7':
        required_attrs = ['Conventions', 'title', 'institution', 'source', 'history', 'references', 'comment']
    else:
        raise(cf_version_error("CF checker module has not been developed for version {}".format(CF_version)))
        
    return required_attrs

def get_var_attrs(CF_version):
    
    if CF_version in ['1.7']:
        required_attrs = ['units', 'long_name', 'standard_name']
    else:
        raise(cf_version_error("CF checker module has not been developed for version {}".format(CF_version)))
        
    return required_attrs

def get_canonical_units():
    """
    NOTE: THESE ARE NOT REAL CF STANDARD NAMES. MANY ARE MADE UP. 

    """

    misc_dict = {
                    'signal_intensity_from_multibeam_acoustic_doppler_velocity_sensor_in_sea_water': ['counts', 'db'], # Echo/bachscatter
                    'proportion_of_acceptable_signal_returns_from_acoustic_instrument_in_sea_water': ['%', 'percent'], # percent_good
                    'indicative_error_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water': ['m/s'],
                    'beam_consistency_indicator_from_multibeam_acoustic_doppler_velocity_profiler_in_sea_water': ['%'],
                    'acoustic_signal_roundtrip_travel_time_in_sea_water': ['s'],
                    'acoustic_surface_track_distance': ['m'], # Not a CF name
                    'acoustic_bottom_track_distance': ['m'], # Not a CF name
                    'frequency': ['Hz', 'rad/s'], # Not a CF name
                    'wave_frequency': ['Hz', 'rad/s'],  # Not a CF name
                    'sea_surface_wave_frequency': ['Hz', 'rad/s'],  # Not a CF name
                    'turbulence_frequency': ['Hz', 'rad/s'],  # Not a CF name
                    'wavenumber': ['cpm', 'rad/m'], # Not a CF name
                    'wave_wavenumber': ['cpm', 'rad/m'],  # Not a CF name
                    'sea_surface_wave_wavenumber': ['cpm', 'rad/m'],  # Not a CF name
                    'turbulence_wavenumber': ['cpm', 'rad/m'],  # Not a CF name
                    'variance_spectral_density': ['m2 Hz-1'],
                    'wave_variance_spectral_density': ['m2 Hz-1'],
                    'sea_surface_wave_variance_spectral_density': ['m2 Hz-1'],
                    'directional_variance_spectral_density': ['m2 Hz-1 deg-1'],
                    'sea_surface_wave_directional_variance_spectral_density': ['m2 Hz-1 deg-1'],
                    'directional_wave_variance_spectral_density': ['m2 Hz-1 deg-1'],
                    'directional_sea_surface_wave_variance_spectral_density': ['m2 Hz-1 deg-1'],
       }

    temperature_list = [
                    'temperature', 
                    'air_temperature', 
                    'seawater_temperature'
                    ]

    temperature_dict = {}
    for i in temperature_list:
            temperature_dict[i] = ['degree', 'deg', 'Kelvin', 'K']

    pressure_list = [
                    'pressure',
                    'sea_water_pressure', 
                    'sea_water_pressure_due_to_sea_water'
                    ]

    pressure_dict = {}
    for i in pressure_list:
            pressure_dict[i] = ['dbar']
            
    velocity_list = [
                    'velocity', 
                    'speed', 
                    'seawater_velocity',
                    'seawater_speed',
                    'eastward_seawater_velocity',
                    'northward_seawater_velocity',
                    'upward_seawater_velocity',
                    'wind_velocity',
                    'wind_speed',
                    'easterly_wind_velocity',
                    'radial_sea_water_velocity_away_from_instrument',
                    'radial_sea_water_velocity_toward_instrument',
                    'radial_velocity_of_scatterers_away_from_instrument',
                    'radial_velocity_of_scatterers_toward_instrument']
    
    velocity_dict = {}
    for i in velocity_list:
            velocity_dict[i] = ['m/s', 'm s-1']

    direction_list = ['seawater_current_to_direction', ]
    
    direction_dict = {}
    for i in direction_list:
            direction_dict[i] = ['degree', 'deg']

    depth_list = ['thermocline_depth_below_sea_surface',
                  'pressure_sensor_depth_below_sea_surface',
                  'sensor_depth_below_sea_surface'
                  'sea_floor_depth_below_sea_surface']

    depth_dict = {}
    for i in depth_list:
            depth_dict[i] = ['m']
    
    distance_list = ['distance',
                  'distance_from_instrument',
                  'distance_from_seabed',
                  'distance_from_surface']
    
    distance_dict = {}
    for i in distance_list:
            distance_dict[i] = ['m']

    surface_wave_types = ['sea_surface_wave', 
                  'sea_surface_wind_wave', 
                  'sea_surface_swell_wave', 
                  'sea_surface_primary_swell_wave', 
                  'sea_surface_secondary_swell_wave', 
                  'sea_surface_tertiary_swell_wave', 
                  'sea_surface_infra_gravity_wave']
    
    ## DEFINED FOR EACH WAVE TYPE
    surface_wave_direction_calc_list = [
        'directional_spread',
        'directional_spread_at_variance_spectral_density_maximum',
        'energy_at_variance_spectral_density_maximum',
        'from_direction',
        'from_mean_direction',
        'from_direction_at_variance_spectral_density_maximum', # Peak period
    ]
    
    surface_wave_height_calc_list = [
            'significant_height',
            'mean_height',
            'mean_height_of_highest_third',
            'mean_height_of_highest_tenth',
            'maximum_height',
            'maximum_crest_height',
            'maximum_trough_depth',
    ]
    
    surface_wave_period_calc_list = [
        'mean_period',
        'mean_period_from_variance_spectral_density_first_frequency_moment',
        'mean_period_from_variance_spectral_density_inverse_frequency_moment',
        'mean_period_from_variance_spectral_density_second_frequency_moment',
        'mean_period_of_highest_third',
        'mean_period_of_highest_tenth',
        'significant_period',
        'period_at_variance_spectral_density_maximum',
        'period_of_highest_wave',
        'maximum_period',
    ]

    surface_wave_slope_calc_list = [
        'maximum_steepness',
    ]
    
    wavedirection_dict = {}
    for wave in surface_wave_types:
        for calc in surface_wave_direction_calc_list:
            wavedirection_dict[wave + '_' + calc] = ['degree', 'deg']
            
    waveheight_dict = {}
    for wave in surface_wave_types:
        for calc in surface_wave_height_calc_list:
            waveheight_dict[wave + '_' + calc] = ['m']
            
    waveperiod_dict = {}
    for wave in surface_wave_types:
        for calc in surface_wave_period_calc_list:
            waveperiod_dict[wave + '_' + calc] = ['s']
            
    waveslope_dict = {}
    for wave in surface_wave_types:
        for calc in surface_wave_slope_calc_list:
            waveslope_dict[wave + '_' + calc] = ['1']
    
    out = {}
    out.update(misc_dict)
    out.update(temperature_dict)
    out.update(pressure_dict)
    out.update(velocity_dict)
    out.update(direction_dict)
    out.update(depth_dict)
    out.update(distance_dict)
    out.update(wavedirection_dict)
    out.update(waveheight_dict)
    out.update(waveperiod_dict)
    out.update(waveslope_dict)

    return out

"""
FORCE_*_COMPLIANCE methods. Prepare Netcdfs for CF compliance. Generally add the members of GET Method lists to some other lists.
"""
def force_compliance(ds, CF_version):
    
    pass
    
def force_global_attrs_compliance(ds, CF_version):
    
    required_attrs = get_global_attrs(CF_version)
    
    for required_attr in required_attrs:
        if not required_attr in ds.attrs.keys():
            ds.attrs[required_attr] = '?'
            
    ds.attrs['Conventions'] = 'CF-' + CF_version
        
def force_var_attrs_compliance(ds, CF_version):
    
    required_attrs = get_var_attrs(CF_version)
    
    for data_var in ds.data_vars:
        for required_attr in required_attrs:
            if not required_attr in ds[data_var].attrs.keys():
                ds[data_var].attrs[required_attr] = '?'
                    
"""
Add methods. Append timestamped strings to attributes.
"""
def add_history(ds, author, string):
        
    attr = 'history'
    add_string(ds, attr, author, string)
        
def add_comment(ds, author, string):
        
    attr = 'comment'
    add_string(ds, attr, author, string)
        
def add_string(ds, attr, author, string):
    
    new_string = datetime.datetime.now().isoformat() + ': ' + '[{}]'.format(author) + ' ' + string
    
    if not attr in ds.attrs.keys():
        ds.attrs[attr] = new_string
    elif ds.attrs[attr] in ['', '?', 'blank']:
        ds.attrs[attr] = new_string
    else:
        ds.attrs[attr] = ds.attrs[attr] + ';' + new_string
        
    return new_string

"""
Misc methods
"""
def stardard_name_to_long_name(ds):
    """
    Every variable with a standard name but no long_name will have the long name set to the standard name.
    """

    keys = list(ds.data_vars.keys()) + list(ds.coords.keys())
    
    for key in keys:
        if 'standard_name' in ds[key].attrs:
            if not 'long_name' in ds[key].attrs:
                ds[key].attrs['long_name'] = ds[key].attrs['standard_name']
            elif  ds[key].attrs['long_name'] in ['', '?', 'blank']:
                ds[key].attrs['long_name'] = ds[key].attrs['standard_name']
