#%%
 
# if True:
#     from ..utils import validation
# else:
#     import validation

import datetime
import numpy as np
import xarray as xr
import pandas as pd
 
# Alias for backward compat
from zutils.xrwrap import *

#%% Some code for adding seasons to xarray datasets

def add_season_custom(ds, month_season_dict, custom_season_name):
    """
    Add custom season to xarray dataset.
    
    Inputs:
        ds - xarray dataset. Time dimension must be called "time".
        month_season_dict - dictionary where the keys are the season names, and the values are the months in that season.
        
    """
    
    allmonths = []
    # Check all 12 months are covered with no overlap
    for season in month_season_dict:
        allmonths += month_season_dict[season]
        
    print(allmonths)
    allmonths = [int(i) for i in allmonths]
    l = len(allmonths)
    
    allmonths = np.unique(np.array(allmonths)).tolist()
    if not len(allmonths)==l:
        raise(Exception('Duplicate months.'))
    
    if len(allmonths)<12:
        raise(Exception('Not all months are covered.'))
    if not min(allmonths) == 1:
        raise(Exception('Invalid months'))
    if not max(allmonths) == 12:
        raise(Exception('Invalid months'))
    
    # Initialise
    custom_season = np.array(['none' for i in ds.time])
    month = ds['time.month']
        
    for season in month_season_dict:
        for m in month_season_dict[season]:
            custom_season[ds['time.month']==m] = season
        
    da = xr.DataArray(custom_season, dims={'time': ds.time})
    ds[custom_season_name] = da
        
def add_season_pilbara_wet_dry(ds):
    """
    Add a Pilbara wet/dry season variable to a dataset. No transition seasons. 
    
    Inputs:
        ds - xarray dataset. Time dimension must be called "time".
        
    """
        
    # Pilbara wet dry with no transitions
    month_season_dict = {'wet': [11, 12, 1, 2, 3, 4],
                         'dry': [5, 6, 7, 8, 9, 10]
    }

    custom_season_name = 'Pilbara wet dry'
    
    add_season_custom(ds, month_season_dict, custom_season_name)
    
def add_season_nw_monsoon(ds):
    """
    Add a northwest monsoon season variable to a dataset. 
    
    Inputs:
        ds - xarray dataset. Time dimension must be called "time".
        
    """
    
    # Pilbara wet dry with no transitions
    month_season_dict = {'NW': [9, 10, 11, 12, 1, 2],
                         't1': [3, 4],
                         'SE': [5, 6, 7],
                         't2': [8]
    }

    custom_season_name = 'NW Monsoon'
    
    add_season_custom(ds, month_season_dict, custom_season_name)
    
# %%
