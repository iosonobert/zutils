import pandas as pd, xarray as xr
import numpy as np, matplotlib.pyplot as plt, datetime
import scipy.signal
import zutils.nan
import zutils.time as ztime
import zutils.tidal as ztide

def dirmath_to_dirgeo(dirmath, degrees=False):
    """
    Mathematical direction to geographical. Direction is in radians. 
    Direction is radians unless set degrees=True
    """
    
    return dirgeo_to_dirmath(dirmath, degrees=degrees)

def dirgeo_to_dirmath(dirgeo, degrees=False):
    """
    Mathematical direction to geographical. 
    Direction is radians unless set degrees=True
    """
    
    if degrees: 
        dirmath = 90 - dirgeo
    else: 
        dirmath = np.pi/2 - dirgeo
        
    return dirmath

def flip_dir(dire, degrees=False):
    """
    Flip direction 180 degrees/pi radians.
    Direction is radians unless set degrees=True
    """
    
    if degrees: 
        dire = 180 + dire
        dire = np.mod(dire, 360)
    else: 
        dire = np.pi + dire
        dire = np.mod(dire, 2*np.pi)
        
    return dire
    
def sd_to_uv(spd, dirgeo, convention, degrees=False):
    """
    Convert speed and direction to U and V. 
    Direction is clockwise from North. 
    Direction is radians unless set degrees=True
    """
        
    if not convention.lower() in ['meteo', 'ocean']:
        raise(Exception("Convention must be 'meteo' or 'ocean'"))
        
    if convention.lower() == 'meteo':
        dirgeo = flip_dir(dirgeo, degrees=degrees)
    
    dirmath = dirgeo_to_dirmath(dirgeo, degrees=degrees)
    
    if degrees:
        dirmath = dirmath*np.pi/180
        
    u = spd * np.cos(dirmath)
    v = spd * np.sin(dirmath)
    
    return u, v

def uv_to_sd(u, v, convention, degrees=False):
    """
    Convert U and V to speed and direction. 
    Direction is clockwise from North. 
    Direction is radians unless set degrees=True
    """
        
    if not convention.lower() in ['meteo', 'ocean']:
        raise(Exception("Convention must be 'meteo' or 'ocean'"))
        
    dirmath = np.arctan2(v, u)
    dirmath = np.mod(dirmath, np.pi*2)
    if degrees:
        dirmath = dirmath*180/np.pi
        
    spd = np.sqrt(u**2 + v**2)
    
    dirgeo = dirmath_to_dirgeo(dirmath, degrees=degrees)
    if convention.lower() == 'meteo':
        dirgeo = flip_dir(dirgeo, degrees=degrees)
        
    return spd, dirgeo
