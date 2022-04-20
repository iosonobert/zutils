import os, pandas as pd, matplotlib.pyplot as plt, xarray as xr, datetime, sfoda
import sfoda.tides.readotps as stides, numpy as np

package_directory = os.path.dirname(os.path.abspath(__file__))

modfile_default = os.path.join(package_directory, r'..\data\tpxo\REGIONAL_SOLNS\IO\Model_IO')

print(modfile_default)

def quick_tidal_pred(lon, lat, time, modfile=modfile_default):
    """
    Very quick tidal prediction code using SFODA ontop of TPXO binaries.

    Rerurns:
        h 
        u 
        v 
    """

    pathfile = os.path.split(modfile)
    path = pathfile[0]

    with open(modfile,'r') as f:
        hfile = path+'/' + f.readline().strip()
        uvfile = path+'/' + f.readline().strip()
        grdfile = path+'/' + f.readline().strip() 

    h, u, v = stides.tide_pred(modfile, lon, lat, time)

    return h, u, v