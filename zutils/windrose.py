#%%
 
# if True:
#     from ..utils import validation
# else:
#     import validation

import datetime
import numpy as np
import xarray as xr
import pandas as pd
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import windrose
 
def mainax_2_windax(main_ax, size):
    """
    Simple wrapper to add a windrose axis into an existing axis. 

        size - axis size in cm [assumed square]
    """

    windax = inset_axes(main_ax,
        width=size/2.54,                             # size in inches
        height=size/2.54,                            # size in inches
        loc='center',                        # center bbox at given position
        bbox_to_anchor=(0, 0), # position of the axe
        bbox_transform=main_ax.transData,    # use data coordinate (not axe coordinate)
        axes_class=windrose.WindroseAxes,    # specify the class of the axe
        )

    main_ax.set_xlim([-1, 1])
    main_ax.set_ylim([-1, 1])
    main_ax.axis('off')

    return windax
            
def scale_windax(windax, rl, rd):
    if not rl is None:
        windax.set_rlim((0, rl))
        windax.set_rticks(np.arange(0, rl+1, rd))
        windax.set_yticklabels(['{}%'.format(x) for x in np.arange(0, rl+1, rd)])