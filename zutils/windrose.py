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

import matplotlib as mpl
import locale

 
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


def pseudolegend(legend_ax, wind_axis, decimal_places=1, slope=0.5, occupy_x=1, text_args={}, **kwargs):
    """
    THIS IS BASICALLY JUST A TWEAK OF THE LEGEND FUNCTION FROM THE ORIGINAL WINDROSE MODULE. 

    Turns a matplotlib into a PSEUDOLEGEND. This is called a 'pseudolegend' because it only imitatesactually  a 
    matplotlib legend, it does not geterate a legend.
    
    Parameters
    ----------
    legend_ax : this is the axis you'll be adding the legend into
    wind_ax : windrose.windrose.WindroseAxes object to get the info from. 
    decimal_places : int, default 1
        The decimal places of the formated legend
    title: str, default Noned  
        Title of the legend - include units
    occupy_x : int, default 1
        Controls how much of the x axis is taken up with the plot. 
    
    Other Parameters
    ----------------
    slope : float, Default = 0.5
        This is the slope that the axis is plotted at
    
    """

    def get_handles():
        handles = list()
        for p in wind_axis.patches_list:
            if isinstance(p, mpl.patches.Polygon) or isinstance(
                p, mpl.patches.Rectangle
            ):
                color = p.get_facecolor()
            elif isinstance(p, mpl.lines.Line2D):
                color = p.get_color()
            else:
                raise AttributeError("Can't handle patches")
            handles.append(
                mpl.patches.Rectangle(
                    (0, 0), 0.2, 0.2, facecolor=color, edgecolor="black"
                )
            )
        return handles

    def get_labels(decimal_places=1):
        _decimal_places = str(decimal_places)

        fmt = "%." + _decimal_places + "f " + "-%0." + _decimal_places + "f"

        labels = np.copy(wind_axis._info["bins"])
        if locale.getlocale()[0] in ["fr_FR"]:
            fmt += "["
        else:
            fmt += ""

        labels_ = [fmt % (labels[i], labels[i + 1]) for i in range(len(labels) - 1)]
        
        if '-inf' in labels_[-1]:
            fmt = ">%." + _decimal_places + "f " 
            labels_[-1] = fmt % (labels[-2])
        
        labels = labels_
        
        return labels

    kwargs.pop("labels", None)
    kwargs.pop("handles", None)
    
    title = kwargs.pop('title', None)

    handles = get_handles()
    labels = get_labels(decimal_places)

    slope *= occupy_x # This keeps the slope scaling even if you ajust occupy x
    n = len(handles)

    for i, handle in enumerate(handles):
        x = np.array([0, 0, 1, 1]) + i
        y = np.array([-1, 1, 1, -1]) + (np.array([-i, 0+i, 1+i, -1-i]))*slope

#         h, = legend_ax.fill(x, y, edgecolor=handle.get_edgecolor())
        h, = legend_ax.fill(x, y, edgecolor='w')
        h.set_facecolor(handle.get_facecolor())

        legend_ax.text(i+0.5, -(1+slope*(n+1)), labels[i], horizontalalignment='center', verticalalignment='top', **text_args)

        
        
#         legend_ax.set_xlim(0, n)
        legend_ax.set_title(title)
      
    dx = n*(1+(1-occupy_x))-n
    xl = (-dx, n+dx)
    print(xl)
    legend_ax.set_xlim(xl)
    legend_ax.axis('off')

    legend_ax.set_ylim((-(1+slope*n+0.1)/occupy_x, (1+slope*n+0.2)/occupy_x))

    return handles, labels
    