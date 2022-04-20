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

# from windrose import WindroseAxes # duplicated below now
# import windrose # duplicated below now

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
        axes_class=WindroseAxes,    # specify the class of the axe
        )

    main_ax.set_xlim([-1, 1])
    main_ax.set_ylim([-1, 1])
    main_ax.axis('off')

    return windax
            
def scale_windax(windax, rl, rd, offset=0, horizontalalignment='center', verticalalignment='center'):


    if not rl is None:
        windax.set_rlim((0, offset+rl))
        windax.set_rticks(np.arange(offset+rd, offset+rl+1, rd))

        windax.set_yticklabels(['{}%'.format(x) for x in np.arange(rd, rl+1, rd)], 
                    horizontalalignment=horizontalalignment,
                    verticalalignment=verticalalignment) 


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



###########################################################
## Below is a duplication of some parts of the windrose library for testing before sending pullrequests to the admins

import locale
import matplotlib as mpl
from matplotlib import docstring
import numpy as np
import random
from matplotlib.projections.polar import PolarAxes
from numpy.lib.twodim_base import histogram2d
import matplotlib.pyplot as plt

ZBASE = -1000  # The starting zorder for all drawing, negative to have the grid on
# ZBASE = 1000  # The starting zorder for all drawing, negative to have the grid on
VAR_DEFAULT = "speed"
DIR_DEFAULT = "direction"
FIGSIZE_DEFAULT = (8, 8)
DPI_DEFAULT = 80

CALM_CIRCLE_EDGE_COLOR = "black"
CALM_CIRCLE_FILL_COLOR = "blue"
CALM_CIRCLE_ALPHA = 0.1

CALM_CIRCLE_EDGE_COLOR = "black"
CALM_CIRCLE_FILL_COLOR = "gray"
CALM_CIRCLE_ALPHA = 1

class WindAxesFactory(object):
    """
    Factory class to create WindroseAxes or WindAxes
    """

    @staticmethod
    def create(typ, ax=None, *args, **kwargs):
        """
        Create
        Mandatory:
        Parameters
        ----------
        typ : string, 'windroseaxes' or 'windaxes'
            Type of axes to create
                * windroseaxes : a WindroseAxes axe
                * windaxe : a WindAxes axe
        ax : matplotlib.Axes, optional
            A matplotlib axe
        """
        typ = typ.lower()
        d = {"windroseaxes": WindroseAxes, "windaxes": WindAxes}
        if typ in d.keys():
            cls = d[typ]
            if isinstance(ax, cls):
                return ax
            else:
                ax = cls.from_ax(ax, *args, **kwargs)
                return ax
        else:
            raise NotImplementedError("typ=%r but it might be in %s" % (typ, d.keys()))


class WindroseAxes(PolarAxes):
    """
    Create a windrose axes
    """

    name = "windrose"

    def __init__(self, *args, **kwargs):
        """
        See Axes base class for args and kwargs documentation
        """

        # Uncomment to have the possibility to change the resolution directly
        # when the instance is created
        # self.RESOLUTION = kwargs.pop('resolution', 100)
        self.rmax = kwargs.pop("rmax", None)
        
        
        self.theta_labels = kwargs.pop("theta_labels", None)
        if self.theta_labels is None:
            self.theta_labels = ["E", "N-E", "N", "N-W", "W", "S-W", "S", "S-E"]
        PolarAxes.__init__(self, *args, **kwargs)
        self.set_aspect("equal", adjustable="box", anchor="C")
        self.radii_angle = 67.5
        self.cla()

    @staticmethod
    def from_ax(ax=None, fig=None, rmax=None, theta_labels=None, rect=None, *args, **kwargs):
        """
        Return a WindroseAxes object for the figure `fig`.
        """
        if ax is None:
            if fig is None:
                fig = plt.figure(
                    figsize=FIGSIZE_DEFAULT,
                    dpi=DPI_DEFAULT,
                    facecolor="w",
                    edgecolor="w",
                )
            if rect is None:
                rect = [0.1, 0.1, 0.8, 0.8]
            ax = WindroseAxes(fig, rect, rmax=rmax, theta_labels=theta_labels, *args, **kwargs)
            fig.add_axes(ax)
            return ax
        else:
            return ax

    def cla(self):
        """
        Clear the current axes
        """
        PolarAxes.cla(self)

        self.theta_angles = np.arange(0, 360, 45)
        self.set_thetagrids(angles=self.theta_angles, labels=self.theta_labels)

        self._info = {"dir": list(), "bins": list(), "table": list()}

        self.patches_list = list()

        self.calm_count = None

    def _colors(self, cmap, n):
        """
        Returns a list of n colors based on the colormap cmap
        """
        return [cmap(i) for i in np.linspace(0.0, 1.0, n)]

    def set_radii_angle(self, **kwargs):
        """
        Set the radii labels angle
        """

        kwargs.pop("labels", None)
        angle = kwargs.pop("angle", None)
        if angle is None:
            angle = self.radii_angle
        self.radii_angle = angle
        N = 5
        rmax = self.get_rmax()
        radii = np.linspace(0, rmax, N + 1)
        if rmax % N == 0:
            fmt = "%d"
        else:
            fmt = "%.1f"
        radii_labels = [fmt % r for r in radii]
        # radii_labels[0] = ""  # Removing label 0
        self.set_rgrids(
            radii=radii[1:], labels=radii_labels[1:], angle=self.radii_angle, **kwargs
        )

    def _update(self):
        if not self.rmax:
            self.rmax = np.max(np.sum(self._info["table"], axis=0))
        calm_count = self.calm_count or 0
        self.set_rmax(rmax=self.rmax + calm_count)
        self.set_radii_angle(angle=self.radii_angle)

    def legend(self, loc="lower left", decimal_places=1, units=None, **kwargs):
        """
        Sets the legend location and her properties.
        Parameters
        ----------
        loc : int, string or pair of floats, default: 'lower left'
            see :obj:`matplotlib.pyplot.legend`.
        decimal_places : int, default 1
            The decimal places of the formated legend
        units: str, default None
        Other Parameters
        ----------------
        isaxes : boolean, default True
            whether this is an axes legend
        prop : FontProperties(size='smaller')
            the font property
        borderpad : float
            the fractional whitespace inside the legend border
        shadow : boolean
            if True, draw a shadow behind legend
        labelspacing : float, 0.005
            the vertical space between the legend entries
        handlelenght : float, 0.05
            the length of the legend lines
        handletextsep : float, 0.02
            the space between the legend line and legend text
        borderaxespad : float, 0.02
            the border between the axes and legend edge
        kwarg
            Every other kwarg argument supported by
            :obj:`matplotlib.pyplot.legend`
        """

        def get_handles():
            handles = list()
            for p in self.patches_list:
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

        def get_labels(decimal_places=1, units=None):
            _decimal_places = str(decimal_places)

            fmt = "[%." + _decimal_places + "f " + ": %0." + _decimal_places + "f"

            labels = np.copy(self._info["bins"])
            if locale.getlocale()[0] in ["fr_FR"]:
                fmt += "["
            else:
                fmt += ")"

            if units:
                fmt += ' ' + units

            labels = [fmt % (labels[i], labels[i + 1]) for i in range(len(labels) - 1)]
            return labels

        kwargs.pop("labels", None)
        kwargs.pop("handles", None)

        # decimal_places = kwargs.pop('decimal_places', 1)

        handles = get_handles()
        labels = get_labels(decimal_places, units)
        self.legend_ = mpl.legend.Legend(self, handles, labels, loc, **kwargs)
        return self.legend_

    def set_legend(self, **pyplot_arguments):
        if "borderaxespad" not in pyplot_arguments:
            pyplot_arguments["borderaxespad"] = -0.10
        legend = self.legend(**pyplot_arguments)
        plt.setp(legend.get_texts(), fontsize=8)
        return legend

    def _init_plot(self, direction, var, calm_var, **kwargs):
        """
        Internal method used by all plotting commands
        Parameters
        ----------
        direction : 1D array,
            directions the wind blows from, North centred
        var : 1D array,
            values of the variable to compute. Typically the wind speeds
        calm_var : 1D array
            values of the variable to which the calms are considered against. 
            Default is the same as the var variable. Use another variable if calms
            are to be judged by a variable other than the plotted variable. 
        Other Parameters
        ----------------
        normed : boolean, default False
        blowto : boolean, default False
        colors : str or list of str, default None
            The colors of the plot.
        cmap : color map, default `jet`
            A :obj:`matplotlib.cm` colormap for the plot.
            Warning! It overrides `colors`.
        weibull_factors :
        mean_values :
        frequency :
        calm_limit : float, default None
        kwarg
            Any argument accepted by :obj:`matplotlib.pyplot.plot`.
        """

        # if weibull factors are entered overwrite direction and var
        if "weibull_factors" in kwargs or "mean_values" in kwargs:
            if "weibull_factors" in kwargs and "mean_values" in kwargs:
                raise TypeError("cannot specify both weibull_factors and mean_values")
            statistic_type = "unset"
            if "weibull_factors" in kwargs:
                statistic_type = "weibull"
                val = kwargs.pop("weibull_factors")
            elif "mean_values" in kwargs:
                statistic_type = "mean"
                val = kwargs.pop("mean_values")
            if val:
                if "frequency" not in kwargs:
                    raise TypeError(
                        "specify 'frequency' argument for statistical input"
                    )
                windFrequencies = kwargs.pop("frequency")
                if len(windFrequencies) != len(direction) or len(direction) != len(var):
                    if len(windFrequencies) != len(direction):
                        raise TypeError("len(frequency) != len(direction)")
                    elif len(direction) != len(var):
                        raise TypeError("len(frequency) != len(direction)")
                windSpeeds = []
                windDirections = []
                for dbin in range(len(direction)):
                    for _ in range(int(windFrequencies[dbin] * 10000)):
                        if statistic_type == "weibull":
                            windSpeeds.append(
                                random.weibullvariate(var[dbin][0], var[dbin][1])
                            )
                        elif statistic_type == "mean":
                            windSpeeds.append(
                                random.weibullvariate(var[dbin] * 2 / np.sqrt(np.pi), 2)
                            )
                        windDirections.append(direction[dbin])
                var, direction = windSpeeds, windDirections

        # self.cla()
        kwargs.pop("zorder", None)

        # Init of the bins array if not set
        bins = kwargs.pop("bins", None)
        if bins is None:
            bins = np.linspace(np.min(var), np.max(var), 6)
        if isinstance(bins, int):
            bins = np.linspace(np.min(var), np.max(var), bins)
        bins = np.asarray(bins)
        nbins = len(bins)

        # Number of sectors
        nsector = kwargs.pop("nsector", None)
        if nsector is None:
            nsector = 16

        # Sets the colors table based on the colormap or the "colors" argument
        colors = kwargs.pop("colors", None)
        cmap = kwargs.pop("cmap", None)
        if colors is not None:
            if isinstance(colors, str):
                colors = [colors] * nbins
            if isinstance(colors, (tuple, list)):
                if len(colors) != nbins:
                    raise ValueError("colors and bins must have same length")
        else:
            if cmap is None:
                cmap = mpl.cm.jet
            colors = self._colors(cmap, nbins)

        # Building the angles list
        angles = np.arange(0, -2 * np.pi, -2 * np.pi / nsector) + np.pi / 2

        normed = kwargs.pop("normed", False)
        blowto = kwargs.pop("blowto", False)

        # Calm condition
        calm_limit = kwargs.pop("calm_limit", None)
        if calm_limit is not None:
            mask = calm_var > calm_limit
            self.calm_count = len(calm_var) - np.count_nonzero(mask)
            if normed:
                self.calm_count = self.calm_count * 100 / len(calm_var)
            var = var[mask]
            direction = direction[mask]

        # Set the global information dictionnary
        self._info["dir"], self._info["bins"], self._info["table"] = histogram(
            direction, var, bins, nsector, normed, blowto
        )

        return bins, nbins, nsector, colors, angles, kwargs

    def _calm_circle(self, default_offset=0):
        """
        Draw the calm centered circle
        and return the initial offset for plots methods
        """
        if self.calm_count and self.calm_count > 0:
            if default_offset==0:
                ccr=self.calm_count
            else:
                ccr = default_offset

            zorder = abs(ZBASE)**2 
            # zorder = ZBASE

            circle = mpl.patches.Circle(
                (0., 0.),
                ccr,
                transform=self.transData._b,
                facecolor=CALM_CIRCLE_FILL_COLOR,
                edgecolor=CALM_CIRCLE_EDGE_COLOR,
                alpha=CALM_CIRCLE_ALPHA,
                zorder=zorder,
                )

            self.add_artist(circle)
            plt.text(0, 0, 'calm', ha='center', va='bottom', zorder=zorder+1)
            plt.text(0, 0, '{:.0f}%'.format(self.calm_count), ha='center', va='top', zorder=zorder+1)
            
        elif default_offset>0:
            ccr = default_offset
            circle = mpl.patches.Circle(
                (0., 0.),
                ccr,
                transform=self.transData._b,
                facecolor='white',
                edgecolor=CALM_CIRCLE_EDGE_COLOR,
                alpha=1,
                zorder=abs(ZBASE)**2,
                )

            self.add_artist(circle)
        
        else: 
            ccr=0

        # return self.calm_count or default_offset
        return ccr


    def bar(self, direction, var, calm_var=None, **kwargs):
        """
        Plot a windrose in bar mode. For each var bins and for each sector,
        a colored bar will be draw on the axes.
        Parameters
        ----------
        direction : 1D array
            directions the wind blows from, North centred
        var : 1D array
            values of the variable to compute. Typically the wind speeds.
        calm_var : 1D array
            values of the variable to which the calms are considered against. 
            Default is the same as the var variable. Use another variable if calms
            are to be judged by a variable other than the plotted variable. 
        Other Parameters
        ----------------
        nsector : integer, optional
            number of sectors used to compute the windrose table. If not set,
            nsectors=16, then each sector will be 360/16=22.5°, and the
            resulting computed table will be aligned with the cardinals points.
        bins : 1D array or integer, optional
            number of bins, or a sequence of bins variable. If not set, bins=6
            between min(`var`) and max(`var`).
        blowto : bool, optional.
            if True, the windrose will be pi rotated, to show where the wind
            blow to (usefull for pollutant rose).
        colors : string or tuple, optional
            one string color ('k' or 'black'), in this case all bins will be
            plotted in this color; a tuple of matplotlib color args (string,
            float, rgb, etc), different levels will be plotted
            in different colors in the order specified.
        cmap : a cm Colormap instance from :obj:`matplotlib.cm`, optional.
            if cmap == None and colors == None, a default Colormap is used.
        edgecolor : string, optional
            The string color each edge box will be plotted.
            Default : no edgecolor
        opening : float, optional
            between 0.0 and 1.0, to control the space between each sector (1.0
            for no space)
        calm_limit : float, optional
            Calm limit for the var parameter. If not None, a centered red
            circle will be draw for representing the calms occurences and all
            datas below this value will be removed from the computation.
        """

        if calm_var is None:
            calm_var = var.copy()

        bins, nbins, nsector, colors, angles, kwargs = self._init_plot(
            direction, var, calm_var, **kwargs
        )
        kwargs.pop("facecolor", None)
        edgecolor = kwargs.pop("edgecolor", None)
        default_offset = kwargs.pop("default_offset", 0)
        
        

        if edgecolor is not None:
            if not isinstance(edgecolor, str):
                raise ValueError("edgecolor must be a string color")
        opening = kwargs.pop("opening", None)
        if opening is None:
            opening = 0.8
        dtheta = 2 * np.pi / nsector
        opening = dtheta * opening

        offs = self._calm_circle(default_offset=default_offset)
        print('Offset is below')
        print(offs)
        print('Offset is above')
        
        for j in range(nsector):
            offset = offs
            for i in range(nbins):
                if i > 0:
                    offset += self._info["table"][i - 1, j]
                val = self._info["table"][i, j]
                zorder = ZBASE + nbins - i
                patch = mpl.patches.Rectangle(
                    (angles[j] - opening / 2, offset),
                    opening,
                    val,
                    facecolor=colors[i],
                    edgecolor=edgecolor,
                    zorder=zorder,
                    clip_on=False,
                    **kwargs
                )
                self.add_patch(patch)
                if j == 0:
                    self.patches_list.append(patch)
        self._update()
        
        return offs


def histogram(direction, var, bins, nsector, normed=False, blowto=False):
    """
    Returns an array where, for each sector of wind
    (centred on the north), we have the number of time the wind comes with a
    particular var (speed, polluant concentration, ...).
    Parameters
    ----------
    direction : 1D array
        directions the wind blows from, North centred
    var : 1D array
        values of the variable to compute. Typically the wind speeds
    bins : list
        list of var category against we're going to compute the table
    nsector : integer
        number of sectors
    Other Parameters
    ----------------
    normed : boolean, default False
        The resulting table is normed in percent or not.
    blowto : boolean, default False
        Normaly a windrose is computed with directions as wind blows from. If
        true, the table will be reversed (usefull for pollutantrose)
    """

    if len(var) != len(direction):
        raise ValueError("var and direction must have same length")

    angle = 360. / nsector

    dir_bins = np.arange(-angle / 2, 360. + angle, angle, dtype=float)
    dir_edges = dir_bins.tolist()
    dir_edges.pop(-1)
    dir_edges[0] = dir_edges.pop(-1)
    dir_bins[0] = 0.

    var_bins = bins.tolist()
    var_bins.append(np.inf)

    if blowto:
        direction = direction + 180.
        direction[direction >= 360.] = direction[direction >= 360.] - 360

    table = histogram2d(x=var, y=direction, bins=[var_bins, dir_bins], normed=False)[0]
    # add the last value to the first to have the table of North winds
    table[:, 0] = table[:, 0] + table[:, -1]
    # and remove the last col
    table = table[:, :-1]
    if normed:
        table = table * 100 / table.sum()

    return dir_edges, var_bins, table

@docstring.copy(WindroseAxes.bar)
def wrbar(direction, var, ax=None, rmax=None, **kwargs):
    ax = WindroseAxes.from_ax(ax, rmax=rmax)
    ax.bar(direction, var, **kwargs)
    ax.set_legend()
    return ax

def wrscatter(direction, var, ax=None, rmax=None, *args, **kwargs):
    """
    Draw scatter plot
    """
    ax = WindroseAxes.from_ax(ax, rmax=rmax)
    direction = -np.array(direction) + np.radians(90)
    ax.scatter(direction, var, *args, **kwargs)
    return ax

# %%
