import math, os
import numpy as np

import cartopy.crs as ccrs
# from pyproj import Proj
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from zutils import plotting as zplot

# returns square of distance b/w two points 
def length_square(X, Y, verbose=True): 
    xDiff = X[0] - Y[0] 
    yDiff = X[1] - Y[1] 

    return xDiff * xDiff + yDiff * yDiff
    
def get_angles(A, B, C, verbose=True): 
      
    # Square of lengths be a2, b2, c2 
    a2 = length_square(B, C) 
    b2 = length_square(A, C) 
    c2 = length_square(A, B) 
  
    # length of sides be a, b, c 
    a = math.sqrt(a2); 
    b = math.sqrt(b2); 
    c = math.sqrt(c2); 
  
    # From Cosine law 
    alpha = math.acos((b2 + c2 - a2) /
                         (2 * b * c)); 
    betta = math.acos((a2 + c2 - b2) / 
                         (2 * a * c)); 
    gamma = math.acos((a2 + b2 - c2) / 
                         (2 * a * b)); 
  
    # Converting to degree 
    alpha = alpha * 180 / math.pi; 
    betta = betta * 180 / math.pi; 
    gamma = gamma * 180 / math.pi; 
  
    # printing all the angles 
    if verbose:
        print("alpha : %f" %(alpha)) 
        print("betta : %f" %(betta))
        print("gamma : %f" %(gamma))

    return [alpha, betta, gamma], [a2, b2, c2]

class Tgrid:
    """
    Class for dealing with triangular grids. 
    
    Borrows heaviliy from mrayson's ugrid (part of sfoda), simplified for triangular only grids.  
    """
    
    def __init__(self, node_xyz, connectivity):
        """
        inputs:
            - node_xyz: n-point list of 3 point lists. One inner list for each node in an n-node mesh, and 
                            each inner list is the [X, Y, Z] of that node.   
            - connectivity: an m-point list where of 3 point lists. One inner list for each eklement in an m-element mesh, 
                            and each inner list are the indixes of the 3 nodes in that triangular element. 
        """
        
        n = [len(con) for con in connectivity]
        
        self.node_xyz = node_xyz
        self.connectivity = connectivity
        self.n = n
        
        if not self.is_triangular:
            raise(Exception('Currently supports triangular elements only'))
        else:
            print('This is a triangular mesh. ')
        
        self.calculate_element_vertices()
        self.calculate_element_xyz()
        self.calculate_t_geom_stats()
        
    @property
    def is_triangular(self):
        
        return np.all(np.array(self.n) == 3)
        
    def calculate_element_vertices(self):
        print('Prepping verts.')
        
        self.element_vertices = [self.get_element_vertex(e) for e in self.connectivity]
        print('Prepped.')

    def calculate_element_xyz(self):
        print('Prepping xyzs.')

        self.element_xyx = [self.get_element_xyz(e) for e in self.connectivity]

        self.element_x = [e[0] for e in self.element_xyx]
        self.element_y = [e[1] for e in self.element_xyx]
        self.element_z = [e[2] for e in self.element_xyx]
        print('Prepped.')

    def get_element_vertex(self, e):
        return [[self.node_xyz[e_, 0], self.node_xyz[e_, 1]] for e_ in e]

    def get_element_xyz(self, e):
        """
        Calculate the z of an element as the mean of the vertices. 
        """
        return np.mean([self.node_xyz[e_, 0] for e_ in e]), np.mean([self.node_xyz[e_, 1] for e_ in e]), np.mean([self.node_xyz[e_, 2] for e_ in e])
    
    def calculate_t_geom_stats(self):
        """
        Calculate largest and smallest angle, as well as largest and smallest sides
        """
        a_large = []
        l_large = []
        a_small = []
        l_small = []

        print('Calculating mesh geometry stats by very inefficient loop')
        for ev in self.element_vertices:
            a, l = get_angles(ev[0], ev[1], ev[2], verbose=False)
            a_large += [max(a)]
            l_large += [max(l)]

            a_small += [min(a)]
            l_small += [min(l)]

        print('Inefficient loop complete [not slow enough to bother making more efficient]')

        self.theta_min = a_small
        self.theta_max = a_large
        self.l_min = l_small
        self.l_max = l_large

    
    def to_shape(self, outfile):
        import shapefile
        from shapely.geometry import Polygon


        with shapefile.Writer(outfile) as w:
            
            w.field('Element', 'C')
            w.field('element_z', 'N', decimal=1)
            w.field('theta_max', 'N', decimal=0)
            w.field('theta_min', 'N', decimal=0)
            w.field('l_max', 'N', decimal=1)
            w.field('l_min', 'N', decimal=1)

            for ii in np.arange(0, len(self.element_z)):
                w.poly([self.element_vertices[ii]])
                
                w.record(ii, 
                        self.element_z[ii], 
                        self.theta_max[ii], 
                        self.theta_min[ii], 
                        self.l_max[ii], 
                        self.l_min[ii])

            w.close()

    def quick_figplot(self, colorby='z', clim = None, cmap = 'gist_earth', edgecolors='k', colorbar_label=''):
        
        def get_zl(xl, yl):
            """
            Return an axis layer object scaled to the aspect ratio of the mesh.
            """
            x = 8
            y = x/np.diff(xl)[0]*np.diff(yl)[0]
            
            print(y)
            
            widths=[x, x, 0.5]     # Width of each column of axes
            heights=[y, y]    # Height of each row of axes
                            # Use list to set each spacing seperately
            vspace = hspace = 0.5           # Vertical spacing between rows. Use single value to set all vertical spaces evenly

            bottom = 0.5
            top = 0.5

            zl = zplot.axis_layer(widths=widths, heights=heights, hspace=hspace, vspace=vspace, 
                                bottom=bottom, top=top, right=3, left=2)
            zl.verbose = False # Reduce printouts

            return zl

        xl = [min(self.element_x), max(self.element_x)]
        yl = [min(self.element_y), max(self.element_y)]

        if type(colorby) == str:
            if colorby.lower() == 'z':
                colour_vector = self.element_z
                colorbar_label = "Still water level [m MSL]"
            elif colorby.lower() == 'theta_min':
                colour_vector = self.theta_min
                colorbar_label = "Min angle"
            elif colorby.lower() == 'theta_max':
                colour_vector = self.theta_max
                colorbar_label = "Max angle"
            elif colorby.lower() == 'l_max':
                colour_vector = self.l_max
                colorbar_label = "Max side length"
            elif colorby.lower() == 'l_min':
                colour_vector = self.l_max
                colorbar_label = "Min side length"
            elif colorby.lower() == 'log10_l_max':
                colour_vector = np.log10(self.l_max)
                colorbar_label = "log10(Max side length)"
            elif colorby.lower() == 'log10_l_min':
                colour_vector = np.log10(self.l_max)
                colorbar_label = "log10(Min side length)"
            else:
                raise(Exception('Unrecognised colorby string'))
        else:
            assert(type(colorby) == list, 'Colorby must be a string or a list')
            assert(len(colorby) == len(self.element_z), 'Colorby list not the right length')
            colour_vector = colorby

        zl = get_zl(xl, yl)

        lower = np.min(colour_vector)
        upper = np.max(colour_vector)

    #     cm = 'seismic_r'
        
        cm = plt.cm.get_cmap(cmap)
        colors = cm((colour_vector-lower)/(upper-lower))

        poly = PolyCollection(self.element_vertices, facecolors = colors, edgecolors=edgecolors, linewidths=(0.05,), cmap=cm)
        poly.set_alpha(1)

        fig = plt.figure(dpi=300)
        ax = zl.lay(1, 0, rowbleed=1, colbleed=1, projection=ccrs.PlateCarree())

        if True:
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                        linewidth=0.5, color='gray', alpha=1, linestyle='--')

            gl.xlabels_bottom = False
            gl.ylabels_right = False

        plt.xlim(xl)
        plt.ylim(yl)

        ax.add_collection(poly)

        poly.set(array=np.array(colour_vector))
        poly.set_clim(clim)

        cax = zl.lay(1, 2)
        fig.colorbar(poly, label=colorbar_label, cax=cax)

        return fig



