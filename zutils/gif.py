import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import os, imageio

import datetime

class gif_maker:
    """
    Simple tool for making gifs from matplotlib figures. 
    
    USAGE:
        (1) initialise the gif_maker giving it a name for the gif [without extension] and a directory [default is '.']
            
                e.g.:
                    gm = gif_maker('Eta with quivers', '../../folder_to_save_gif')
                    
        (2) run some loop that creates your frames and pass the matplotlib figure handle to the gif_maker each time. 
            Note you can update the same figure continually or you can 
            
                e.g.:
                    for i in np.arange(0, nloops):
                        fig = plt.figure(figsize=(10, 8)) # Create the figure
                        
                        # SOME USER-DEFINED FUNCTION # SOME USER-DEFINED FUNCTION
                        update_figure(fig, i) # SOME USER-DEFINED FUNCTION
                        # SOME USER-DEFINED FUNCTION # SOME USER-DEFINED FUNCTION
                        
                        gm.capture_fig(fig) # This is the gif_maker capture step. It saves a png image to file and tracks 
                                              it for splicing into the gif atr a later stage
                                              
        (3) make the gif once the loop is complete. frame rate can be spoecified at this point.
        
                e.g.:
                    
                                              
    """
    def __init__(self, gif_name, gif_dir='.'):

        self.gif_name = gif_name
        
        if not os.path.exists(gif_dir):
            os.mkdir(gif_dir)

        gif_stills_dir = os.path.join(gif_dir, 'stills')
        if not os.path.exists(gif_stills_dir):
            os.mkdir(gif_stills_dir)

        self.gif_dir = gif_dir
        self.gif_stills_dir = gif_stills_dir
        
        self.stills = []
        
    def capture_fig(self, fig):

        still_number = len(self.stills)
        savename = os.path.join(self.gif_stills_dir, '{}{}.png'.format(self.gif_name, still_number))
        
        fig.savefig(savename)
        self.stills += [savename]
        pass
    
    @property
    def gif_fullpath(self):
        
        return os.path.join(self.gif_dir, '{}.gif'.format(self.gif_name))
    
    def make_gif(self, fps=5):
        
        with imageio.get_writer(self.gif_fullpath, mode='I', fps=fps) as writer:
            for still in self.stills:
                image = imageio.imread(still)
                writer.append_data(image)
