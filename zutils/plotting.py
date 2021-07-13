import numpy as np
import matplotlib.pyplot as plt
    
class axis_layer():
    """
    A class to assist in laying out axes for publication. 

    ...

    Attributes
    ----------
    widths : int/float or 1xn list of int/float 
        The widths of each column of axes in cm. If a single numeric value is entered, only a single column will be created. 
        If an 1xn list of int/float is entered, n columns will be created. 
    
    heights : int/float or 1xm list of int/float 
        The heights of each row of axes in cm. If a single numeric value is entered, only a single row will be created. 
        If an 1xm list of int/float is entered, m rows will be created. 

    hspace: int/float or 1x(n-1) list of int/float | default = 1.
        Horizontal spacing in cm  between columns of axes. Specify a list 1x(n-1) for each column space, or a single int/float for equal spacing. 

    vspace: int/float or 1x(n-1) list of int/float | default = 1.
        Vertical spacing in cm between rows of axes. Specify a list 1x(m-1) for each row space, or a single int/float for equal spacing. 

    right/left/top/bottom: int/float | default = 1.
        right/left/top/bottom margin in cm.

    Methods
    -------
    lay(posx, posy, **kwargs)
        Creates an axis at the position posx, posy.
    
    """

    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        widths : int/float or 1xn list of int/float 
            The widths of each column of axes in cm. If a single numeric value is entered, only a single column will be created. 
            If an 1xn list of int/float is entered, n columns will be created. 
        
        heights : int/float or 1xm list of int/float 
            The heights of each row of axes in cm. If a single numeric value is entered, only a single row will be created. 
            If an 1xm list of int/float is entered, m rows will be created. 

        hspace: int/float or 1x(n-1) list of int/float | default = 1.
            Horizontal spacing in cm  between columns of axes. Specify a list 1x(n-1) for each column space, or a single int/float for equal spacing. 

        vspace: int/float or 1x(n-1) list of int/float | default = 1.
            Vertical spacing in cm between rows of axes. Specify a list 1x(m-1) for each row space, or a single int/float for equal spacing. 

        right/left/top/bottom: int/float | default = 1.
            right/left/top/bottom margin in cm.
        """

        self._left = 1.0
        self._right = 1.0
        self._top = 1.0
        self._bottom = 1.0
        
        self._hspace = [1.0]
        self._vspace = [1.0]
        
        self._widths = [5.0]
        self._heights = [5.0]
        
        self.update_class_with_loop(kwargs)

        self.verbose = True

        self._top_to_bottom = True # Never tested with this set to False
        
        fsx, fsy = self.get_figsize_cm()
        print('Figure size is {} x {} cm'.format(fsx, fsy))

    def lay(self, row, col, rowbleed=0, colbleed=0, force_new=False, frameon=False, resize=True, **kwargs):
        """
        Creates an axis at the position posx, posy.

        
        Parameters
        ----------
        row: int
            x position (column) of the axis
        col: int
            y position (row) of the axis

            Advanced:
            Rowbleed: int 
                This allows the axis to spread over multiple rows. To use this, specify the lowest row as the row input (above). Bleeding over outside the range will raise an exception.
            Colbleed: int 
                This allows the axis to spread over multiple columns. To use this, specify the leftmost column as the column input (above). Bleeding over outside the range will raise an exception.
        
        """
        if 'figure' in kwargs.keys():
            print('Figure specified')
            f = kwargs['figure']
            kwargs.pop("figure")
        else: 
            f = plt.gcf()
            
        fsx, fsy = self.get_figsize_inches()
        if resize:
            f.set_size_inches(fsx, fsy)
        else:
            print('Not resizing figure.')
        
        rect = self.get_pos_norm(col, row, rowbleed=rowbleed, colbleed=colbleed)
    
        if force_new:
            ax = f.add_axes(rect, frameon=frameon, **kwargs)
        else:
            ax = plt.axes(rect, **kwargs)
        # ax.set_title('row: {} | col: {}'.format(row, col))
        
        return ax
        
    def update_class_with_update(self, kwargs):
        
        self.__dict__.update(kwargs) # This won't work with object setters
        
    def update_class_with_loop(self, kwargs):
        
        for key in kwargs.keys():

            setattr(self, key, kwargs[key])
            
    def get_w(self):
        
        w = [self.left]
        for i in np.arange(0, len(self.widths)-1):
            w.append(self.widths[i])
            w.append(self.hspace[i])
        w.append(self.widths[-1])
        w.append(self.right)
        
        return w
        
    def get_h(self):
        
        h = [self.bottom]
        for i in np.arange(0, len(self.heights)-1):
            h.append(self.heights[i])
            h.append(self.vspace[i])
        h.append(self.heights[-1])
        h.append(self.top)
        
        if self.top_to_bottom:
            h.reverse()
            
        return h
        
    def get_figsize_cm(self):
        
        w = self.get_w()
        h = self.get_h()
        
        fsx = sum(w)
        fsy= sum(h)

        return fsx, fsy
    
    def get_figsize_inches(self):
        
        fsx, fsy = self.get_figsize_cm()
        
        fsx = fsx/2.54
        fsy = fsy/2.54
        
        return fsx, fsy
    
    def get_matrix(self):
    
        w = self.get_w()
        h = self.get_h()
        
        h = [self.top]
        for i in np.arange(0, len(self.heights)-1):
            h.append(self.heights[i])
            h.append(self.vspace[i])
        h.append(self.heights[-1])
        h.append(self.bottom)
        
        W, H = np.meshgrid(w, h)
        
        return W, H
        
    def get_pos(self, posx, posy, rowbleed=0, colbleed=0):
        
        W, H = self.get_matrix()
        fsx, fsy = self.get_figsize_cm()
        
        Ws = np.cumsum(W, axis=1)
        Hs = fsy-np.cumsum(H, axis=0)
            
        if not self.top_to_bottom:
            posy = len(self.heights)-posy-1
            
        pullx = posx*2+1
        pully = posy*2+1

        x = Ws[pully, pullx-1]
        y = Hs[pully, pullx-1]

        pullx_b = np.arange(posx*2+1, (posx+colbleed)*2+2)
        pully_b = np.arange((posy-rowbleed)*2+1, (posy)*2+2)
        
        if self.verbose:
            print('pullx_b {}'.format(pullx_b))
            print('pully_b {}'.format(pully_b))

        w = sum(W[pully, pullx_b])
        h = sum(H[pully_b, pullx])
                
        rect = [x, y, w, h]
        
        if self.verbose:
            
            print('colbleed {}'.format(colbleed))
            print('rowbleed {}'.format(rowbleed))

            print('posx {}'.format(posx))
            print('posy {}'.format(posy))
        
            print('pullx {}'.format(pullx))
            print('pully {}'.format(pully))
            
            print('pullx_b {}'.format(pullx_b))
            print('pully_b {}'.format(pully_b))
        
            print('W')
            print(W)
        
            print('H')
            print(H)
        
            print('Ws')
            print(Ws)
        
            print('Hs')
            print(Hs)
            
            print('rect')
            print(rect)
            
        return rect
    
    def get_pos_norm(self, posx, posy, rowbleed=0, colbleed=0):
    
        fsx, fsy = self.get_figsize_cm()
        
        rect = self.get_pos(posx, posy, rowbleed=rowbleed, colbleed=colbleed)
        
        rect_norm = [rect[0]/fsx, rect[1]/fsy, rect[2]/fsx, rect[3]/fsy]
        
        if self.verbose:
            
            print('rect_norm')
            print(rect_norm)
            
        return rect_norm
        
    @property
    def left(self):
        return self._left
    @property
    def right(self):
        return self._right
    @property
    def top(self):
        return self._top
    @property
    def bottom(self):
        return self._bottom
    
    @property
    def widths(self):
        return self._widths
    @property
    def heights(self):
        return self._heights

    @property
    def top_to_bottom(self):
        return self._top_to_bottom
    
    @property
    def hspace(self):
        v = self._hspace
        if type(v)==int or type(v)==float:
            v = [v]
            
        if type(v)==list:
            v = [float(e) for e in v]
            if len(v) == 1:
                v = v*(len(self.widths)-1)
            elif len(v) == len(self.widths)-1:
                pass
            else:
                raise(Exception('widths and hspace not compatible lengths'))
        else:
            raise(Exception('hspace must be a list, an int or a float'))
                
        return v
    
    @property
    def vspace(self):
        v = self._vspace
        if type(v)==int or type(v)==float:
            v = [v]
            
        if type(v)==list:
            v = [float(e) for e in v]
            if len(v) == 1:
                v = v*(len(self.heights)-1)
            elif len(v) == len(self.heights)-1:
                pass
            else:
                raise(Exception('heights and vspace not compatible lengths'))
        else:
            raise(Exception('vspace must be a list, an int or a float'))
                
        return v                
    
    @left.setter
    def left(self, v):
        if not (type(v)==int or type(v)==float):
            raise(Exception('must be int'))
        self._left = float(v)
    @right.setter
    def right(self, v):
        if not (type(v)==int or type(v)==float):
            raise(Exception('must be int'))
        self._right = float(v)
    @top.setter
    def top(self, v):
        if not (type(v)==int or type(v)==float):
            raise(Exception('must be int'))
        self._top = float(v)
    @bottom.setter
    def bottom(self, v):
        if not (type(v)==int or type(v)==float):
            raise(Exception('must be int'))
        self._bottom = float(v)
    @hspace.setter
    def hspace(self, v):
        if not (type(v)==int or type(v)==float or type(v)==list):
            raise(Exception('must be int or list'))
        self._hspace = v
    @vspace.setter
    def vspace(self, v):
        if not (type(v)==int or type(v)==float or type(v)==list):
            raise(Exception('must be int or list'))
        self._vspace = v
        
    @widths.setter
    def widths(self, v):
        if not (type(v)==int or type(v)==float  or type(v)==list):
            raise(Exception('widths must be int or list not {}'.format(type(v))))
            
        v = [float(e) for e in v]
        self._widths = v
    @heights.setter
    def heights(self, v):
        if not (type(v)==int or type(v)==float  or type(v)==list):
            raise(Exception('heights must be int or list not {}'.format(type(v))))
        v = [float(e) for e in v]
        self._heights = v

    @top_to_bottom.setter
    def top_to_bottom(self, v):
        if not type(v)==bool:
            raise(Exception('must be bool (logical)'))
            
        self._top_to_bottom = v

def hide_xticks(ax):
    """
    This just hides the x ticks of an axis. I'm surprised there isn't a simple matplotlib command for this.
        Inputs: 
                ax - single axis object for which the ticks must be written.  
    """
    for tick in ax.axes.get_xticklabels():
        tick.set_visible(False)

def hide_yticks(ax):
    """
    This just hides the y ticks of an axis. I'm surprised there isn't a simple matplotlib command for this.
        Inputs: 
                ax - single axis object for which the ticks must be written.  
    """
    for tick in ax.axes.get_yticklabels():
        tick.set_visible(False)