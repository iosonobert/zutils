# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 15:17:42 2018

@author: 20270917
"""
#%%
import numpy as np
import sys
import matplotlib.pyplot as plt 

from PyQt5.QtWidgets import (QWidget, QSlider, QHBoxLayout,
                             QLabel, QApplication, QPushButton,
                             QScrollBar)

from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.dates as mdates

from PyQt5.QtWidgets import (QWidget, QSlider, QHBoxLayout,
                             QLabel, QApplication, QPushButton,
                             QScrollBar, QAction)

from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

"""
Simple customisable netcdf viewer GUI for the pIMOS library. This is an attempted generalisation of the mrayson/mycurrents ocean mooring GUI. 
Hope is to generalise to:
    - Work with pIMOS/xrwrap objects
    - Plot any number of axes

Example Usage:
    params=[{'name':'U',
             'vmin': -0.4,
             'vmax': 0.4,
             'cmap': 'RdBu',
             'title': 'U vel [no QC applied]'},
            {'name':'U',
             'vmin': -0.4,
             'vmax': 0.4,
             'cmap': 'RdBu',
             'qc': True,
             'title': 'U vel [QC applied]'},
            {'name':'echo',
             'operation':'max',
             'operation_ax':2,
             'cmap': 'inferno',
             'vmin': 30,
             'vmax': 130,
             'title': 'Maximum Echo'},
           ]

    app = QtWidgets.QApplication([])
    foo = MyMainWindow(gui_data, params)
    foo.show()
    sys.exit(app.exec_())

"""

def on_xlims_change(event_ax, self):
    """
    Handler for axis changes by zooming. Disconnect handlers first. Then reconnect at the end. Or you'll loop forever. 
    """
    
    xl = event_ax.get_xlim()
    self.disconnect_lim_changed_callbacks()
    
    for i, ax in enumerate(self.window_widget.ax):
        if ax == event_ax:
            continue
        ax.set_xlim(xl)
        self.window_widget.canvas[i].draw()
        
    self.connect_lim_changed_callbacks()
    
    print('hi')
    print(event_ax)
    
class Window(QtWidgets.QWidget):
    def __init__(self, parent, n):
        super(Window, self).__init__(parent)

        vertical = False
        
        self.figure  = []
        self.canvas = []
        self.ax = []
        self.cax = []
        self.toolbar = []
        
        print('MAKING FIGURE AND AXIS')
        
        for i in np.arange(0, n):
            
            figure = Figure()
            gs = figure.add_gridspec(1, 20)
            
            canvas = FigureCanvas(figure)
            ax = figure.add_subplot(gs[0, :-1])
            cax = figure.add_subplot(gs[0, -1])
            toolbar = NavigationToolbar(canvas, self)
            if vertical:
                toolbar.setOrientation(QtCore.Qt.Vertical)
            toolbar.setIconSize(QtCore.QSize(16, 16))
            
#             mainLayout.addWidget(canvas,1,0)
            
            self.figure += [figure]
            self.canvas += [canvas]
            self.ax += [ax]
            self.cax += [cax]
            self.toolbar += [toolbar]
        
        mainLayout = QtWidgets.QGridLayout()
        if vertical:
            for i in np.arange(0, n):
                mainLayout.addWidget(self.toolbar[i],i, 0)
            for i in np.arange(0, n):
                mainLayout.addWidget(self.canvas[i],i, 1)
        else:
            for i in np.arange(0, n):
                mainLayout.addWidget(self.toolbar[i],i*2, 0)
            for i in np.arange(0, n):
                mainLayout.addWidget(self.canvas[i],i*2+1,0)
            
        self.setLayout(mainLayout)
        self.setWindowTitle("Flow Layout")

class MyMainWindow(QtWidgets.QMainWindow):

    start = 0
    length = 10
    n = 5000
    block_days=7
    
    def __init__(self, gui_data, params, parent=None):

        super(MyMainWindow, self).__init__(parent)
        
        n_axes = len(params)
        self.n_axes = n_axes
        self.gui_data = gui_data
        self.params = params
        
        self.ds_full = gui_data.ds.copy()  # We will trick the gui into thinking it has less data
        
        self.exitAction = QAction("&Exit", parent)
        self.exitAction.setStatusTip('Close applocation')
        self.exitAction.triggered.connect(self.close)

        self.menubar0 = self.menuBar()
        self.fileMenu = self.menubar0.addMenu('&File')
        # self.fileMenu.addAction(exitAct)
        self.fileMenu.addAction(self.exitAction)

        self.toolbar0 = self.addToolBar("toolbar")
        self.toolbar0.setStyleSheet("padding:10px")
        
        self.plotFullRecord_button = QPushButton('Plot Full [not recommended for large datasets]', self)
        self.plotFullRecord_button.setToolTip('Set plot limits to cover the full data')
        self.plotFullRecord_button.clicked.connect(self.plotFullRecord)
        self.plotFullRecord_button.setStyleSheet("padding:5px")
        self.toolbar0.addWidget(self.plotFullRecord_button)
        
        self.toolbar0.addSeparator()
        
        self.plotOneBlock_button = QPushButton('Plot First Week', self)
        self.plotOneBlock_button.setToolTip('Set plot limits to cover the first week')
        self.plotOneBlock_button.clicked.connect(self.plotOneBlock)
        self.plotOneBlock_button.setStyleSheet("padding:5px")
        self.toolbar0.addWidget(self.plotOneBlock_button)
        
        self.plotOneBlock_button = QPushButton('Back One Week', self)
        self.plotOneBlock_button.setToolTip('Set plot limits to cover the previous week')
        self.plotOneBlock_button.clicked.connect(self.backOneBlock)
        self.plotOneBlock_button.setStyleSheet("padding:5px")
        self.toolbar0.addWidget(self.plotOneBlock_button)
        
        self.plotOneBlock_button = QPushButton('Advance One Week', self)
        self.plotOneBlock_button.setToolTip('Set plot limits to cover the next week')
        self.plotOneBlock_button.clicked.connect(self.advanceOneBlock)
        self.plotOneBlock_button.setStyleSheet("padding:5px")
        self.toolbar0.addWidget(self.plotOneBlock_button)
        
        self.toolbar0.addSeparator()
        
        self.printLimits_button = QPushButton('Print Limits', self)
        self.printLimits_button.setToolTip('Print Limits of each axis.')
        self.printLimits_button.clicked.connect(self.printLimits)
        self.printLimits_button.setStyleSheet("padding:5px")
        self.toolbar0.addWidget(self.printLimits_button)
        
        

        self.toolbar1 = self.addToolBar("toolbar")

        self.start_slider = QSlider(QtCore.Qt.Horizontal, self)
        self.start_slider.setRange(0, 500)
        self.start_slider.setFocusPolicy(QtCore.Qt.NoFocus)
        self.start_slider.setPageStep(5)
        self.start_slider.valueChanged.connect(self.updateStart)
        self.start_slider.setTickPosition

        self.start_label = QLabel('0', self)
        self.start_label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.start_label.setMinimumWidth(200)

        self.length_slider = QSlider(QtCore.Qt.Horizontal, self)
        self.length_slider.setRange(self.block_days, self.block_days)
        self.length_slider.setValue(self.block_days)
        self.length_slider.setFocusPolicy(QtCore.Qt.NoFocus)
        self.length_slider.setPageStep(1)
        self.length_slider.valueChanged.connect(self.updateLength)
        self.length_slider.setSingleStep(1)

        self.length_label = QLabel('0', self)
        self.length_label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
        self.length_label.setMinimumWidth(200)

        self.toolbar1.addWidget(self.start_slider)
        self.toolbar1.addWidget(self.start_label)
        self.toolbar1.addWidget(self.length_slider)
        self.toolbar1.addWidget(self.length_label)

        QtWidgets.QMainWindow.insertToolBarBreak(self, self.toolbar1)
        
        if False:
            self.form_widget = FormWidget(self) 
            self.setCentralWidget(self.form_widget) 
        else:
            self.window_widget = Window(self, n_axes) 
            self.setCentralWidget(self.window_widget) 
            
        self.connect_lim_changed_callbacks()
        
    def close(self):

        print('Closing application')
        self.close() 

    def plotFullRecord(self):
        
        print('IM IN THE FULL RECORD PLOT BUTTON')
        
        xl = self.gui_data.ds.time
        xl = [min(xl), max(xl)]
        
        for i in np.arange(0, self.n_axes):
            self.window_widget.ax[i].set_xlim(xl)
            self.window_widget.canvas[i].draw()
            
        self.actually_plot()   
            
    def update_plot_lims(self):
        
        self.disconnect_lim_changed_callbacks()
        
        which_block = self.start_slider.value() 
        
        start = self.gui_data.ds.time.values[0]
        xl = [which_block*self.block_days, (which_block+1)*self.block_days]
        xl = [start + np.timedelta64(x, 'D') for x in xl]
            
        for i in np.arange(0, self.n_axes):
            print("UPDATING PLOT LIMS")
            self.window_widget.ax[i].set_xlim(xl)
            self.window_widget.canvas[i].draw()
    
        self.refresh_plot()
        
        self.connect_lim_changed_callbacks()
        
#     @staticmethod
        
    def connect_lim_changed_callbacks(self):
        """
        Only do the xlims for now, as we may have different plots with very different y lims.
        Would need some in deth handling of this. Not interested. 
        """
#         print('Connecting xlim_changed event handler')
        for ax in self.window_widget.ax:
            
            ax.callbacks.connect('xlim_changed', lambda x: on_xlims_change(x, self))
#                 ax.callbacks.connect('ylim_changed', lambda x: on_ylims_change(x, self))

                    
    def disconnect_lim_changed_callbacks(self):
        print('Disonnecting xlim_changed event handler')
        for ax in self.window_widget.ax:
#             
            if 'xlim_changed' in ax.callbacks.callbacks.keys():
                eids = list(ax.callbacks.callbacks['xlim_changed'].keys())

                for eid in eids:
                    ax.callbacks.disconnect(eid)
            if 'ylim_changed' in ax.callbacks.callbacks.keys():
                eids = list(ax.callbacks.callbacks['ylim_changed'].keys())

                for eid in eids:
                    ax.callbacks.disconnect(eid)
    
    @property
    def n_blocks(self):
        print(self.block_days)
        
        start = self.gui_data.ds.time.values[0]
        end = self.gui_data.ds.time.values[-1]

        start = self.gui_data.ds.time.values[0]
        end = self.gui_data.ds.time.values[-1]
        
        length = end-start
        length_s = length / np.timedelta64(1, 's')

        length_week = np.int(np.ceil(length_s/(self.block_days*86400)))

        return length_week-1

    def plotOneBlock(self):
                
        print(self.block_days)
        
        start = self.gui_data.ds.time.values[0]
        end = self.gui_data.ds.time.values[-1]
        
        length = end-start
        length_s = length / np.timedelta64(1, 's')

        length_week = np.int(np.ceil(length_s/(self.block_days*86400)))
        length_week
        
        self.start_slider.setRange(0, length_week-1)
        self.start_slider.setValue(0)
        
        self.update_plot_lims()
        
    def advanceOneBlock(self):
        
        self.moveOneBlock(1)
    
    def backOneBlock(self):
        
        self.moveOneBlock(-1)
        
    def moveOneBlock(self, step):
        
        which_block = self.start_slider.value() 
        try:
            which_block += 1*step
            self.start_slider.setValue(which_block)
            self.update_plot_lims()
        except:
            print("Failed to move block")
            
    def printLimits(self):
        
        print("")
        print("")

        for i in np.arange(0, self.n_axes):
            xl = self.window_widget.ax[i].get_xlim() 
            yl = self.window_widget.ax[i].get_ylim() 
            
            print('AXIS {}'.format(i))
            print('     X LIM {} to {}'.format(mdates.num2date(xl[0]), mdates.num2date(xl[1])))
            print('     Y LIM {}'.format(yl))
            
    def updateStart(self, value):

        self.start = value
        self.start_label.setText('Week: ' + str(value) + ' of ' + str(self.n_blocks))
#         self.update_plot()
        self.update_plot_lims()
        
    def updateLength(self, value):

        self.length = value
        self.length_label.setText('Length: ' + str(value))
        
    def prep_dataarray_vals(self, i):
    
        prep_dict = self.params[i]
        
        if 'qc' in prep_dict.keys():
            qc = prep_dict['qc']
        else:
            qc = False
            
        if qc:
            da = self.gui_data.get_qaqc_var(prep_dict['name'])
        else:
            da = self.gui_data.ds[prep_dict['name']]
            
        val = da.values
        
        if 'operation' in prep_dict.keys() and 'operation_ax' in prep_dict.keys():
            operation = prep_dict['operation']
            operation_ax = prep_dict['operation_ax']

            if operation.lower() == 'max':
                val = np.nanmax(val, axis=operation_ax)
            if operation.lower() == 'mean':  
                val = np.nanmean(val, axis=operation_ax)
            if operation.lower() == 'min':  
                val = np.nanmin(val, axis=operation_ax)

        return da, val

    def refresh_plot(self):
        print('Refreshing plot')
        self.disconnect_lim_changed_callbacks()
        
        self.gui_data.ds = self.ds_full.copy()
    
        xl = self.window_widget.ax[0].get_xlim()
        start_date = mdates.num2date(xl[0]).replace(tzinfo=None)
        end_date = mdates.num2date(xl[1]).replace(tzinfo=None)

        print(start_date)
        print(end_date)

#         print(len(self.ds_full.time))
#         print(len(self.gui_data.ds.time))
        self.gui_data.ds = self.gui_data.ds.sel(time=slice(start_date, end_date))
#         print(len(self.gui_data.ds.time))
        
#         for i in np.arange(0, self.n_axes):
#             print(xl)
#             self.window_widget.ax[i].cla()  # Clear the canvas.
        self.actually_plot()
            
        self.connect_lim_changed_callbacks()
        
    def actually_plot(self):
        
        self.disconnect_lim_changed_callbacks()
        gui_data = self.gui_data

        for i in np.arange(0, self.n_axes):
            self.window_widget.ax[i].cla()  # Clear the canvas.
             
            da, vals = self.prep_dataarray_vals(i)
            
            # cmap handling
            try:
                cmap = self.params[i]['cmap']
            except:
                cmap = plt.get_cmap('cool')
            
            coordlist = [key for key in da.coords.keys()]
            if coordlist == ['time']:
                self.window_widget.ax[i].plot(gui_data.ds.time, vals)
                self.window_widget.ax[i].set_xlim([gui_data.ds.time[0], gui_data.ds.time[-1]])
                
                self.window_widget.cax[i].get_xaxis().set_visible(False)
                self.window_widget.cax[i].get_yaxis().set_visible(False)

            elif 'time' in coordlist and  'distance'  in coordlist:
                # vmin/vmax handling
                try:
                    mb = self.window_widget.ax[i].pcolor(gui_data.ds.time, gui_data.ds.distance, vals, cmap=cmap,
                                                   vmin=self.params[i]['vmin'],
                                                   vmax=self.params[i]['vmax'])
                    plt.colorbar(mb, cax=self.window_widget.cax[i])
                except:
                    self.window_widget.ax[i].pcolor(gui_data.ds.time, gui_data.ds.distance, vals, cmap=cmap,
                                                   cax=self.window_widget.cax[i])
                    plt.colorbar(mb, cax=self.window_widget.cax[i])

                try:
                    self.window_widget.ax[i].plot(gui_data.ds.time, gui_data.ds.depth, 'k')
                except:
                    pass
            elif 'time' in coordlist and  'freqs'  in coordlist:
                # vmin/vmax handling
                try:
                    mb = self.window_widget.ax[i].pcolor(gui_data.ds.time, gui_data.ds.freqs, vals.T, cmap=cmap,
                                                   vmin=self.params[i]['vmin'],
                                                   vmax=self.params[i]['vmax'])
                    plt.colorbar(mb, cax=self.window_widget.cax[i])
                except:
                    self.window_widget.ax[i].pcolor(gui_data.ds.time, gui_data.ds.freqs, vals.T, cmap=cmap,
                                                   cax=self.window_widget.cax[i])
                    plt.colorbar(mb, cax=self.window_widget.cax[i])

            else:
                print(coordlist)
                error
                
            try:
                self.window_widget.ax[i].set_title(self.params[i]['title'])
            except:
                pas
                
            self.window_widget.ax[i].grid()
            self.window_widget.canvas[i].draw()
    
        self.gui_data.ds = self.ds_full.copy()

        self.connect_lim_changed_callbacks()
        
def main(gui_data, params):
    
    app = QtWidgets.QApplication([])
    foo = MyMainWindow(gui_data, params)
    foo.show()
    sys.exit(app.exec_())
    
# %%
