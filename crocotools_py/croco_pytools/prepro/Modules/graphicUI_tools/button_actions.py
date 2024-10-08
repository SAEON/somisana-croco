import numpy as np
import matplotlib
# We want matplotlib to use a wxPython backend
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from threading import Thread
from time import sleep
import wx

from traits.api import *
from traitsui.api import View, Item, Group, HSplit, Handler, EnumEditor, FileEditor,DirectoryEditor
from traitsui.menu import NoButtons
from traitsui.wx.editor import Editor
from traitsui.basic_editor_factory import BasicEditorFactory

import matplotlib.pyplot as plt
import tools_make_grid

class ComputeGridThread(Thread):
    """
    This is the worker thread that 
    spawns the easygrid processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the easygrid computing loop
        """
        self.easy(self.inputs, self.outputs)
        if self.inputs.zview == 'topo':
            self.topo(self.outputs, self.topo_file,self.shp_file)
        if self.inputs.zview == 'mask':
            self.mask(self.outputs,self.shp_file)

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())
        easyparam = (self.inputs.nx,      self.inputs.ny,
                     self.inputs.tra_lon, self.inputs.tra_lat,
                     self.inputs.size_x,  self.inputs.size_y,
                     self.inputs.rot)

        self.display('----- min=dy %.2f,  max=dy %.2f\n' % (dy))
        self.display('----- min=dx %.2f,  max=dx %.2f'   % (dx))
        self.display(''.join(('--- computing grid (nx=%i, ny=%i, ',
                              'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                              'rot=%.2f)')) % easyparam)

        print('Done computing grid')
class ComputeSmthThread(Thread):
    """
    This is the worker thread that
    spawns the smoothing processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the smoothing computing loop
        """
        self.easy(self.inputs, self.outputs)
        self.topo(self.outputs, self.topo_file,self.shp_file,smooth=self.inputs_smth,
                  sgl_connect=self.single_connect)

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())

        easyparam = (self.inputs.nx,      self.inputs.ny,
                     self.inputs.tra_lon, self.inputs.tra_lat,
                     self.inputs.size_x,  self.inputs.size_y,
                     self.inputs.rot)

        self.display('----- min=dy %.2f,  max=dy %.2f\n' % (dy))
        self.display('----- min=dx %.2f,  max=dx %.2f'   % (dx))
        self.display(''.join(('--- computing grid (nx=%i, ny=%i, ',
                              'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                              'rot=%.2f)')) % easyparam)
        print('Done smoothing')

class ComputeZmThread(Thread):
    """
    This is the worker thread that 
    spawns the easygrid processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the easygrid computing loop
        """
 
        self.easy(self.inputs, self.outputs)

        self.topo(self.outputs, self.topo_file,self.shp_file,
                  smooth=self.inputs_smth,
                  sgl_connect=self.single_connect,
                  prt_grd=self.topo_prt)
        self.match_topo(self.topo_prt,self.outputs,self.openb)

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())
        easyparam = (self.inputs.nx,      self.inputs.ny,
                     self.inputs.tra_lon, self.inputs.tra_lat,
                     self.inputs.size_x,  self.inputs.size_y,
                     self.inputs.rot)

        self.display('----- min=dy %.2f,  max=dy %.2f\n' % (dy))
        self.display('----- min=dx %.2f,  max=dx %.2f'   % (dx))
        self.display(''.join(('--- computing grid (nx=%i, ny=%i, ',
                              'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                              'rot=%.2f)')) % easyparam)
        print('Done making zoom')

class ComputeC2cThread(Thread):
    """
    This is the worker thread that
    spawns the smoothing processing job
    """
    wants_abort = False

    def run(self):
        """
        Runs the smoothing computing loop
        """

        self.nest(self.topo_prt,self.inputs,self.outputs)
        self.topo(self.outputs, self.topo_file,self.shp_file,
                  smooth=self.inputs_smth,
                  hmin=np.nanmin(self.topo_prt.h),
                  hmax=np.nanmax(self.topo_prt.h),
                  sgl_connect=self.single_connect,
                  prt_grd=self.topo_prt,
                  coef=self.inputs.coef)
        self.match_topo(self.topo_prt,self.outputs,self.openb) 

        dx = (1 / self.outputs.pm.max(), 1 / self.outputs.pm.min())
        dy = (1 / self.outputs.pn.max(), 1 / self.outputs.pn.min())

        print('Done making AGRIF zoom')
class SaveGridThread(Thread):
    """
    
    """
    wants_abort = False

    def run(self):
        """
        Runs the easygrid save to nc loop
        """
        self.save2netcdf.create_grid_nc(self.outputs_file,self.inputs, self.outputs,prt_grd=self.prt_grd)

        if self.prt_grd is None or self.prt_grd[0]==False:
            easyparam = (self.inputs.nx, self.inputs.ny,
                         self.inputs.tra_lon, self.inputs.tra_lat,
                         self.inputs.size_x, self.inputs.size_y,
                         self.inputs.rot)
            self.display(''.join(('--- saving grid (nx=%i, ny=%i, ',
                                  'lon=%.2f, lat=%.2f, sx=%.2f, sy=%.2f, ',
                                  'rot=%.2f)\n')) % easyparam)
        else:
            easyparam = (self.outputs.h.shape[1]-2, self.outputs.h.shape[0]-2,
                         self.prt_grd[2], 
                         self.prt_grd[3],self.prt_grd[4],self.prt_grd[5],self.prt_grd[6])
            self.display(''.join(('--- saving grid (nx=%i, ny=%i,',
                                  'coef=%i, imin=%i, imax=%i, jmin=%i, ',
                                  'jmax=%i)\n')) % easyparam)


