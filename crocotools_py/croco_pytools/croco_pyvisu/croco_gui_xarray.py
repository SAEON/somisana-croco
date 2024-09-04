#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#

import sys
import os
from time import sleep
from shutil import copyfile
from threading import Thread
import wx
import xarray as xr
import numpy as np
from numpy.matlib import repmat
import numpy.ma as ma
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
from CrocoXarray import Croco
from derived_variables import get_pv, get_zetak, get_dtdz, get_richardson
from myplot import plotCurv, mypcolor
from matplotlib.animation import FFMpegWriter
from observer import Publisher, Subscriber

wildcard = "Netcdf Files (*.nc)|*.nc"
figsize = [6, 5]


########################################################################

class SectionFrame(wx.Frame):
    """
    Window class to plot longitude and latitude sections.
    The window contains a canvas to plot, several buttons (animation, zoom in, zoom out,
    print and reset color) and text control (start time and end time for animation,
    min and max color for the colorbar)

    Attributes:
    croco : croco instance to study
    variableName : Name of the variable to plot in the canvas
    variable : 3D dataarray (x,y,t) of the variable to plot
    x : numpy array of x coordinates
    y : numpy array of y coordinates
    typSection: type of slice XZ or YZ
    sliceCoord: coordinate of the slice (latitude for XZ, latitude for YZ)
    timeIndex: current time index to plot
    subscriber: instance of class to which send notification on click event
    """

    def __init__(self, croco=None, variableName=None, variable=None,
                 x=None, y=None, typSection=None, sliceCoord=None,
                 sliceIndex=None, timeIndex=None, subscriber=None):
        """ return a SectioFrame instance """

        # Initialize the variables of the class
        self.croco = croco
        self.variable = variable
        self.x = x
        self.y = y
        self.variableName = variableName
        self.typSection = typSection
        self.sliceCoord = sliceCoord
        self.sliceIndex = sliceIndex
        self.timeIndex = timeIndex

        # Create the window
        wx.Frame.__init__(self, None, wx.ID_ANY, title='Section')

        # Now create the Panel to put the other controls on.
        self.panel = wx.Panel(self, wx.ID_ANY)

        # and a few controls
        self.figure = Figure()
        self.axes = self.figure.add_axes([0, 0, 1, 1])
        self.canvas = FigureCanvas(self.panel, -1, self.figure)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Hide()

        self.TimeLabel = wx.StaticText(self.panel, -1, label="Time", style=wx.ALIGN_CENTER)
        self.TimeMinusBtn = wx.Button(self.panel, wx.ID_ANY, "<")
        self.TimeTxt = wx.TextCtrl(self.panel, wx.ID_ANY, " ", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.TimePlusBtn = wx.Button(self.panel, wx.ID_ANY, ">")
        self.ZoomBtn = wx.Button(self.panel, wx.ID_ANY, "Zoom")
        self.PanBtn = wx.Button(self.panel, wx.ID_ANY, "Pan")
        self.HomeBtn = wx.Button(self.panel, wx.ID_ANY, "Home")
        self.SavePlotBtn = wx.Button(self.panel, wx.ID_ANY, "Save Plot")

        self.AnimationBtn = wx.Button(self.panel, wx.ID_ANY, "Animation")
        self.startTimeTxt = wx.TextCtrl(self.panel, wx.ID_ANY, "1", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.endTimeTxt = wx.TextCtrl(self.panel, wx.ID_ANY, "1", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.AbortBtn = wx.Button(self.panel, wx.ID_ANY, "Abort")
        self.SaveAnimBtn = wx.Button(self.panel, wx.ID_ANY, "Save Anim")

        self.ResetColorBtn = wx.Button(self.panel, wx.ID_ANY, "Reset Color")
        self.MinColorTxt = wx.TextCtrl(self.panel, wx.ID_ANY, "Min Color", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.MaxColorTxt = wx.TextCtrl(self.panel, wx.ID_ANY, "Max Color", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)

        if self.typSection == "XY":
            self.TopoBtn = wx.ToggleButton(self.panel, wx.ID_ANY, "Topo", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
            self.TopoTxt = wx.TextCtrl(self.panel, wx.ID_ANY, "10", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)

        # bind the menu event to an event handler
        self.canvas.mpl_connect('button_press_event', self.onFigureClick)
        self.AnimationBtn.Bind(wx.EVT_BUTTON, self.onAnimationBtn)
        self.AbortBtn.Bind(wx.EVT_BUTTON, self.onAbortBtn)
        self.startTimeTxt.Bind(wx.EVT_TEXT_ENTER, self.onstartTimeTxt)
        self.endTimeTxt.Bind(wx.EVT_TEXT_ENTER, self.onendTimeTxt)
        self.ZoomBtn.Bind(wx.EVT_BUTTON, self.onZoomBtn)
        self.PanBtn.Bind(wx.EVT_BUTTON, self.onPanBtn)
        self.HomeBtn.Bind(wx.EVT_BUTTON, self.onHomeBtn)
        self.SavePlotBtn.Bind(wx.EVT_BUTTON, self.onSavePlotBtn)
        self.ResetColorBtn.Bind(wx.EVT_BUTTON, self.onResetColorBtn)
        self.MinColorTxt.Bind(wx.EVT_TEXT_ENTER, self.onMinColorTxt)
        self.MaxColorTxt.Bind(wx.EVT_TEXT_ENTER, self.onMaxColorTxt)
        self.AbortBtn.Bind(wx.EVT_BUTTON, self.onAbortBtn)
        self.SaveAnimBtn.Bind(wx.EVT_BUTTON, self.onSaveAnimBtn)
        self.TimeMinusBtn.Bind(wx.EVT_BUTTON, self.onTimeMinusBtn)
        self.TimeTxt.Bind(wx.EVT_TEXT_ENTER, self.onTimeTxt)
        self.TimePlusBtn.Bind(wx.EVT_BUTTON, self.onTimePlusBtn)

        if self.typSection == "XY":
            self.TopoBtn.Bind(wx.EVT_TOGGLEBUTTON, self.onTopoBtn)
            self.TopoTxt.Bind(wx.EVT_TEXT_ENTER, self.onTopoTxt)

        self.showPosition = self.CreateStatusBar(2)
        self.showPosition.SetStatusText("x=   , y=  ", 1)
        self.showPosition.SetStatusWidths([-1, 150])

        self.__do_layout()

        # Register main window subscriber to send new location level/latitude/longitude 
        # on window click
        self.sub = subscriber
        self.pub = Publisher()
        self.pub.register(self.sub, self.sub.update)

        # Set  variables of the class
        if croco is not None:
            self.time = self.croco.wrapper._get_date(self.timeIndex)
            self.TimeTxt.SetValue(str(self.time))
        if typSection == "XY":
            self.xlabel = "Longitude"
            self.ylabel = "Latitude"
        elif typSection == "XZ":
            self.xlabel = "Longitude"
            self.ylabel = "Depth"
            self.slice = "Latitude"
        elif typSection == "YZ":
            self.xlabel = "Latitude"
            self.ylabel = "Depth"
            self.slice = "Longitude"
        if croco is not None:
            timeMin = self.croco.wrapper._get_date(0)
            timeMax = self.croco.wrapper._get_date(self.croco.wrapper.ntimes - 1)
            self.startTimeTxt.SetValue(str(timeMin))
            self.startTime = timeMin
            self.startTimeIndex = 0
            self.endTimeTxt.SetValue(str(timeMax))
            self.endTime = timeMax
            self.endTimeIndex = self.croco.wrapper.ntimes - 1

    def __do_layout(self):
        """
        Use a sizer to layout the controls, stacked vertically or horizontally
        """
        topSizer = wx.BoxSizer(wx.VERTICAL)
        canvasSizer = wx.BoxSizer(wx.VERTICAL)
        timeSizer = wx.BoxSizer(wx.HORIZONTAL)
        buttonsSizer = wx.BoxSizer(wx.HORIZONTAL)
        colorSizer = wx.BoxSizer(wx.HORIZONTAL)

        canvasSizer.Add(self.canvas, 0, wx.ALL, 5)

        timeSizer.Add(self.TimeLabel, 0, wx.ALL, 5)
        timeSizer.Add(self.TimeMinusBtn, 0, wx.ALL, 5)
        timeSizer.Add(self.TimeTxt, 0, wx.ALL, 5)
        timeSizer.Add(self.TimePlusBtn, 0, wx.ALL, 5)
        timeSizer.Add(self.ZoomBtn, 0, wx.ALL, 5)
        timeSizer.Add(self.PanBtn, 0, wx.ALL, 5)
        timeSizer.Add(self.HomeBtn, 0, wx.ALL, 5)
        timeSizer.Add(self.SavePlotBtn, 0, wx.ALL, 5)

        buttonsSizer.Add(self.AnimationBtn, 0, wx.ALL, 5)
        buttonsSizer.Add(self.startTimeTxt, 0, wx.ALL, 5)
        buttonsSizer.Add(self.endTimeTxt, 0, wx.ALL, 5)
        buttonsSizer.Add(self.AbortBtn, 0, wx.ALL, 5)
        buttonsSizer.Add(self.SaveAnimBtn, 0, wx.ALL, 5)

        colorSizer.Add(self.ResetColorBtn, 0, wx.ALL, 5)
        colorSizer.Add(self.MinColorTxt, 0, wx.ALL, 5)
        colorSizer.Add(self.MaxColorTxt, 0, wx.ALL, 5)

        if self.typSection == "XY":
            colorSizer.Add(self.TopoBtn, 0, wx.ALL, 5)
            colorSizer.Add(self.TopoTxt, 0, wx.ALL, 5)

        topSizer.Add(canvasSizer, 0, wx.CENTER)
        topSizer.Add(timeSizer, 0, wx.ALL | wx.EXPAND, 5)
        topSizer.Add(buttonsSizer, 0, wx.ALL | wx.EXPAND, 5)
        topSizer.Add(colorSizer, 0, wx.ALL | wx.EXPAND, 5)

        self.panel.SetSizer(topSizer)
        topSizer.Fit(self)

        self.Layout()
        
        if self.typSection == "XY":
            self.toggletopo = self.TopoBtn.GetValue()
            self.nbtopo = int(self.TopoTxt.GetValue())
        else:
            self.toggletopo = None
            self.nbtopo = None

    # Event handler

    # Event handler on plot canvas
    def onFigureClick(self, event):
        """
        Event handler for the button click on plot
        send new location to main window
        """
        # Retreive new location on canvas click
        self.xPress, self.yPress = event.xdata, event.ydata
        level = longitude = latitude = None
        if self.typSection == "XY":
            longitude = self.xPress
            latitude = self.yPress
        elif self.typSection == "XZ":
            longitude = self.xPress
            level = self.yPress
        elif self.typSection == "YZ":
            latitude = self.xPress
            level = self.yPress
        # send new location to main window
        self.pub.dispatch(level, longitude, latitude)

    def ShowPosition(self, event):
        """Show the current cursor position at the bottom right of the window """
        if event.inaxes:
            self.showPosition.SetStatusText(
                "x={:5.1f}  y={:5.1f}".format(event.xdata, event.ydata), 1)

    def onAnimationBtn(self, event):
        """
        Launch the animation in a separate thread to keep the GUI alive
        Prepare tools to save the movie
        """
        # Initialize creation of movie
        self.movie = FFMpegWriter()
        self.movie.setup(self.figure, 'moviez.mp4', dpi=100)

        self.shouldAbort = False

        # Launch animate in a thread
        self.thread = Thread(target=self.animate)
        self.thread.start()

    def onAbortBtn(self, event):
        self.shouldAbort = True

    def animate(self):
        """
        Calculate and draw the variable from start Time to end Time
        """
        for i in range(self.startTimeIndex, self.endTimeIndex + 1):
            sleep(0.5)
            self.timeIndex = i
            # Get the variable at the new timeIndex
            self.updateVariableZ(setlim=False)
            # Draw the variable in the main thread
            wx.CallAfter(self.drawz, setlim=False, setcol=False, anim=True)
            if self.shouldAbort:
                break
        # Go to animateDone after aborting or ending the animation
        wx.CallAfter(self.animateDone, self.shouldAbort)

    def animateDone(self, abort):
        """
        method called at the end of the animation or after aborting it
        """
        # Close the movie
        self.movie.finish()

        # Remove the movie if aborting
        if self.shouldAbort:
            os.system('rm -rf ' + 'moviez.mp4')
        self.shouldAbort = False
        self.time = self.croco.wrapper._get_date(self.timeIndex)
        self.TimeTxt.SetValue(str(self.time))

    # Event handler for Save animation
    def onSaveAnimBtn(self, event):
        """
        Copy file movie*.mp4 in the right place with the right name
        """
        if os.path.isfile('moviez.mp4'):
            printDir = self.croco.startDir + "/Figures_" + self.croco.get_run_name() + "/"
            if not os.path.isdir(printDir):
                os.mkdir(printDir)
            time1 = str(self.croco.wrapper._get_date(self.startTimeIndex))
            time2 = str(self.croco.wrapper._get_date(self.endTimeIndex))
            filename = "{:s}_{:s}{:4.1f}_Time{:s}-{:s}.mp4".format(self.variableName, self.slice, self.sliceCoord,
                                                                   time1, time2).replace(" ", "")
            copyfile('moviez.mp4', printDir + filename)
            os.system('rm -rf ' + 'moviez.mp4')

    def onstartTimeTxt(self, event):
        """Event handler for Enter key in start time text """
        self.startTime = float(self.startTimeTxt.GetValue())
        times = self.croco.get_coord(self.variableName, direction='t')
        # find nearest index corresponding to instant time to plot
        self.startTimeIndex = min(range(len(times)), key=lambda j: abs(self.startTime - times[j]))
        self.startTime = self.croco.wrapper._get_date(self.startTimeIndex)
        self.startTimeTxt.SetValue(str(self.startTime))

    # Event handler for Time dialog
    def onendTimeTxt(self, event):
        """Event handler for Enter key in end time text """
        self.endTime = float(self.endTimeTxt.GetValue())
        times = self.croco.get_coord(self.variableName, direction='t')
        # find nearest index corresponding to instant time to plot
        self.endTimeIndex = min(range(len(times)), key=lambda j: abs(self.endTime - times[j]))
        self.endTime = self.croco.wrapper._get_date(self.endTimeIndex)
        self.endTimeTxt.SetValue(str(self.endTime))

    def onTimeTxt(self, event):
        """Event handler for Enter key in end time text """
        time = float(self.TimeTxt.GetValue())
        times = self.croco.get_coord(self.variableName, direction='t')
        # find index corresponding to the nearest instant time to plot
        self.timeIndex = min(range(len(times)), key=lambda j: abs(time - times[j]))
        self.time = self.croco.wrapper._get_date(self.timeIndex)
        self.TimeTxt.SetValue(str(self.time))
        self.updateVariableZ(setlim=False)
        self.drawz(setlim=False)

    def onTimeMinusBtn(self, event):
        self.timeIndex = max(self.timeIndex - 1, 0)
        self.time = self.croco.wrapper._get_date(self.timeIndex)
        self.TimeTxt.SetValue(str(self.time))
        self.updateVariableZ()
        self.drawz(setlim=False)

    def onTimePlusBtn(self, event):
        self.timeIndex = min(self.timeIndex + 1, self.croco.wrapper.ntimes - 1)
        self.time = self.croco.wrapper._get_date(self.timeIndex)
        self.TimeTxt.SetValue(str(self.time))
        self.updateVariableZ()
        self.drawz(setlim=False)

    def notify(self, ax):
        self.xlim = ax.get_xlim()
        self.ylim = ax.get_ylim()
        # print("notify:",ax.get_xlim(),ax.get_ylim())

    # Event handler for zoom
    def onZoomBtn(self, event):
        """
        Event handler for the button click Zoom in button
        """
        # self.figure.RS.set_active(True)
        self.toolbar.zoom()

    # Event handler for zoom
    def onPanBtn(self, event):
        """
        Event handler for the button click Zoom in button
        """
        self.toolbar.pan()

    def onHomeBtn(self, event):
        """
        Event handler for the button click Zoom out button
        """
        self.drawz(setlim=True, setcol=False)
        # self.toolbar.home()

    # Event handler for Print
    def onSavePlotBtn(self, event):
        """
        Event handler for the button click Print button
        """
        printDir = self.croco.startDir + "/Figures_" + self.croco.get_run_name() + "/"
        if not os.path.isdir(printDir):
                os.mkdir(printDir)
        filename = self.title.replace(',', '_').replace(" ", "") + ".png"
        os.system('rm -rf ' + printDir + filename)
        self.figure.savefig(printDir + filename, dpi=self.figure.dpi)
        # self.toolbar.save_figure(None)

    # Event handler for Color setup
    def onResetColorBtn(self, event):
        """
        Event handler for the button click Reset Color button
        """
        self.drawz(setlim=False, setcol=True)

    def onMinColorTxt(self, event):
        """Event handler for Enter key in Min Color text """
        self.clim[0] = float(self.MinColorTxt.GetValue())
        self.drawz(setlim=False, setcol=False)

    def onMaxColorTxt(self, event):
        """Event handler for Enter key in Max Color text """
        self.clim[1] = float(self.MaxColorTxt.GetValue())
        self.drawz(setlim=False, setcol=False)

    # Event handler for topography toggle button
    def onTopoBtn(self, evnet):
        self.toggletopo = self.TopoBtn.GetValue()
        self.drawz(setlim=False, setcol=False)

    def onTopoTxt(self, evnet):
        self.nbtopo = int(self.TopoTxt.GetValue())
        self.drawz(setlim=False, setcol=False)

    #  Methods of class

    def updateVariableZ(self, setlim=True):
        """ reload current variable depending on the time and plot it """

        # Level section
        if self.typSection == "XY":

            try:
                dims = self.croco.variables[self.variableName].dims
            except Exception:
                dims = []

            if "x_u" in dims:
                mask = self.croco.umask()
                self.topo = self.croco.rho2u_2d(self.croco.wrapper.metrics["h"].values)
            elif "y_v" in dims:
                mask = self.croco.vmask()
                self.topo = self.croco.rho2v_2d(self.croco.wrapper.metrics["h"].values)
            else:
                mask = self.croco.rmask()
                self.topo = self.croco.wrapper.metrics["h"]
            mask = np.where(mask == 0., np.nan, mask)

            # Level plot
            if self.sliceCoord > 0:
                self.slice = "Level"
                if self.variableName in self.croco.ListOfVariables:
                    # 3D variable rho level
                    try:
                        self.variableZ = self.croco.variables[self.variableName]\
                            .isel(t=self.timeIndex, z_r=self.sliceIndex)
                    except Exception:
                        # 3D variable w level
                        try:
                            self.variableZ = self.croco.variables[self.variableName]\
                                .isel(t=self.timeIndex, z_w=self.sliceIndex)
                        # 2D variable
                        except Exception:
                            self.variableZ = self.croco.variables[self.variableName].isel(t=self.timeIndex)

                elif self.variableName in self.croco.ListOfDerived:
                    if 'pv' in self.variableName:
                        var = get_pv(self.croco, self.timeIndex, depth=self.sliceCoord, typ=self.variableName)
                        self.variableZ = xr.DataArray(data=var)
                    elif self.variableName == 'zeta_k':
                        var = get_zetak(self.croco, self.timeIndex, depth=self.sliceCoord)
                        self.variableZ = xr.DataArray(data=var)
                    elif self.variableName == 'dtdz':
                        try:
                            self.croco.variables['temp']
                            var = get_dtdz(self.croco, self.timeIndex, depth=self.sliceCoord)
                            self.variableZ = xr.DataArray(data=var)
                        except Exception:
                            print("No variable temp")
                            return
                    elif self.variableName == 'log(Ri)':
                            var = get_richardson(self.croco, self.timeIndex, depth=self.sliceCoord)
                            self.variableZ = xr.DataArray(data=var)
                else:
                    print("unknown variable ", self.variableName)
                    return

            # Depth plot
            elif self.sliceCoord <= 0:
                self.slice = "Depth"
                # Calculate depths
                z = self.croco.get_coord(self.variableName, direction='z', timeIndex=self.timeIndex)
                minlev, maxlev = self.croco.zslice(None, mask, z, self.sliceCoord, findlev=True)
                # Variable from croco file
                if self.variableName in self.croco.ListOfVariables:
                    try:
                        var = self.croco.variables[self.variableName]\
                            .isel(t=self.timeIndex, z_r=slice(minlev, maxlev + 1))
                    except Exception:
                        var = self.croco.variables[self.variableName]\
                            .isel(t=self.timeIndex, z_w=slice(minlev, maxlev + 1))
                    try:
                        zslice = self.croco.zslice(var.values, mask, z[minlev:maxlev + 1, :, :], self.sliceCoord)[0]
                        self.variableZ = xr.DataArray(data=zslice)
                    except Exception:
                        print("Not enough points")

                # Derived variable
                elif self.variableName in self.croco.ListOfDerived:
                    if 'pv' in self.variableName:
                        var = get_pv(self.croco, self.timeIndex, depth=self.sliceCoord, minlev=minlev,
                                     maxlev=maxlev, typ=self.variableName)
                        try:
                            zslice = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], self.sliceCoord)[0]
                            self.variableZ = xr.DataArray(data=zslice)
                        except Exception:
                            print("Not enough points")
                            pass

                    elif self.variableName == 'zeta_k':
                        var = get_zetak(self.croco, self.timeIndex, depth=self.sliceCoord, minlev=minlev, maxlev=maxlev)
                        try:
                            zslice = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], self.sliceCoord)[0]
                            self.variableZ = xr.DataArray(data=zslice)
                        except Exception:
                            print("Not enough points")
                            pass

                    elif self.variableName == 'dtdz':
                        var = get_dtdz(self.croco, self.timeIndex, depth=self.sliceCoord, minlev=minlev, maxlev=maxlev)
                        try:
                            zslice = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], self.sliceCoord)[0]
                            self.variableZ = xr.DataArray(data=zslice)
                        except Exception:
                            print("Not enough points")
                            pass

                    elif self.variableName == 'log(Ri)':
                        var = get_richardson(self.croco, self.timeIndex, depth=self.sliceCoord, minlev=minlev, maxlev=maxlev)
                        try:
                            zslice = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], self.sliceCoord)[0]
                            self.variableZ = xr.DataArray(data=zslice)
                        except Exception:
                            print("Not enough points")
                            pass
            # Draw the new self.variableZ
            self.variableZ.values = mask * self.variableZ.values

        # Latitude section
        elif self.typSection == "XZ":

            # Variable from croco file
            if self.variableName in self.croco.ListOfVariables:
                self.variableZ = self.croco.get_variable(self.variableName,
                                                         tindex=self.timeIndex,
                                                         yindex=self.sliceIndex)
            # Derived Variable
            elif self.variableName in self.croco.ListOfDerived:
                # ntimes = self.croco.wrapper.ntimes
                N = self.croco.wrapper.N
                L = self.croco.wrapper.L
                self.variableZ = self.croco.create_DataArray(data=np.full((N, L),np.nan),
                                                             dimstyp='xz')
                if 'pv' in self.variableName:
                    self.variableZ[1:, :] = get_pv(self.croco, self.timeIndex,
                                                   minlev=0,
                                                   maxlev=self.croco.wrapper.N - 1,
                                                   latindex=self.sliceIndex,
                                                   typ=self.variableName)

                elif self.variableName == 'zeta_k':
                    self.variableZ[1:, :] = get_zetak(self.croco, self.timeIndex,
                                                      minlev=0,
                                                      maxlev=self.croco.wrapper.N - 1,
                                                      latindex=self.sliceIndex)

                elif self.variableName == 'dtdz':
                    self.variableZ[1:, :] = get_dtdz(self.croco, self.timeIndex,
                                                     minlev=0,
                                                     maxlev=self.croco.wrapper.N - 1,
                                                     latindex=self.sliceIndex)

                elif self.variableName == 'log(Ri)':
                    self.variableZ[1:, :] = get_richardson(self.croco, self.timeIndex,
                                                     minlev=0,
                                                     maxlev=self.croco.wrapper.N - 1,
                                                     latindex=self.sliceIndex)

        # Longitude section
        elif self.typSection == "YZ":

            # Variable from croco file
            if self.variableName in self.croco.ListOfVariables:
                self.variableZ = self.croco.get_variable(self.variableName,
                                                         tindex=self.timeIndex,
                                                         xindex=self.sliceIndex)
            # Derived Variable
            elif self.variableName in self.croco.ListOfDerived:
                # ntimes = self.croco.wrapper.ntimes
                N = self.croco.wrapper.N
                M = self.croco.wrapper.M
                self.variableZ = self.croco.create_DataArray(data=np.full((N, M),np.nan),
                                                             dimstyp='yz')
                if 'pv' in self.variableName:
                    self.variableZ[1:, :] = get_pv(self.croco, self.timeIndex,
                                                   minlev=0,
                                                   maxlev=self.croco.wrapper.N - 1,
                                                   lonindex=self.sliceIndex,
                                                   typ=self.variableName)

                elif self.variableName == 'zeta_k':
                    self.variableZ[1:, :] = get_zetak(self.croco, self.timeIndex,
                                                      minlev=0,
                                                      maxlev=self.croco.wrapper.N - 1,
                                                      lonindex=self.sliceIndex)

                elif self.variableName == 'dtdz':
                    self.variableZ[1:, :] = get_dtdz(self.croco, self.timeIndex,
                                                     minlev=0,
                                                     maxlev=self.croco.wrapper.N - 1,
                                                     lonindex=self.sliceIndex)

                elif self.variableName == 'log(Ri)':
                    self.variableZ[1:, :] = get_richardson(self.croco, 
                                            self.timeIndex,
                                            minlev=0,
                                            maxlev=self.croco.wrapper.N - 1,
                                            lonindex=self.sliceIndex)

    def drawz(self, setlim=False, setcol=False, anim=False):
        """ plot the current variable in the canvas """

        # Don't plot if variable full of Nan
        if np.count_nonzero(~np.isnan(self.variableZ.values)) == 0: return

        self.figure.clf()
        self.canvas.mpl_connect('button_press_event', self.onFigureClick)
        self.canvas.mpl_connect('motion_notify_event', self.ShowPosition)

        # if self.variableZ.values is None: return
        variableZ = ma.masked_invalid(self.variableZ.values)
        if setcol:
            if self.variableName == 'log(Ri)':
                self.mincolor = -3.2
                self.maxcolor = 2.    
            else:
                self.mincolor = np.min(variableZ)
                self.maxcolor = np.max(variableZ)
            self.MinColorTxt.SetValue('%.2E' % self.mincolor)
            self.MaxColorTxt.SetValue('%.2E' % self.maxcolor)
            self.clim = [self.mincolor, self.maxcolor]
        if setlim:
            self.xlim = [np.min(self.x), np.max(self.x)]
            self.ylim = [np.min(self.y), np.max(self.y)]
        if self.typSection == "XY" and self.toggletopo == True:
            topo = self.topo
        else:
            topo = None
        time = str(self.croco.wrapper._get_date(self.timeIndex))
        self.title = "{:s}, {:s}={:4.1f}, Time={:s}".\
            format(self.variableName, self.slice, self.sliceCoord, time)
        mypcolor(self, self.x, self.y, variableZ,
                 title=self.title,
                 xlabel=self.xlabel,
                 ylabel=self.ylabel,
                 xlim=self.xlim,
                 ylim=self.ylim,
                 clim=self.clim,
                 topo=topo,
                 nbtopo=self.nbtopo)

        if anim:
            self.movie.grab_frame()

        self.canvas.draw()
        self.canvas.Refresh()
        self.Show()

# end of SectionFrame Class
########################################################################


class ProfileFrame(wx.Frame):
    """
    Window class to plot time series or depth profile.
    The window contains a canvas to plot, several buttons (zoom in, zoom out
    and print )

    Attributes:
    croco : croco instance to study
    x : numpy array of x coordinates
    y : numpy array of y coordinates
    variableName : Name of the variable to plot in the canvas
    title : title of the plot
    xlabel : label of x axis
    ylabel : label of y axis
    """

    def __init__(self, croco=None, x=None, y=None,
                 variableName=None, title=None, xlabel=None, ylabel=None):

        # Create the window
        wx.Frame.__init__(self, None, wx.ID_ANY, title='Profile')

        # Now create the Panel to put the other controls on.
        self.panel = wx.Panel(self, wx.ID_ANY)

        # and a few controls
        self.figure = Figure()
        self.canvas = FigureCanvas(self.panel, -1, self.figure)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Hide()

        self.ZoomBtn = wx.Button(self.panel, wx.ID_ANY, "Zoom")
        self.HomeBtn = wx.Button(self.panel, wx.ID_ANY, "Home")
        self.PanBtn = wx.Button(self.panel, wx.ID_ANY, "Pan")
        self.SavePlotBtn = wx.Button(self.panel, wx.ID_ANY, "Save Plot")

        # bind the menu event to an event handler
        self.ZoomBtn.Bind(wx.EVT_BUTTON, self.onZoomBtn)
        self.HomeBtn.Bind(wx.EVT_BUTTON, self.onHomeBtn)
        self.PanBtn.Bind(wx.EVT_BUTTON, self.onPanBtn)
        self.SavePlotBtn.Bind(wx.EVT_BUTTON, self.onSavePlotBtn)

        self.showPosition = self.CreateStatusBar(2)
        self.showPosition.SetStatusText("x=   , y=  ", 1)
        self.showPosition.SetStatusWidths([-1, 150])

        self.__do_layout()

        # Initialize the variables of the class
        self.croco = croco
        self.x = x
        self.y = y
        self.variableName = variableName
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel

    def __do_layout(self):

        """
        Use a sizer to layout the controls, stacked vertically or horizontally
        """
        topSizer = wx.BoxSizer(wx.VERTICAL)
        canvasSizer = wx.BoxSizer(wx.VERTICAL)
        buttonsSizer = wx.BoxSizer(wx.HORIZONTAL)

        canvasSizer.Add(self.canvas, 0, wx.ALL, 5)
        buttonsSizer.Add(self.ZoomBtn, 0, wx.ALL, 5)
        buttonsSizer.Add(self.PanBtn, 0, wx.ALL, 5)
        buttonsSizer.Add(self.HomeBtn, 0, wx.ALL, 5)
        buttonsSizer.Add(self.SavePlotBtn, 0, wx.ALL, 5)

        topSizer.Add(canvasSizer, 0, wx.CENTER)
        topSizer.Add(buttonsSizer, 0, wx.ALL | wx.EXPAND, 5)

        self.panel.SetSizer(topSizer)
        topSizer.Fit(self)

        self.Layout()

    # ------------ Event handler

    def rect_select_callback(self, eclick, erelease):
        """Event handler for rectangle selector on plot"""
        self.xPress, self.yPress = eclick.xdata, eclick.ydata
        self.xRelease, self.yRelease = erelease.xdata, erelease.ydata
        self.xlim = [min(self.xPress, self.xRelease), max(self.xPress, self.xRelease)]
        self.ylim = [min(self.yPress, self.yRelease), max(self.yPress, self.yRelease)]
        self.draw(setlim=False)

    def ShowPosition(self, event):
        if event.inaxes:
            self.showPosition.SetStatusText(
                "x={:5.1f}  y={:5.1f}".format(event.xdata, event.ydata), 1)

    def onZoomBtn(self, event):
        """
        Event handler for the button click Zoom in button
        """
        # self.figure.RS.set_active(True)
        self.toolbar.zoom()

    def onHomeBtn(self, event):
        """
        Event handler for the button click Zoom out button
        """
        # self.draw()
        self.toolbar.home()

    def onPanBtn(self, event):
        """
        Event handler for the button click Zoom out button
        """
        # self.draw()
        self.toolbar.pan()

    def onSavePlotBtn(self, event):
        """
        Event handler for the button click Print button
        """
        printDir = self.croco.startDir + "/Figures_" \
                                       + self.croco.get_run_name() + "/"
        if not os.path.isdir(printDir):
                os.mkdir(printDir)
        filename = self.title.replace(",", "_").replace(" ", "") + ".png"
        os.system('rm -Rf ' + printDir + filename)
        self.figure.savefig(printDir + filename, dpi=self.figure.dpi)
        # self.toolbar.save_figure(None)

    # Methods of class

    def draw(self, setlim=True):
        """ plot the current variable in the canvas """

        # Don't plot if variable full of Nan
        if np.count_nonzero(~np.isnan(self.x)) == 0 or \
           np.count_nonzero(~np.isnan(self.y)) == 0: return

        self.canvas.mpl_connect('motion_notify_event', self.ShowPosition)

        self.x = ma.masked_invalid(self.x)
        self.y = ma.masked_invalid(self.y)
        if setlim:
            self.xlim = [np.min(self.x), np.max(self.x)]
            self.ylim = [np.min(self.y), np.max(self.y)]
        title = self.title
        plotCurv(self, x=self.x, y=self.y, title=title, xlabel=self.xlabel,
                 ylabel=self.ylabel, xlim=self.xlim, ylim=self.ylim)
        self.canvas.draw()
        self.canvas.Refresh()
        self.Show()

# end of ProfileFrame Class
########################################################################


class CrocoGui(wx.Frame):
    """
    Window class to plot the XY sections, manage variables, times, levels and
    create other windows for vertical sections and profiles

    Attributes:
    title : name of the window
    """

    def __init__(self):

        # Create the window
        wx.Frame.__init__(self, None, wx.ID_ANY, title='Main Window')
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        # Now create the Panel to put the other controls on.
        self.Panel = wx.Panel(self, wx.ID_ANY)

        # and a few controls
        self.CrocoVariableChoice = wx.Choice(self.Panel, wx.ID_ANY,
                                             choices=["Croco Variables ..."])
        self.CrocoVariableChoice.SetSelection(0)

        self.DerivedVariableChoice = wx.Choice(self.Panel, wx.ID_ANY,
                                               choices=["Derived Variables ..."])
        self.DerivedVariableChoice.SetSelection(0)

        self.LabelTime = wx.StaticText(self.Panel, -1, label="Choose Time",
                                       style=wx.ALIGN_CENTER)
        self.LabelMinMaxTime = wx.StaticText(self.Panel, wx.ID_ANY, " ", style=wx.ALIGN_LEFT)
        self.TimeTxt = wx.TextCtrl(self.Panel, wx.ID_ANY, "Time", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)

        self.LabelLevel = wx.StaticText(self.Panel, -1, label="Choose level (level>0, depth<0)",
                                        style=wx.ALIGN_CENTER)
        self.LabelMinMaxLevel = wx.StaticText(self.Panel, wx.ID_ANY, " ", style=wx.ALIGN_LEFT)
        self.LabelMinMaxDepth = wx.StaticText(self.Panel, wx.ID_ANY, " ", style=wx.ALIGN_LEFT)
        self.LevSectionBtn = wx.Button(self.Panel, wx.ID_ANY, "Level/Depth Section")
        self.LevSectionTxt = wx.TextCtrl(self.Panel, wx.ID_ANY, "Level", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.LonSectionBtn = wx.Button(self.Panel, wx.ID_ANY, "Longitude Section")
        self.LonSectionTxt = wx.TextCtrl(self.Panel, wx.ID_ANY, "Longitude", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.LatSectionBtn = wx.Button(self.Panel, wx.ID_ANY, "Latitude Section")
        self.LatSectionTxt = wx.TextCtrl(self.Panel, wx.ID_ANY, "Latitude", style=wx.TE_CENTRE | wx.TE_PROCESS_ENTER)
        self.TimeSeriesBtn = wx.Button(self.Panel, wx.ID_ANY, "Time Series")
        self.VerticalProfileBtn = wx.Button(self.Panel, wx.ID_ANY, "Vertical Profile")

        # bind the menu event to an event handler
        self.CrocoVariableChoice.Bind(wx.EVT_CHOICE, self.onCrocoVariableChoice)
        self.DerivedVariableChoice.Bind(wx.EVT_CHOICE, self.onDerivedVariableChoice)
        self.TimeTxt.Bind(wx.EVT_TEXT_ENTER, self.onTimeTxt)
        self.LevSectionBtn.Bind(wx.EVT_BUTTON, self.onLevSectionBtn)
        self.LonSectionBtn.Bind(wx.EVT_BUTTON, self.onLonSectionBtn)
        self.LatSectionBtn.Bind(wx.EVT_BUTTON, self.onLatSectionBtn)
        self.TimeSeriesBtn.Bind(wx.EVT_BUTTON, self.onTimeSeriesBtn)
        self.VerticalProfileBtn.Bind(wx.EVT_BUTTON, self.onVerticalProfileBtn)

        # self.__set_properties()
        self.__do_layout()

        # Subscribe for future publication of sections. To retrieve notifications 
        # when location (level/longitude/latitude) change after a click on a plot
        self.sub = Subscriber(self)

        # Initialize the croco instance
        self.openCroco()

    def __do_layout(self):

        """
        Use a sizer to layout the controls, stacked vertically or horizontally
        """

        topSizer = wx.BoxSizer(wx.HORIZONTAL)
        leftSizer = wx.BoxSizer(wx.VERTICAL)
        rightSizer = wx.BoxSizer(wx.VERTICAL)
        chooseVariablesSizer = wx.BoxSizer(wx.HORIZONTAL)
        labelTimeSizer = wx.BoxSizer(wx.HORIZONTAL)
        timeSizer = wx.BoxSizer(wx.HORIZONTAL)
        labelLevelSizer = wx.BoxSizer(wx.HORIZONTAL)
        labelMinMaxLevelSizer = wx.BoxSizer(wx.HORIZONTAL)
        labelMinMaxDepthSizer = wx.BoxSizer(wx.HORIZONTAL)
        levelSizer = wx.BoxSizer(wx.HORIZONTAL)
        longitudeSizer = wx.BoxSizer(wx.HORIZONTAL)
        latitudeSizer = wx.BoxSizer(wx.HORIZONTAL)
        timeSeriesSizer = wx.BoxSizer(wx.HORIZONTAL)
        profileSizer = wx.BoxSizer(wx.HORIZONTAL)

        chooseVariablesSizer.Add(self.CrocoVariableChoice, 0, wx.ALL, 5)
        chooseVariablesSizer.Add(self.DerivedVariableChoice, 0, wx.ALL, 5)

        labelTimeSizer.Add(self.LabelMinMaxTime, 0, wx.ALL | wx.EXPAND, 5)
        timeSizer.Add(self.LabelTime, 0, wx.ALL, 5)
        timeSizer.Add(self.TimeTxt, 0, wx.ALL, 5)

        labelLevelSizer.Add(self.LabelLevel, 0, wx.ALL | wx.EXPAND, 5)
        labelMinMaxLevelSizer.Add(self.LabelMinMaxLevel, 0, wx.ALL | wx.EXPAND, 5)
        labelMinMaxDepthSizer.Add(self.LabelMinMaxDepth, 0, wx.ALL | wx.EXPAND, 5)

        levelSizer.Add(self.LevSectionBtn, 0, wx.ALL, 5)
        levelSizer.Add(self.LevSectionTxt, 0, wx.ALL, 5)

        longitudeSizer.Add(self.LonSectionBtn, 0, wx.ALL, 5)
        longitudeSizer.Add(self.LonSectionTxt, 0, wx.ALL, 5)

        latitudeSizer.Add(self.LatSectionBtn, 0, wx.ALL, 5)
        latitudeSizer.Add(self.LatSectionTxt, 0, wx.ALL, 5)

        timeSeriesSizer.Add(self.TimeSeriesBtn, 0, wx.ALL, 5)

        profileSizer.Add(self.VerticalProfileBtn, 0, wx.ALL, 5)

        leftSizer.Add(chooseVariablesSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(labelTimeSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(timeSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(labelLevelSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(labelMinMaxLevelSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(labelMinMaxDepthSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(levelSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(longitudeSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(latitudeSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(timeSeriesSizer, 0, wx.ALL | wx.EXPAND, 5)
        leftSizer.Add(profileSizer, 0, wx.ALL | wx.EXPAND, 5)

        topSizer.Add(leftSizer, 0, wx.ALL | wx.EXPAND, 5)

        self.Panel.SetSizer(topSizer)
        self.Panel.SetAutoLayout(True)
        topSizer.Fit(self)

        self.Layout()

    # ------------ Event handler

    def update(self, level=None, longitude=None, latitude=None):
        """
        method to be applied when a section publisher send a new location 
        (level/longitude/latitude) after a mouse click the  plot
        """
        if level is not None:
            if level > 0:
                self.LevSectionTxt.SetValue('%d' % level)
            else:
                self.LevSectionTxt.SetValue('%.1F' % level)
        if longitude is not None:
            self.LonSectionTxt.SetValue('%.2F' % longitude)
        if latitude is not None:
            self.LatSectionTxt.SetValue('%.2F' % latitude)

    def OnClose(self, event):
        self.Destroy()
        sys.exit()

    def openCroco(self):
        """
        Create and show the Open FileDialog to select file name
        Initialize few outputs
        """

        startDir = os.getcwd()
        self.croco = Croco()
        self.croco.startDir = startDir

        # Fill the different text of the main window
        timeMin = self.croco.wrapper._get_date(0)
        timeMax = self.croco.wrapper._get_date(self.croco.wrapper.ntimes - 1)
        self.LabelMinMaxTime.SetLabel("Min/Max Time = " + str(timeMin) + " ... " +
                                      str(timeMax) + " days")
        self.TimeTxt.SetValue(str(timeMin))
        self.timeIndex = 0
        self.time = timeMin
        # minLevel = 1
        maxLevel = int(self.croco.wrapper.N)
        minDepth = - int(self.croco.wrapper.metrics['h'].max())
        maxDepth = 0
        self.LabelMinMaxLevel.SetLabel("Min/Max Level = 1 ... " + str(maxLevel))
        self.LabelMinMaxDepth.SetLabel("Min/Max Depth = " + str(minDepth) + " ... " + str(maxDepth))
        self.levelIndex = self.croco.wrapper.N - 1
        self.LevSectionTxt.SetValue(str(self.levelIndex + 1))
        self.depth = self.levelIndex + 1
        self.CrocoVariableChoice.AppendItems(self.croco.ListOfVariables)
        self.DerivedVariableChoice.AppendItems(self.croco.ListOfDerived)

        lon = self.croco.get_coord('ssh', direction='x')
        self.lon = lon[int(0.5 * self.croco.wrapper.M), int(0.5 * self.croco.wrapper.L)]
        lat = self.croco.get_coord('ssh', direction='y')
        self.lat = lat[int(0.5 * self.croco.wrapper.M), int(0.5 * self.croco.wrapper.L)]
        self.latIndex, self.lonIndex = self.findLatLonIndex(self.lon, self.lat)
        self.LonSectionTxt.SetValue('%.2F' % self.lon)
        self.LatSectionTxt.SetValue('%.2F' % self.lat)

    def onCrocoVariableChoice(self, event):
        ''' Choose variable from croco file to plot '''
        self.variableName = self.CrocoVariableChoice.GetString(self.CrocoVariableChoice.GetSelection())
        self.DerivedVariableChoice.SetSelection(0)

    def onDerivedVariableChoice(self, event):
        ''' Choose a computed variable to plot '''
        self.variableName = self.DerivedVariableChoice.GetString(self.DerivedVariableChoice.GetSelection())
        self.CrocoVariableChoice.SetSelection(0)

    def onTimeTxt(self, event):
        time = float(self.TimeTxt.GetValue())
        times = self.croco.get_coord(self.variableName, direction='t')
        # find index corresponding to the nearest instant time to plot
        self.timeIndex = min(range(len(times)), key=lambda j: abs(time - times[j]))
        self.time = self.croco.wrapper._get_date(self.timeIndex)
        self.TimeTxt.SetValue(str(self.time))

    def onLevSectionBtn(self, event):

        # if 2D variable, reset level to N
        try:
            dims = self.croco.variables[self.variableName].dims
        except Exception:
            dims = []
        if len(dims) == 3:
                self.level = self.croco.wrapper.N
                self.LevSectionTxt.SetValue('%d' % self.level)

        self.level = float(self.LevSectionTxt.GetValue())
        if self.level > 0:
            self.levelIndex = min(int(self.level - 1), self.croco.wrapper.N - 1)
            # self.LevSectionTxt.SetValue(str(self.levelIndex + 1))
        else:
            z = self.croco.get_coord(self.variableName, direction='z', timeIndex=self.timeIndex)
            self.levelIndex = np.argmax(z[:, self.latIndex, self.lonIndex] >= self.level)
        self.drawz(typSection="XY")

    def onLonSectionBtn(self, event):
        # if variable without z dimension
        try:
            ndims = len(self.croco.variables[self.variableName].dims)
        except Exception:
            ndims = 4
        if ndims < 4:
            print("Not 3D variable")
            return
        self.lon = float(self.LonSectionTxt.GetValue())
        # Find nearest indices of selected point
        self.latIndex, self.lonIndex = self.findLatLonIndex(self.lon, self.lat)
        self.drawz(typSection="YZ")

    def onLatSectionBtn(self, event):
        # if variable without z dimension
        try:
            ndims = len(self.croco.variables[self.variableName].dims)
        except Exception:
            ndims = 4
        if ndims < 4:
            print("Not 3D variable")
            return
        self.lat = float(self.LatSectionTxt.GetValue())
        # Find nearest indices of selected point
        self.latIndex, self.lonIndex = self.findLatLonIndex(self.lon, self.lat)
        self.drawz(typSection="XZ")

    def onTimeSeriesBtn(self, event):
        depth = float(self.LevSectionTxt.GetValue())

        # Get the mask at the rigth point
        try:
            dims = self.croco.variables[self.variableName].dims
        except Exception:
            dims = []

        if "x_u" in dims:
            mask = self.croco.umask()
        elif "y_v" in dims:
            mask = self.croco.vmask()
        else:
            mask = self.croco.rmask()
        # mask = np.where(mask == 0., np.nan, mask)

        # Get x coordinate: time
        x = self.croco.get_coord(self.variableName, direction='t').astype('timedelta64[D]').astype('float')

        # Time series on level
        if depth > 0:
            # Variable from croco file
            if self.variableName in self.croco.ListOfVariables:
                y = self.croco.get_variable(self.variableName, xindex=self.lonIndex,
                                            yindex=self.latIndex, zindex=self.levelIndex).values

            # Derived variable
            elif self.variableName in self.croco.ListOfDerived:
                y = np.full_like(x, np.nan)
                for it in range(len(x)):
                    if 'pv' in self.variableName:
                        y[it] = get_pv(self.croco, it, depth=self.levelIndex,
                                       typ=self.variableName)[self.latIndex, self.lonIndex]
                    elif self.variableName == 'zeta_k':
                        y[it] = get_zetak(self.croco, it, depth=self.levelIndex)[self.latIndex, self.lonIndex]
                    elif self.variableName == 'dtdz':
                        y[it] = get_dtdz(self.croco, it, depth=self.levelIndex)[self.latIndex, self.lonIndex]
                    elif self.variableName == 'log(Ri)':
                        y[it] = get_richardson(self.croco, it, depth=self.levelIndex)[self.latIndex, self.lonIndex]

            title = "{:s}, Lon={:4.1f}, Lat={:4.1f}, Level={:4.0f}".\
                format(self.variableName, self.lon, self.lat, self.depth)

        # Time series on depth
        else:

            y = np.zeros_like(x)
            # recalculate the depth slice at each time step
            for it in range(len(x)):
                # Calculate z
                z = self.croco. get_coord(self.variableName, direction='z', timeIndex=self.timeIndex)

                if self.variableName in self.croco.ListOfVariables:
                    # Find levels around depth
                    maxlev = np.argmax(z[:, self.latIndex, self.lonIndex] >= depth)
                    minlev = maxlev - 1
                    z1 = z[minlev, self.latIndex, self.lonIndex]
                    z2 = z[maxlev, self.latIndex, self.lonIndex]
                    # read variable and do interpolation
                    try:
                        var = self.croco.variables[self.variableName].isel(t=it,
                              z_r=slice(minlev, maxlev + 1))[:, self.latIndex, self.lonIndex]
                    except Exception:
                        var = self.croco.variables[self.variableName].isel(t=it,
                              z_w=slice(minlev,maxlev+1))[:,self.latIndex,self.lonIndex]
                    y[it] = ((var[0] - var[1]) * depth + var[1] * z1 - var[0] * z2) / (z1 - z2)

                elif self.variableName in self.croco.ListOfDerived:
                    # Find all the level corresponding to depth
                    minlev, maxlev = self.croco.zslice(None, mask, z, depth, findlev=True)
                    if 'pv' in self.variableName:
                        # Compute pv between these levels
                        var = get_pv(self.croco, it, depth=depth, minlev=minlev, maxlev=maxlev,
                                     typ=self.variableName)
                        # Extract the slice of pv corresponding at the depth
                        try:
                            varz = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], depth)[0]
                        except Exception:
                            print("Not enough points")
                            pass
                        y[it] = varz[self.latIndex, self.lonIndex]

                    elif self.variableName == 'zeta_k':
                        # Compute pv between these levels
                        var = get_zetak(self.croco, it, depth=depth, minlev=minlev, maxlev=maxlev)
                        # Extract the slice of pv corresponding at the depth
                        try:
                            varz = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], depth)[0]
                        except Exception:
                            print("Not enough points")
                            pass
                        y[it] = varz[self.latIndex, self.lonIndex]

                    elif self.variableName == 'dtdz':
                        # Compute pv between these levels
                        var = get_dtdz(self.croco, it, depth=depth, minlev=minlev, maxlev=maxlev)
                        # Extract the slice of pv corresponding at the depth
                        try:
                            varz = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], depth)[0]
                        except Exception:
                            print("Not enough points")
                            pass
                        y[it] = varz[self.latIndex, self.lonIndex]

                    elif self.variableName == 'log(Ri)':
                        # Compute pv between these levels
                        var = get_richardson(self.croco, it, depth=depth, minlev=minlev, maxlev=maxlev)
                        # Extract the slice of pv corresponding at the depth
                        try:
                            varz = self.croco.zslice(var, mask, z[minlev:maxlev, :, :], depth)[0]
                        except Exception:
                            print("Not enough points")
                            pass
                        y[it] = varz[self.latIndex, self.lonIndex]

            title = "{:s}, Lon={:4.1f}, Lat={:4.1f}, depth={:4.1f}".format(self.variableName, self.lon, self.lat, depth)

        # Plot the time series
        self.timeFrame = ProfileFrame(croco=self.croco,
                                      x=x, y=y,
                                      variableName=self.variableName,
                                      title=title,
                                      xlabel="Time (days)")
        self.timeFrame.draw()

    def onVerticalProfileBtn(self, event):
        # Dimension must have z coordinate
        try:
            ndims = len(self.croco.variables[self.variableName].dims)
        except Exception:
            ndims = 4
        if ndims < 4:
            print("Not 3D variable")
            return

        time = str(self.croco.wrapper._get_date(self.timeIndex))
        title = "{:s}, Lon={:4.1f}, Lat={:4.1f}, Time={:s}".\
            format(self.variableName, self.lon, self.lat, time)
        # Get depths coordinate
        z = self.croco. get_coord(self.variableName, direction='z',
                                  timeIndex=self.timeIndex)[:, self.latIndex, self.lonIndex]

        # Get variable profile
        if self.variableName in self.croco.ListOfVariables:
            x = self.croco.get_variable(self.variableName,
                                        xindex=self.lonIndex, yindex=self.latIndex,
                                        tindex=self.timeIndex).values

        elif self.variableName in self.croco.ListOfDerived:
            if 'pv' in self.variableName:
                x = np.full_like(z, np.nan)
                var = get_pv(self.croco, self.timeIndex,
                             minlev=0, maxlev=self.croco.wrapper.N - 1,
                             lonindex=self.lonIndex, typ=self.variableName)
                x[1:] = var[:, self.latIndex]

            elif self.variableName == 'zeta_k':
                x = np.full_like(z, np.nan)
                var = get_zetak(self.croco, self.timeIndex,
                                minlev=0, maxlev=self.croco.wrapper.N - 1,
                                lonindex=self.lonIndex,)
                x[1:] = var[:, self.latIndex]

            elif self.variableName == 'dtdz':
                x = np.full_like(z, np.nan)
                var = get_dtdz(self.croco, self.timeIndex,
                               minlev=0, maxlev=self.croco.wrapper.N - 1,
                               lonindex=self.lonIndex,)
                x[1:] = var[:, self.latIndex]

            elif self.variableName == 'log(Ri)':
                x = np.full_like(z, np.nan)
                var = get_richardson(self.croco, self.timeIndex,
                               minlev=0, maxlev=self.croco.wrapper.N - 1,
                               lonindex=self.lonIndex,)
                x[1:] = var[:, self.latIndex]

        # Plot the profile
        self.profileFrame = ProfileFrame(croco=self.croco,
                                         x=x, y=z,
                                         variableName=self.variableName,
                                         title=title,
                                         ylabel="Depth (m)")
        self.profileFrame.draw()

    def notify(self, ax):
        self.xlim = ax.get_xlim()
        self.ylim = ax.get_ylim()

    def findLatLonIndex(self, lonValue, latValue):
        ''' Find nearest  grid point of  click value '''
        a = abs(self.croco.wrapper.coords['lon_r'] - lonValue) + \
            abs(self.croco.wrapper.coords['lat_r'] - latValue)
        return np.unravel_index(a.argmin(), a.shape)

    def drawz(self, typSection=None):
        '''
        Extract the  rigth section for the current variable and
        plot in a new window
        '''
        # Level/Depth section
        if typSection == "XY":

            # Get Latitude coordinates
            y = self.croco.get_coord(self.variableName, direction='y')

            # Get Longitude coordinates
            x = self.croco.get_coord(self.variableName, direction='x')
            # x = repmat(x[self.latIndex, :].squeeze(), y.shape[0], 1)

            # Create new window
            self.sectionXY = SectionFrame(croco=self.croco,
                                          variableName=self.variableName,
                                          # variable=variable,
                                          x=x, y=y,
                                          typSection="XY",
                                          sliceCoord=self.level,
                                          sliceIndex=self.levelIndex,
                                          timeIndex=self.timeIndex,
                                          subscriber=self.sub)
            # Draw the plot
            self.sectionXY.updateVariableZ()
            self.sectionXY.drawz(setlim=True, setcol=True)

        # Latitude section
        elif typSection == "XZ":

            # Get depths coordinates
            y = self.croco.get_coord(self.variableName, direction='z',
                                     timeIndex=self.timeIndex)[:, self.latIndex, :]
            # Get Longitude coordinates
            x = self.croco.get_coord(self.variableName, direction='x')
            x = repmat(x[self.latIndex, :].squeeze(), y.shape[0], 1)

            # Create new window
            self.sectionXZ = SectionFrame(croco=self.croco,
                                          variableName=self.variableName,
                                          # variable=variable,
                                          x=x, y=y,
                                          typSection="XZ",
                                          sliceCoord=self.lat,
                                          sliceIndex=self.latIndex,
                                          timeIndex=self.timeIndex,
                                          subscriber=self.sub)
            # Draw the plot
            self.sectionXZ.updateVariableZ()
            self.sectionXZ.drawz(setlim=True, setcol=True)

        # Longitude section
        elif typSection == "YZ":

            # Get depths coordinates
            y = self.croco.get_coord(self.variableName, direction='z',
                                     timeIndex=self.timeIndex)[:, :, self.lonIndex]
            # Get Latitude coordinates
            x = self.croco.get_coord(self.variableName, direction='y')
            x = repmat(x[:, self.lonIndex].squeeze(), y.shape[0], 1)
            # Create new window
            self.sectionYZ = SectionFrame(croco=self.croco,
                                          variableName=self.variableName,
                                          # variable=variable,
                                          x=x, y=y,
                                          typSection="YZ",
                                          sliceCoord=self.lon,
                                          sliceIndex=self.lonIndex,
                                          timeIndex=self.timeIndex,
                                          subscriber=self.sub)
            # Draw the plot
            self.sectionYZ.updateVariableZ()
            self.sectionYZ.drawz(setlim=True, setcol=True)


# end of class CrocoGui


# Run the program
if __name__ == "__main__":
    import ctypes
    if sys.platform.startswith('linux'):
        try:
            print(ctypes.x)
            sys.exit()
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except Exception:
            print("Warning: failed to XInitThreads()")
    app = wx.App(False)
    frame = CrocoGui()
    frame.Show()
    app.MainLoop()
