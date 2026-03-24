"""
In this file, one can find the used workers.
"""

import sys

from PySide6.QtGui import QFont
from PySide6.QtCore import QObject, Signal, Slot
import pyqtgraph as pg
import numpy as np
import Code_Gui.Gui_General_Code.General_Functions_Library as GFL
import Code_Gui.Gui_General_Code.Gas_Mixtures_Spectra_Library as GMSL
import h5py
import os
import os.path as op
import matplotlib as mpl
from matplotlib.figure import Figure
from lmfit import Model
from lmfit import minimize, Parameters, fit_report
import matplotlib.pyplot as plt
from radis.io.hitran import fetch_hitran


class WorkerPlotter(QObject):
    """
    Worker used for plotting.
    """
    def __init__(self, parent):
        super().__init__()
        self.dict_plots = {}
        self.parent = parent
        self.font_label = QFont("Times", 11)
        self.font_title = QFont("Times", 16)

    @Slot()
    def create_empty_plot(self, name):
        """
        Function that is capable of creating a directory where the different graphs are stored and shown within our GUI.
        Using a directory also makes it easier to call upon the created graphs later.

        :param name: The name we give to each plot within our "saving"-directory
        :return:
        """
        if name not in self.dict_plots.keys():
            self.dict_plots[name] = pg.GraphicsLayoutWidget()
        else:
            del self.dict_plots[name]
            self.dict_plots[name] = pg.GraphicsLayoutWidget()

        return self.dict_plots[name]

    @Slot()
    def single_transmission_plot(self, data):
        """
        This function creates a single data plot. Can not be used to make multiple

        :param data: a list that contains both the name of the graph as well as the data one wants to plot
        :return: -
        """
        name = data[0]

        self.dict_plots[name].p1 = self.dict_plots[name].addPlot(0, 0)

        # Set plot parameters for plot 1
        self.dict_plots[name].p1.getViewBox().setMouseMode(pg.ViewBox.RectMode)
        self.dict_plots[name].p1.setTitle("Simulated Spectra")
        self.dict_plots[name].p1.titleLabel.item.setFont(self.font_title)
        self.dict_plots[name].p1.setLabel('left', "Intensity / (A.U.)")
        self.dict_plots[name].p1.getAxis('left').label.setFont(self.font_label)
        self.dict_plots[name].p1.setLabel('bottom', " Wavenumber / (cm^-1)")
        self.dict_plots[name].p1.getAxis('bottom').label.setFont(self.font_label)
        self.dict_plots[name].p1.pen = pg.mkPen("red")
        self.dict_plots[name].p1.legend = self.dict_plots[name].p1.addLegend()
        self.dict_plots[name].p1.legend.mouseDragEvent = lambda *args, **kwargs: None
        self.dict_plots[name].p1.legend.anchor((1, 0), (1, 0))

        # Plot the data as well as save the x- and y-data
        self.dict_plots[name].p1.x = data[1]
        if self.parent.tab_ftir_simulator.absorbance_bool:
            self.dict_plots[name].p1.y = GFL.tr_to_ab(data[2])
        else:
            self.dict_plots[name].p1.y = data[2]
        self.dict_plots[name].p1.enableAutoRange('y', enable=True)
        self.dict_plots[name].p1.dataline = \
            self.dict_plots[name].p1.plot(self.dict_plots[name].p1.x, self.dict_plots[name].p1.y,
                                          pen=self.dict_plots[name].p1.pen)

    @Slot()
    def update_single_transmission_plot(self, data):
        """
        Function that updates previously created graphs, by rewriting the data inside of it. We also change the range of
        the y-range, so that all the data is present and as clearly visible as possible.

        :param data: a list that contains both the name of the graph as well as the data one wants to plot
        :return: -
        """
        name = ""
        if len(data) == 2:
            name = data[0]
            if data[1]:
                self.dict_plots[name].p1.y = GFL.tr_to_ab(self.dict_plots[name].p1.y)
            else:
                self.dict_plots[name].p1.y = GFL.ab_to_tr(self.dict_plots[name].p1.y)
        if len(data) == 3:
            name = data[0]
            self.dict_plots[name].p1.x = data[1]
            if self.parent.tab_ftir_simulator.absorbance_bool:
                self.dict_plots[name].p1.y = GFL.tr_to_ab(data[2])
            else:
                self.dict_plots[name].p1.y = data[2]

        self.dict_plots[name].p1.dataline.setData(self.dict_plots[name].p1.x, self.dict_plots[name].p1.y)
        self.dict_plots[name].p1.enableAutoRange('y', enable=True)

    @Slot()
    def update_fitting_and_residual_plot(self, data):
        """
        Function that updates previously created graphs, by rewriting the data inside of it. We also change the range of
        the y-range, so that all the data is present and as clearly visible as possible.

        :param data: a list that contains both the name of the graph as well as the data one wants to plot
        :return: -
        """

        name = data[0]

        if len(data) == 2:
            if data[1]:
                self.dict_plots[name].p1.y_exp = GFL.tr_to_ab(self.dict_plots[name].p1.y_exp)
                self.dict_plots[name].p1.y_fit = GFL.tr_to_ab(self.dict_plots[name].p1.y_fit)
                if all(i == 0 for i in self.dict_plots[name].p2.y_res):
                    self.dict_plots[name].p2.y_res = self.dict_plots[name].p2.y_res
                else:
                    self.dict_plots[name].p2.y_res = \
                        (self.dict_plots[name].p1.y_fit - self.dict_plots[name].p1.y_exp).tolist()
            else:
                self.dict_plots[name].p1.y_exp = GFL.ab_to_tr(self.dict_plots[name].p1.y_exp)
                self.dict_plots[name].p1.y_fit = GFL.ab_to_tr(self.dict_plots[name].p1.y_fit)
                if all(i == 0 for i in self.dict_plots[name].p2.y_res):
                    self.dict_plots[name].p2.y_res = self.dict_plots[name].p2.y_res
                else:
                    self.dict_plots[name].p2.y_res = \
                        (self.dict_plots[name].p1.y_fit - self.dict_plots[name].p1.y_exp).tolist()

        if len(data) == 5:
            self.dict_plots[name].p1.x_fit = data[1]
            if self.parent.tab_ftir_fitting.absorbance_bool:
                self.dict_plots[name].p1.y_exp = GFL.tr_to_ab(data[2])
                self.dict_plots[name].p1.y_fit = GFL.tr_to_ab(data[3])
                self.dict_plots[name].p2.y_res = \
                    (self.dict_plots[name].p1.y_fit - self.dict_plots[name].p1.y_exp).tolist()
            else:
                self.dict_plots[name].p1.y_exp = data[2]
                self.dict_plots[name].p1.y_fit = data[3]
                self.dict_plots[name].p2.y_res = data[4]

        self.dict_plots[name].p1.dataline_exp.setData(self.dict_plots[name].p1.x_fit, self.dict_plots[name].p1.y_exp)
        self.dict_plots[name].p1.dataline_fit.setData(self.dict_plots[name].p1.x_fit, self.dict_plots[name].p1.y_fit)
        self.dict_plots[name].p2.dataline_res.setData(self.dict_plots[name].p1.x_fit, self.dict_plots[name].p2.y_res)
        self.dict_plots[name].p1.enableAutoRange('y', enable=True)
        self.dict_plots[name].p2.enableAutoRange('y', enable=True)

    @Slot()
    def fitting_and_residual_plot(self, data):
        """
        This function creates a single figure, split into two plots. Into the first, both the experimental data as well
        as the theoretical fit to this experimental data is plotted. Into the second, the residual is plotted, to show
        how big the error is.

        :param data: a list that contains both the name of the graph as well as the data (for exp, fit and res) one
        wants to plot
        :return: -
        """

        name = data[0]
        x = data[1]
        y = data[2]

        # Create plot for 1)fitting and 2)residual
        self.dict_plots[name].p1 = self.dict_plots[name].addPlot(0, 0)
        self.dict_plots[name].p2 = self.dict_plots[name].addPlot(1, 0)
        self.dict_plots[name].p1.setXLink(self.dict_plots[name].p2)

        # Set needed parameters for plot 1
        self.dict_plots[name].p1.getViewBox().setMouseMode(pg.ViewBox.RectMode)
        self.dict_plots[name].p1.setTitle("Experimental Spectra vs Fitted Spectra")
        self.dict_plots[name].p1.titleLabel.item.setFont(self.font_title)
        self.dict_plots[name].p1.setLabel('left', "Intensity / (A.U.)")
        self.dict_plots[name].p1.getAxis('left').label.setFont(self.font_label)
        self.dict_plots[name].p1.enableAutoRange('y', enable=True)
        self.dict_plots[name].p1.legend = self.dict_plots[name].p1.addLegend()
        self.dict_plots[name].p1.legend.mouseDragEvent = lambda *args, **kwargs: None
        self.dict_plots[name].p1.legend.anchor((1, 0), (1, 0))

        # Set needed parameters for plot 2
        self.dict_plots[name].p2.setTitle("Residual")
        self.dict_plots[name].p2.titleLabel.item.setFont(self.font_title)
        self.dict_plots[name].p2.setLabel('left', "Intensity / (A.U.)")
        self.dict_plots[name].p2.getAxis('left').label.setFont(self.font_label)
        self.dict_plots[name].p2.setLabel('bottom', " Wavenumber / (cm^-1)")
        self.dict_plots[name].p2.getAxis('bottom').label.setFont(self.font_label)
        self.dict_plots[name].p2.getViewBox().setMouseMode(pg.ViewBox.RectMode)
        self.dict_plots[name].p2.enableAutoRange('y', enable=True)
        self.dict_plots[name].p2.legend = self.dict_plots[name].p2.addLegend()
        self.dict_plots[name].p2.legend.mouseDragEvent = lambda *args, **kwargs: None
        self.dict_plots[name].p2.legend.anchor((1, 0), (1, 0))

        self.dict_plots[name].p1.x_exp = x
        self.dict_plots[name].p1.x_fit = x
        self.dict_plots[name].p2.x_res = x

        if self.parent.tab_ftir_fitting.absorbance_bool:
            self.dict_plots[name].p1.y_exp = GFL.tr_to_ab(y)
            self.dict_plots[name].p1.y_fit = GFL.tr_to_ab(np.ones(len(self.dict_plots[name].p1.x_fit)))
            self.dict_plots[name].p2.y_res = GFL.tr_to_ab(np.zeros(len(self.dict_plots[name].p2.x_res)))
        else:
            self.dict_plots[name].p1.y_exp = y
            self.dict_plots[name].p1.y_fit = np.ones(len(self.dict_plots[name].p1.x_fit))*np.round(np.mean(y))
            self.dict_plots[name].p2.y_res = np.zeros(len(self.dict_plots[name].p2.x_res))

        self.dict_plots[name].p1.pen_exp = pg.mkPen("red")
        self.dict_plots[name].p1.pen_fit = pg.mkPen("blue")
        self.dict_plots[name].p2.pen_res = pg.mkPen("green")
        self.dict_plots[name].p1.dataline_exp = \
            self.dict_plots[name].p1.plot(self.dict_plots[name].p1.x_exp, self.dict_plots[name].p1.y_exp,
                                          name="experimental data", pen=self.dict_plots[name].p1.pen_exp)
        self.dict_plots[name].p1.dataline_fit = \
            self.dict_plots[name].p1.plot(self.dict_plots[name].p1.x_fit, self.dict_plots[name].p1.y_fit,
                                          name="fit", pen=self.dict_plots[name].p1.pen_fit)
        self.dict_plots[name].p2.dataline_res = \
            self.dict_plots[name].p2.plot(self.dict_plots[name].p2.x_res, self.dict_plots[name].p2.y_res,
                                          name="residual", pen=self.dict_plots[name].p2.pen_res)

class WorkerFTIRFitter(QObject):
    """
    Worker used for fitting
    """
    signal_fitting_plot = Signal(list)
    signal_save_plot_as_png = Signal(str)
    signal_save_data_in_txt = Signal(str)

    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.dict_molecules = {}
        self.dict_molecules_errors = {}


    def ftir_fit_one(self, s):
        self.tab = self.parent.tab_ftir_fitting
        self.inner_tab = self.parent.tab_ftir_fitting.inner_tab
        files_in_directory = self.inner_tab.files_in_directory
        print(files_in_directory)
        len_files_in_dir = len(files_in_directory)
        index = self.inner_tab.currentIndex()
        current_file = files_in_directory[index]
        print(self.tab.directory_save_InvenioR_processed)
        data_text_file_selected = self.tab.directory_save_InvenioR_processed + "/Data (in text files)" + "/" + \
                                  current_file + ".txt"
        data = self.inner_tab.array_data[index]
        list_molecules = self.tab.molecule_storage_for_fitting
        dict_molecules = {}
        dict_molecules_errors = {}

        if os.path.exists(data_text_file_selected) and not s[0] and os.stat(data_text_file_selected).st_size > 0:
            w_exp = []
            t_exp = []
            t_fit = []
            t_residual = []

            with open(data_text_file_selected, "r") as f:
                lines = f.readlines()
                for line_i in range(1, len(lines)):
                    w_exp.append(float(lines[line_i].split()[0]))
                    t_exp.append(float(lines[line_i].split()[1]))
                    t_fit.append(float(lines[line_i].split()[2]))
                    t_residual.append(float(lines[line_i].split()[3]))
            self.signal_fitting_plot.emit(["ftir_fitting_" + current_file, w_exp, t_exp, t_fit, t_residual])
            self.tab.layout.label_info_fit.setText("Fitting for " + current_file + " made")
            self.tab.layout.label_info_fit.setStyleSheet("background-color:white; font-size:11pt; font:bold; color:blue")


        else:
            bool_molecules = False
            bool_pressure = False
            bool_temperature = False

            if len(list_molecules) != 0:
                bool_molecules = True

            list_pressure = self.tab.list_ftir_pressure
            if len(list_pressure) == 1:
                pressure = list_pressure[0]
                list_pressure = np.zeros(len_files_in_dir)
                bool_pressure = True
            elif len(list_pressure) == len_files_in_dir and len(list_pressure) != 0:
                pressure = list_pressure[index]
                bool_pressure = True

            list_temperature = self.tab.list_ftir_temperature
            if len(list_temperature) == 1:
                temperature = list_temperature[0]
                list_temperature = np.zeros(len_files_in_dir)
                bool_temperature = True
            elif len(list_temperature) == len_files_in_dir and len(list_temperature) != 0:
                temperature = list_temperature[index]
                bool_temperature = True


            list_pathlength = self.tab.list_ftir_pathlength 
            if len(list_pathlength) == 1:
                pathlength = list_pathlength[0]
                list_pathlength = np.zeros(len_files_in_dir)
                bool_pathlength = True
            elif len(list_pathlength) == len_files_in_dir and len(list_pathlength) != 0:
                pathlength = list_pathlength[index]
                bool_pathlength= True


            if bool_molecules and bool_temperature and bool_pressure:
                if not op.exists(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5"):
                    hf = h5py.File(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5",
                                   "w")
                self.tab.layout.label_info_fit.setText("Fitting for "
                                                       + current_file
                                                       + " ....")
                self.tab.layout.label_info_fit.setStyleSheet("background-color:white; font-size:11pt; font:bold; color:blue")

                list_time = []
                hour_0 = 0
                minute_0 = 0
                second_0 = 0
                for i in range(0, len_files_in_dir):
                    try:
                        hour = self.inner_tab.files_in_directory[i][:2]
                        minute = self.inner_tab.files_in_directory[i][3:5]
                        second = self.inner_tab.files_in_directory[i][6:8]
                        if i == 0:
                            hour_0 = int(hour)
                            minute_0 = int(minute)
                            second_0 = int(second)
                            time_spend = 0
                            list_time.append(time_spend)
                        else:
                            time_spend = 60 * 60 * (int(hour) - hour_0) + 60 * (int(minute) - minute_0) + (
                                    int(second) - second_0)
                            list_time.append(time_spend)
                    except:
                        list_time.append(0)

                with h5py.File(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5", "r+") as hf:
                    keys = list(hf.keys())

                    if "File Names" not in keys:
                        hf.create_dataset("File Names", data=files_in_directory)
                    elif not len(hf["File Names"][:] == len_files_in_dir):
                        del hf["File Names"]
                        hf.create_dataset("File Names", data=files_in_directory)

                    if "Time" not in keys:
                        hf.create_dataset("Time", data=list_time)
                    elif not len(hf["Time"][:]) == len_files_in_dir:
                        del hf["Time"]
                        hf.create_dataset("Time", data=list_time)

                    if "Pressure" not in keys:
                        hf.create_dataset("Pressure", data=list_pressure)
                    elif not len(hf["Pressure"][:]) == len_files_in_dir:
                        del hf["Pressure"]
                        hf.create_dataset("Pressure", data=list_pressure)

                    if "Temperature" not in keys:
                        hf.create_dataset("Temperature", data=list_temperature)
                    elif not len(hf["Temperature"][:]) == len_files_in_dir:
                        del hf["Temperature"]
                        hf.create_dataset("Temperature", data=list_temperature)

                    if "Pathlength" not in keys:
                        hf.create_dataset("Pathlength", data=list_pathlength)
                    elif not len(hf["Pathlength"][:]) == len_files_in_dir:
                        del hf["Pathlength"]
                        hf.create_dataset("Pathlength", data=list_pathlength)

                    if "Slitsize" not in keys:
                        hf.create_dataset("Slitsize", data=np.ones(len_files_in_dir)*self.tab.slit_size)
                    elif not len(hf["Slitsize"][:]) == len_files_in_dir:
                        del hf["Slitsize"]
                        hf.create_dataset("Slitsize", data=np.ones(len_files_in_dir)*self.tab.slit_size)

                    if "k0" not in keys:
                        hf.create_dataset("k0", data=np.ones(len_files_in_dir)*self.tab.k0)
                    elif not len(hf["k0"][:]) == len_files_in_dir:
                        del hf["k0"]
                        hf.create_dataset("k0", data=np.ones(len_files_in_dir)*self.tab.k0)

                    if "k1" not in keys:
                        hf.create_dataset("k1", data=np.ones(len_files_in_dir)*self.tab.k1)
                    elif not len(hf["k1"][:]) == len_files_in_dir:
                        del hf["k1"]
                        hf.create_dataset("k1", data=np.ones(len_files_in_dir)*self.tab.k1)

                    if "wmin" not in keys:
                        hf.create_dataset("wmin", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_min)
                    elif not len(hf["wmin"][:]) == len_files_in_dir:
                        del hf["wmin"]
                        hf.create_dataset("wmin", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_min)

                    if "wmax" not in keys:
                        hf.create_dataset("wmax", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_max)
                    elif not len(hf["wmax"][:]) == len_files_in_dir:
                        del hf["wmax"]
                        hf.create_dataset("wmax", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_max)                        

                    for mol_i in range(len(list_molecules)):
                        if list_molecules[mol_i] not in list(dict_molecules.keys()):
                            dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                            dict_molecules_errors[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                        elif not len(dict_molecules[list_molecules[mol_i]]) == len_files_in_dir:
                            dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                            dict_molecules_errors[list_molecules[mol_i]] = np.zeros(len(files_in_directory))

                        if list_molecules[mol_i] + "_concentration" not in keys:
                            hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                              data=np.zeros(len(files_in_directory)))
                        elif not len(hf[list_molecules[mol_i] + "_concentration"]) == len_files_in_dir:
                            del hf[list_molecules[mol_i] + "_concentration"]
                            hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                              data=np.zeros(len(files_in_directory)))

                        if list_molecules[mol_i] + "_concentration_errors" not in keys:
                            hf.create_dataset(list_molecules[mol_i] + "_concentration_errors",
                                              data=np.zeros(len(files_in_directory)))
                        elif not len(hf[list_molecules[mol_i] + "_concentration_errors"]) == len_files_in_dir:
                            del hf[list_molecules[mol_i] + "_concentration_errors"]
                            hf.create_dataset(list_molecules[mol_i] + "_concentration_errors",
                                              data=np.zeros(len(files_in_directory)))


                try:
                    data_wavelength = np.array(data.get_range("AB")[0:-1])
                    data_transmission = np.array(data["AB"][0:-1])
                except:
                    data_wavelength = data[0]
                    data_transmission = data[1]

                self.wlmin = self.parent.tab_ftir_fitting.wavenumber_min
                self.wlmax = self.parent.tab_ftir_fitting.wavenumber_max
                self.wlmin_j = np.where(data_wavelength == GFL.take_closest(data_wavelength, self.wlmin))[0][0]
                self.wlmax_j = np.where(data_wavelength == GFL.take_closest(data_wavelength, self.wlmax))[0][0]


                if self.parent.tab_ftir_fitting.background_bool and self.inner_tab.background == "" :
                    self.inner_tab.background = files_in_directory[0]

                if self.inner_tab.background != "":
                    self.w_background, self.t_background = self.inner_tab.w_dict[files_in_directory[0]], \
                        self.inner_tab.t_dict[files_in_directory[0]]
                    if self.wlmax_j > self.wlmin_j:
                        self.w_background = self.w_background[self.wlmin_j:self.wlmax_j]
                        self.t_background = self.t_background[self.wlmin_j:self.wlmax_j]
                    else:
                        self.w_background = self.w_background[self.wlmax_j: self.wlmin_j]
                        self.t_background = self.t_background[self.wlmax_j: self.wlmin_j]
                    self.a_background = GFL.tr_to_ab(self.t_background)
                else:
                    if self.wlmax_j > self.wlmin_j:
                        self.a_background = np.zeros(len(data_wavelength))[self.wlmin_j:self.wlmax_j]
                    else:
                        self.a_background = np.zeros(len(data_wavelength))[self.wlmax_j: self.wlmin_j]

                if self.inner_tab.background != current_file:
                    print(self.inner_tab.w_dict.keys())
                    w_exp, t_exp = self.inner_tab.w_dict[current_file], self.inner_tab.t_dict[current_file]
                    print(t_exp)

                    if self.wlmax_j > self.wlmin_j:

                        w_exp = w_exp[self.wlmin_j:self.wlmax_j]
                        t_exp = t_exp[self.wlmin_j:self.wlmax_j]
                    else:
                        w_exp = w_exp[self.wlmax_j: self.wlmin_j]
                        t_exp = t_exp[self.wlmax_j: self.wlmin_j]
                    a_exp = GFL.tr_to_ab(t_exp)


                    bas = GFL.fit_baseline(w_exp, a_exp, self.a_background, self.tab.lnc_mct_detector_bool)


                    a_exp_baseline_corrected = a_exp - bas
                    t_exp_baseline_corrected = GFL.ab_to_tr(a_exp_baseline_corrected)
                    


                    dict_molecules = {}
                    for mol_i in range(len(list_molecules)):
                        dict_molecules[list_molecules[mol_i]] = 0.001
                    
                    dict_spec={}

                    fitParameters = Parameters()
                    fitParameters.add('k0', value=self.tab.k0, min=-5, max=5, vary=True)
                    #fitParameters.add('k1', value=self.tab.k1, vary=False)
                    
                    for molecule in dict_molecules.keys():

                        dict_spec[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule],
                                                                                        pres=pressure, temp=temperature,
                                                                                        path_l=pathlength,
                                                                                        wl_min=self.wlmin, wl_max=self.wlmax, step=0.001)
                        fitParameters.add('c_' + molecule, value=dict_molecules[molecule], min=0.001, max=1, vary=True)

                    GMSL.storage_for_dict(temperature, pressure, pathlength, self.tab.slit_size, self.tab.k1)

                    print(fitParameters)
                    
                    
                    out = minimize(GMSL.spectra_molecules_s, fitParameters, args=(w_exp, t_exp, dict_spec), method='leastq')
                    print(fit_report(out))
                    fit = GMSL.spectra_molecules_s(out.params, w_exp, t_exp, dict_spec, test=True)

                    dict_c_old = {}
                    dict_c_errors_old = {}
                    
                    for molecule in dict_molecules.keys():
                        if out.params['c_' + molecule].value < 0:
                            dict_c_old[molecule] = 0
                            dict_c_errors_old[molecule] = 0
                        else:
                            dict_c_old[molecule] = out.params['c_' + molecule].value
                            if out.params['c_' + molecule].stderr is not None:
                                dict_c_errors_old[molecule] = (out.params['c_' + molecule].stderr/out.params['c_' + molecule].value)
                            else:
                                dict_c_errors_old[molecule] = 0


                    # Re-fit if needed
                    refit_bool = True
                    for fit_i in range(5):
                        print("Refit: " + str(fit_i))
                        if refit_bool:
                            dict_spec_refit = {}
                            fitParameters_new = Parameters()
                            fitParameters_new.add('k0', value=out.params['k0'].value, min=-5, max=5, vary=True)
                            #fitParameters_new.add('k1', value=self.tab.k1, vary=False)

                            for molecule in dict_c_old.keys():
                                dict_spec_refit[molecule] = (
                                    GFL.spectrum_in_air_creator(mol=molecule,
                                                                mf=dict_c_old[molecule],
                                                                pres=pressure, temp=temperature,
                                                                path_l=pathlength,
                                                                wl_min=self.wlmin,
                                                                wl_max=self.wlmax,
                                                                step=0.001))
                                fitParameters_new.add('c_' + molecule, value=dict_c_old[molecule], min=0.8*dict_c_old[molecule], max=1.2*dict_c_old[molecule], vary=True)

                            GMSL.storage_for_dict(temperature, pressure, pathlength, self.tab.slit_size, self.tab.k1)

                            out_refit = minimize(GMSL.spectra_molecules_s, fitParameters_new, args=(w_exp, t_exp, dict_spec_refit), method='leastq')
                            print(fit_report(out_refit))
                            refit = GMSL.spectra_molecules_s(out_refit.params, w_exp, t_exp, dict_spec_refit, test=True)

                            dict_c_new = {}
                            dict_c_errors_new={}
                            for molecule in dict_c_old.keys():
                                if out_refit.params['c_' + molecule].value < 0:
                                    dict_c_new[molecule] = 0
                                    dict_c_errors_new[molecule] = 0
                                else:
                                    dict_c_new[molecule] = out_refit.params['c_' + molecule].value

                                    if out_refit.params['c_' + molecule].stderr is not None:
                                        dict_c_errors_new[molecule] = (
                                            out_refit.params['c_' + molecule].stderr/out_refit.params['c_' + molecule].value)
                                    else:
                                        dict_c_errors_new[molecule] = 0
                            
                            correctness_count = 0
                            for molecule in dict_molecules.keys():
                                if (np.abs(
                                        dict_c_new[molecule] - dict_c_old[molecule]) /
                                    dict_c_old[molecule]) > 0.002: # was 0.005 and < 1*10**-6
                                    correctness_count += 1

                                dict_c_old[molecule] = dict_c_new[molecule]
                                dict_c_errors_old[molecule] = dict_c_errors_new[molecule]
                            if correctness_count == 0:
                                refit_bool = False
                    residual = (refit - t_exp) / t_exp
                    for molecule in dict_molecules.keys():
                        dict_molecules[molecule] = dict_c_new[molecule]
                        dict_molecules_errors[molecule] = dict_c_errors_old[molecule]

                        if dict_molecules[molecule] < 1 * 10 ** -6:
                            dict_molecules[molecule] = 0
                            print("Mole fraction " + str(molecule) + ": " + str(0))
                        else:
                            dict_molecules[molecule] = dict_molecules[molecule]
                            print(
                                "Mole fraction " + str(molecule) + ": " + str(dict_molecules[molecule]))

                    self.signal_fitting_plot.emit(
                        ["ftir_fitting_" + current_file, w_exp, t_exp, refit, residual])
                    with (h5py.File(
                            self.parent.tab_ftir_fitting.directory_save_InvenioR_processed + "/" + "Measurement_file.h5",
                            "r+") as hf):
                        if hf["Temperature"][index] != temperature:
                            hf["Temperature"][index] = temperature
                        if hf["Pressure"][index] != pressure:
                            hf["Pressure"][index] = pressure
                        if hf['Slitsize'][index] != self.parent.tab_ftir_fitting.slit_size:
                            hf['Slitsize'][index] = self.parent.tab_ftir_fitting.slit_size
                        if hf['wmin'][index] != self.parent.tab_ftir_fitting.wavenumber_min:
                            hf['wmin'][index] = self.parent.tab_ftir_fitting.wavenumber_min
                        if hf['wmax'][index] != self.wlmax:
                            hf['wmax'][index] = self.wlmax
                        if hf['k0'][index] != self.tab.k0:
                            hf['k0'][index] = self.tab.k0
                        if hf['k1'][index] != self.tab.k1:
                            hf['k1'][index] = self.tab.k1


                        for molecule in dict_molecules.keys():
                            hf[molecule + "_concentration"][index] = dict_molecules[molecule]
                            hf[molecule + "_concentration_errors"][index] = \
                                dict_molecules_errors[molecule]

                    self.tab.layout.label_info_fit.setText("Fitting for " + current_file + " made")
                    self.tab.layout.label_info_fit.setStyleSheet("background-color:white; font-size:11pt; font:bold; color:green")
                else:
                    self.tab.layout.label_info_fit.setText("No fit, background-file")
                    self.tab.layout.label_info_fit.setStyleSheet("background-color:white; font-size:11pt; font:bold; color:black")
            else:
                self.tab.layout.label_info_fit.setText("Can't fit, "
                                                       "incorrect/no molecules, pressure, "
                                                       "temperature or pathlength")
                self.tab.layout.label_info_fit.setStyleSheet("background-color:red; font-size:11pt; font:bold")
    
    def ftir_fit_all(self, s):
        self.tab = self.parent.tab_ftir_fitting
        self.inner_tab = self.tab.inner_tab
        files_in_directory = self.inner_tab.files_in_directory
        len_files_in_dir = len(files_in_directory)
        list_molecules = self.tab.molecule_storage_for_fitting
        dict_molecules = {}
        dict_molecules_errors = {}
        all_fits_made_bool = False

        list_time = []
        hour_0 = 0
        minute_0 = 0
        second_0 = 0
        for i in range(0, len_files_in_dir):
            try:
                hour = self.inner_tab.files_in_directory[i][:2]
                minute = self.inner_tab.files_in_directory[i][3:5]
                second = self.inner_tab.files_in_directory[i][6:8]
                if i == 0:
                    hour_0 = int(hour)
                    minute_0 = int(minute)
                    second_0 = int(second)
                    time_spend = 0
                    list_time.append(time_spend)
                else:
                    time_spend = 60 * 60 * (int(hour) - hour_0) + 60 * (int(minute) - minute_0) + (
                            int(second) - second_0)
                    list_time.append(time_spend)
            except:
                list_time.append(0)

        bool_molecules = False
        bool_pressure = False
        bool_temperature = False
        bool_pathlength = False

        if len(list_molecules) != 0:
            bool_molecules = True

        list_pressure = self.tab.list_ftir_pressure
        if len(list_pressure) == 1:
            pressure = list_pressure[0]
            list_pressure = np.ones(len_files_in_dir) * pressure
            bool_pressure = True
        elif len(list_pressure) == len_files_in_dir and len(list_pressure) != 0:
            bool_pressure = True

        list_temperature = self.tab.list_ftir_temperature
        if len(list_temperature) == 1:
            temperature = list_temperature[0]
            list_temperature = np.ones(len_files_in_dir) * temperature
            bool_temperature = True
        elif len(list_temperature) == len_files_in_dir and len(list_temperature) != 0:
            bool_temperature = True

        list_pathlength = self.tab.list_ftir_pathlength
        if len(list_pathlength) == 1:
            pathlength = list_pathlength[0]
            list_pathlength = np.ones(len_files_in_dir) * pathlength
            bool_pathlength = True
        elif len(list_pathlength) == len_files_in_dir and len(list_pathlength) != 0:
            bool_pathlength = True

        self.wlmin = self.parent.tab_ftir_fitting.wavenumber_min
        self.wlmax = self.parent.tab_ftir_fitting.wavenumber_max

        data_temp = self.inner_tab.array_data[0]
        try:
            data_wavelength_temp = np.array(data_temp.get_range("AB")[0:-1])
        except:
            data_wavelength_temp = data_temp[0]

        self.wlmin_j = np.where(data_wavelength_temp == GFL.take_closest(data_wavelength_temp, self.wlmin))[0][0]
        self.wlmax_j = np.where(data_wavelength_temp == GFL.take_closest(data_wavelength_temp, self.wlmax))[0][0]

        if self.parent.tab_ftir_fitting.background_bool and self.inner_tab.background == "":
            self.inner_tab.background = files_in_directory[0]

        if self.inner_tab.background != "":
            self.w_background, self.t_background = self.inner_tab.w_dict[files_in_directory[0]], \
                self.inner_tab.t_dict[files_in_directory[0]]
            if self.wlmax_j > self.wlmin_j:
                self.w_background = self.w_background[self.wlmin_j:self.wlmax_j]
                self.t_background = self.t_background[self.wlmin_j:self.wlmax_j]
            else:
                self.w_background = self.w_background[self.wlmax_j: self.wlmin_j]
                self.t_background = self.t_background[self.wlmax_j: self.wlmin_j]
            self.a_background = GFL.tr_to_ab(self.t_background)
        else:
            if self.wlmax_j > self.wlmin_j:
                self.a_background = np.zeros(len(data_wavelength_temp))[self.wlmin_j:self.wlmax_j]
            else:
                self.a_background = np.zeros(len(data_wavelength_temp))[self.wlmax_j: self.wlmin_j]


        for index in range(len_files_in_dir):
            current_file = files_in_directory[index]
            if self.inner_tab.background != current_file:
                data_text_file_selected = self.tab.directory_save_InvenioR_processed + \
                                          "/Data (in text files)" + "/" + current_file + ".txt"
                data = self.inner_tab.array_data[index]

                pressure = list_pressure[index]
                temperature = list_temperature[index]
                pathlength = list_pathlength[index]

                if os.path.exists(data_text_file_selected) and not s[0] and os.stat(data_text_file_selected).st_size > 0:
                    w_exp = []
                    t_exp = []
                    t_fit = []
                    t_residual = []

                    with open(data_text_file_selected, "r") as f:
                        lines = f.readlines()
                        for line_i in range(1, len(lines)):
                            w_exp.append(float(lines[line_i].split()[0]))
                            t_exp.append(float(lines[line_i].split()[1]))
                            t_fit.append(float(lines[line_i].split()[2]))
                            t_residual.append(float(lines[line_i].split()[3]))
                    self.signal_fitting_plot.emit(
                        ["ftir_fitting_" + current_file, w_exp, t_exp, t_fit, t_residual])
                    self.tab.layout.label_info_fit.setText("Fitting for "
                                                           + current_file
                                                           + " made")
                    all_fits_made_bool = True

                else:
                    if bool_molecules == True and bool_pathlength == True and bool_temperature == True and bool_pressure == True:

                        self.tab.layout.label_info_fit.setText("Fitting for "
                                                               + current_file
                                                               + " ....")
                        if not op.exists(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5"):
                            hf = h5py.File(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5",
                                           "w")

                        with h5py.File(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5",
                                       "r+") as hf:
                            keys = list(hf.keys())
                            if "File Names" not in keys:
                                hf.create_dataset("File Names", data=files_in_directory)
                            elif not len(hf["File Names"][:] == len_files_in_dir):
                                del hf["File Names"]
                                hf.create_dataset("File Names", data=files_in_directory)

                            if "Time" not in keys:
                                hf.create_dataset("Time", data=list_time)
                            elif not len(hf["Time"][:]) == len_files_in_dir:
                                del hf["Time"]
                                hf.create_dataset("Time", data=list_time)

                            if "Pressure" not in keys:
                                hf.create_dataset("Pressure", data=list_pressure)
                            elif not len(hf["Pressure"][:]) == len_files_in_dir:
                                del hf["Pressure"]
                                hf.create_dataset("Pressure", data=list_pressure)

                            if "Temperature" not in keys:
                                hf.create_dataset("Temperature", data=list_temperature)
                            elif not len(hf["Pressure"][:]) == len_files_in_dir:
                                del hf["Temperature"]
                                hf.create_dataset("Temperature", data=list_temperature)

                            if "Pathlength" not in keys:
                                hf.create_dataset("Pathlength", data=list_pathlength)
                            elif not len(hf["Pathlength"][:]) == len_files_in_dir:
                                del hf["Pathlength"]
                                hf.create_dataset("Pathlength", data=list_pathlength)
                                
                            if "Slitsize" not in keys:
                                hf.create_dataset("Slitsize", data=np.ones(len_files_in_dir)*self.tab.slit_size)
                            elif not len(hf["Slitsize"][:]) == len_files_in_dir:
                                del hf["Slitsize"]
                                hf.create_dataset("Slitsize", data=np.ones(len_files_in_dir)*self.tab.slit_size)

                            if "k0" not in keys:
                                hf.create_dataset("k0", data=np.ones(len_files_in_dir)*self.tab.k0)
                            elif not len(hf["k0"][:]) == len_files_in_dir:
                                del hf["k0"]
                                hf.create_dataset("k0", data=np.ones(len_files_in_dir)*self.tab.k0)

                            if "k1" not in keys:
                                hf.create_dataset("k1", data=np.ones(len_files_in_dir)*self.tab.k1)
                            elif not len(hf["k1"][:]) == len_files_in_dir:
                                del hf["k1"]
                                hf.create_dataset("k1", data=np.ones(len_files_in_dir)*self.tab.k1)

                            if "wmin" not in keys:
                                hf.create_dataset("wmin", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_min)
                            elif not len(hf["wmin"][:]) == len_files_in_dir:
                                del hf["wmin"]
                                hf.create_dataset("wmin", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_min)

                            if "wmax" not in keys:
                                hf.create_dataset("wmax", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_max)
                            elif not len(hf["wmax"][:]) == len_files_in_dir:
                                del hf["wmax"]
                                hf.create_dataset("wmax", data=np.ones(len_files_in_dir)*self.parent.tab_ftir_fitting.wavenumber_max)   

                            for mol_i in range(len(list_molecules)):
                                if list_molecules[mol_i] not in list(dict_molecules.keys()):
                                    dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                    dict_molecules_errors[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                elif not len(dict_molecules[list_molecules[mol_i]]) == len_files_in_dir:
                                    dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                    dict_molecules_errors[list_molecules[mol_i]] = np.zeros(len(files_in_directory))

                                if list_molecules[mol_i] + "_concentration" not in keys:
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                      data=np.zeros(len(files_in_directory)))
                                elif not len(hf[list_molecules[mol_i] + "_concentration"]) == len_files_in_dir:
                                    del hf[list_molecules[mol_i] + "_concentration"]
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                      data=np.zeros(len(files_in_directory)))

                                if list_molecules[mol_i] + "_concentration_errors" not in keys:
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration_errors",
                                                      data=np.zeros(len(files_in_directory)))
                                elif not len(hf[list_molecules[mol_i] + "_concentration_errors"]) == len_files_in_dir:
                                    del hf[list_molecules[mol_i] + "_concentration_errors"]
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration_errors",
                                                      data=np.zeros(len(files_in_directory)))

                        w_exp, t_exp = self.inner_tab.w_dict[current_file], self.inner_tab.t_dict[current_file]

                        if self.wlmax_j > self.wlmin_j:
                            w_exp = w_exp[self.wlmin_j:self.wlmax_j]
                            t_exp = t_exp[self.wlmin_j:self.wlmax_j]
                        else:
                            w_exp = w_exp[self.wlmax_j: self.wlmin_j]
                            t_exp = t_exp[self.wlmax_j: self.wlmin_j]
                        a_exp = GFL.tr_to_ab(t_exp)

                        bas = GFL.fit_baseline(w_exp, a_exp, self.a_background, self.tab.lnc_mct_detector_bool)

                        a_exp_baseline_corrected = a_exp - bas
                        t_exp_baseline_corrected = GFL.ab_to_tr(a_exp_baseline_corrected)


                        for mol_i in range(len(list_molecules)):
                            dict_molecules[list_molecules[mol_i]][index] = 0.001
                        
                        dict_spec={}

                        fitParameters = Parameters()
                        if self.tab.k0 > 0:
                            fitParameters.add('k0', value=self.tab.k0, min=0.5*self.tab.k0, max=2*self.tab.k0, vary=False)
                        elif self.tab.k0 < 0:
                            fitParameters.add('k0', value=self.tab.k0, min=2*self.tab.k0, max=0.5*self.tab.k0, vary=False)
                    
                        for molecule in dict_molecules.keys():

                            dict_spec[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule][index],
                                                                                        pres=pressure, temp=temperature,
                                                                                        path_l=pathlength,
                                                                                        wl_min=self.wlmin, wl_max=self.wlmax, step=0.001)
                            fitParameters.add('c_' + molecule, value=dict_molecules[molecule][index], min=0.001, max=1, vary=True)

                        GMSL.storage_for_dict(temperature, pressure, pathlength, self.tab.slit_size, self.tab.k1)

                        print(fitParameters)
                        
                        
                        out = minimize(GMSL.spectra_molecules_s, fitParameters, args=(w_exp, t_exp_baseline_corrected, dict_spec), method='leastq')
                        print(fit_report(out))
                        fit = GMSL.spectra_molecules_s(out.params, w_exp, t_exp_baseline_corrected, dict_spec, test=True)

                        dict_c_old = {}
                        dict_c_errors_old = {}
                        
                        for molecule in dict_molecules.keys():
                            if out.params['c_' + molecule].value < 0:
                                dict_c_old[molecule] = 0
                                dict_c_errors_old[molecule] = 0
                            else:
                                dict_c_old[molecule] = out.params['c_' + molecule].value
                                if out.params['c_' + molecule].stderr is not None:
                                    dict_c_errors_old[molecule] = (out.params['c_' + molecule].stderr/out.params['c_' + molecule].value)
                                else:
                                    dict_c_errors_old[molecule] = 0
                        
                        # Re-fit if needed
                        refit_bool = True
                        for fit_i in range(5):
                            print("Refit: " + str(fit_i))
                            if refit_bool:
                                dict_spec_refit = {}
                                fitParameters_new = Parameters()
                                if self.tab.k0 > 0:
                                    fitParameters_new.add('k0', value=self.tab.k0, min=0.5*self.tab.k0, max=2*self.tab.k0, vary=False)
                                elif self.tab.k0 < 0:
                                    fitParameters_new.add('k0', value=self.tab.k0, min=2*self.tab.k0, max=0.5*self.tab.k0, vary=False)

                                for molecule in dict_molecules.keys():
                                    dict_spec_refit[molecule] = (
                                        GFL.spectrum_in_air_creator(mol=molecule,
                                                                    mf=dict_c_old[molecule],
                                                                    pres=pressure, temp=temperature,
                                                                    path_l=pathlength,
                                                                    wl_min=self.wlmin,
                                                                    wl_max=self.wlmax,
                                                                    step=0.001))
                                    fitParameters_new.add('c_' + molecule, value=dict_c_old[molecule], min=0.8*dict_c_old[molecule], max=1.2*dict_c_old[molecule], vary=True)

                                GMSL.storage_for_dict(temperature, pressure, pathlength, self.tab.slit_size, self.tab.k1)

                                out_refit = minimize(GMSL.spectra_molecules_s, fitParameters_new, args=(w_exp, t_exp_baseline_corrected, dict_spec_refit), method='leastq')
                                print(fit_report(out_refit))
                                refit = GMSL.spectra_molecules_s(out_refit.params, w_exp, t_exp_baseline_corrected, dict_spec_refit, test=True)

                                dict_c_new = {}
                                dict_c_errors_new={}
                                for molecule in dict_molecules.keys():
                                    if out_refit.params['c_' + molecule].value < 0:
                                        dict_c_new[molecule] = 0
                                        dict_c_errors_new[molecule] = 0
                                    else:
                                        dict_c_new[molecule] = out_refit.params['c_' + molecule].value

                                        if out_refit.params['c_' + molecule].stderr is not None:
                                            dict_c_errors_new[molecule] = (
                                                out_refit.params['c_' + molecule].stderr/out_refit.params['c_' + molecule].value)
                                        else:
                                            dict_c_errors_new[molecule] = 0
                                
                                correctness_count = 0
                                for molecule in dict_molecules.keys():
                                    if (np.abs(
                                            dict_c_new[molecule] - dict_c_old[molecule]) /
                                        dict_c_old[molecule]) > 0.002: # was 0.005 and < 1*10**-6
                                        correctness_count += 1

                                    dict_c_old[molecule] = dict_c_new[molecule]
                                    dict_c_errors_old[molecule] = dict_c_errors_new[molecule]
                                if correctness_count == 0:
                                    refit_bool = False
                        residual = (refit - t_exp_baseline_corrected) / t_exp_baseline_corrected
                        for molecule in dict_molecules.keys():
                            dict_molecules[molecule][index] = dict_c_new[molecule]
                            dict_molecules_errors[molecule][index] = dict_c_errors_old[molecule]

                            if dict_molecules[molecule][index] < 1 * 10 ** -6:
                                dict_molecules[molecule][index] = 0
                                print("Mole fraction " + str(molecule) + ": " + str(0))
                            else:
                                dict_molecules[molecule][index] = dict_molecules[molecule][index]
                                print(
                                    "Mole fraction " + str(molecule) + ": " + str(dict_molecules[molecule]))

                        self.signal_fitting_plot.emit(
                            ["ftir_fitting_" + current_file, w_exp, t_exp, refit, residual])
                        with h5py.File(self.tab.directory_save_InvenioR_processed + "/" + "Measurement_file.h5",
                                       "r+") as hf:
                            if hf["Temperature"][index] != temperature:
                                hf["Temperature"][index] = temperature
                            if hf["Pressure"][index] != pressure:
                                hf["Pressure"][index] = pressure
                            if hf['Slitsize'][index] != self.parent.tab_ftir_fitting.slit_size:
                                hf['Slitsize'][index] = self.parent.tab_ftir_fitting.slit_size
                            if hf['wmin'][index] != self.wlmin:
                                hf['wmin'][index] = self.wlmin
                            if hf['wmax'][index] != self.wlmax:
                                hf['wmax'][index] = self.wlmax
                            if hf['k0'][index] != self.tab.k0:
                                hf['k0'][index] = self.tab.k0
                            if hf['k1'][index] != self.tab.k1:
                                hf['k1'][index] = self.tab.k1
                            for molecule in dict_molecules.keys():
                                hf[molecule + "_concentration"][index] = dict_molecules[molecule][index]
                                hf[molecule + "_concentration_errors"][index] = dict_molecules_errors[molecule][index]

                        self.tab.layout.label_info_fit.setText("Fitting for " + current_file + " made")
                        all_fits_made_bool = True
                    else:
                        self.tab.layout.label_info_fit.setText("Can't fit for "
                                                               + current_file + ", ""incorrect/no molecules, pressure, "
                                                                                "temperature or pathlength")
                        all_fits_made_bool = False
                        break
            else:
                self.tab.layout.label_info_fit.setText("No fit, background-file")
        if all_fits_made_bool:
            self.tab.layout.label_info_fit.setText("Fitting for all files made")

class WorkerFTIRSimulation(QObject):
    """
    Worker used for simuling FTIR spectra
    """
    signal_ftir_simulation_create_plot = Signal(list)
    signal_ftir_simulation_update_plot = Signal(list)

    def __init__(self, parent):
        super().__init__()
        self.dict_molecules_old = {}

        self.dict_w_old = {}
        self.dict_t_old = {}
        self.dict_a_old = {}

        self.parent = parent
        self.temperature_old = 0
        self.temperature_new = 0
        self.pressure_old = 0
        self.pressure_new = 0
        self.pathlength_old = 0
        self.pathlength_new = 0
        self.slitsize_old = ""
        self.wlmin_old = 0
        self.wlmax_old = 0
        self.wlmin_new = 0
        self.wlmax_new = 0
        self.dict_molecules_new = {}
        self.bool_remake = False
        print(self.parent)

    @Slot()
    def start_simulation(self, s):
        try:  
            self.signal_ftir_simulation_create_plot.connect(self.parent.worker_plotting.single_transmission_plot)
            self.signal_ftir_simulation_update_plot.connect(self.parent.worker_plotting.update_single_transmission_plot)  
            self.dict_molecules_new = {}
            self.tab = self.parent.tab_ftir_simulator

            self.wlmin_new = self.tab.wavenumber_min_for_simulation
            self.wlmax_new = self.tab.wavenumber_max_for_simulation
            self.slitsize_new = self.tab.slitsize_for_simulation

            for mol_i in self.tab.molecule_storage_for_simulation_final.keys():
                self.dict_molecules_new[mol_i] = self.tab.molecule_storage_for_simulation_final[mol_i]
            self.temperature_new = self.tab.temperature_for_simulation
            self.pressure_new = self.tab.pressure_for_simulation
            self.pathlength_new = self.tab.pathlength_for_simulation

            self.dict_w_new = {}
            self.dict_t_new = {}
            self.dict_a_new = {}
            self.list_simulated_molecules = []
            self.dict_spectrum = {}
            for mol_i in self.dict_molecules_new.keys():
                self.list_simulated_molecules.append(mol_i)
                if mol_i not in self.dict_molecules_old.keys() or self.temperature_new != self.temperature_old \
                        or self.pressure_new != self.pressure_old or self.pathlength_new != self.pathlength_old \
                        or self.wlmin_new != self.wlmin_old or self.wlmax_new != self.wlmax_old \
                        or self.slitsize_new != self.slitsize_old:
                    self.dict_spectrum[mol_i] = GFL.spectrum_in_air_creator(mol=mol_i,
                                                                                  mf=self.dict_molecules_new[mol_i],
                                                                                  pres=self.pressure_new,
                                                                                  temp=self.temperature_new,
                                                                                  path_l=self.pathlength_new,
                                                                                  wl_min=self.wlmin_new,
                                                                                  wl_max=self.wlmax_new, step=0.002)
                    if self.slitsize_new != 0:
                        self.dict_spectrum[mol_i].apply_slit(self.slitsize_new, unit="cm-1",
                                                             norm_by="area", inplace=True, shape="gaussian")
                        self.dict_w_new[mol_i] = self.dict_spectrum[mol_i].get("transmittance")[0]
                        self.dict_t_new[mol_i] = self.dict_spectrum[mol_i].get("transmittance")[1]
                        self.dict_a_new[mol_i] = GFL.tr_to_ab(self.dict_t_new[mol_i])
                    else:
                        self.dict_w_new[mol_i] = self.dict_spectrum[mol_i].get("transmittance_noslit")[0]
                        self.dict_t_new[mol_i] = self.dict_spectrum[mol_i].get("transmittance_noslit")[1]
                        self.dict_a_new[mol_i] = GFL.tr_to_ab(self.dict_t_new[mol_i])
                    self.tab.layout.label_simulation.setText(
                        "Simulated data for: " + str(self.list_simulated_molecules))

                if mol_i in self.dict_molecules_old.keys() and self.temperature_new == self.temperature_old \
                        and self.pressure_new == self.pressure_old and self.pathlength_new == self.pathlength_old \
                        and self.wlmin_new == self.wlmin_old and self.wlmax_new == self.wlmax_old \
                        and self.slitsize_new == self.slitsize_old:
                    c = self.dict_molecules_new[mol_i] / self.dict_molecules_old[mol_i]
                    self.dict_w_new[mol_i] = self.dict_w_old[mol_i]
                    self.dict_t_new[mol_i] = GFL.trold_to_trnew(self.dict_t_old[mol_i], c)
                    self.dict_a_new[mol_i] = GFL.tr_to_ab(self.dict_t_new[mol_i])
                    self.tab.layout.label_simulation.setText(
                        "Simulated data for: " + str(self.list_simulated_molecules))

            list_molecules_for_deletion = []
            for mol_i in range(len(self.dict_molecules_old.keys())):
                if not list(self.dict_molecules_old.keys())[mol_i] in self.dict_w_new.keys():
                    list_molecules_for_deletion.append(list(self.dict_molecules_old.keys())[mol_i])
            for i in list_molecules_for_deletion:
                del self.dict_molecules_old[i]

            self.count_molecules = 0
            self.a_full = 0
            for mol_i in self.dict_molecules_new.keys():
                if self.count_molecules == 0:
                    self.w_full = self.dict_w_new[mol_i]
                    self.a_full = self.dict_a_new[mol_i]
                else:
                    self.a_full += self.dict_a_new[mol_i]
                self.dict_molecules_old[mol_i] = self.dict_molecules_new[mol_i]
                self.dict_w_old[mol_i] = self.dict_w_new[mol_i]
                self.dict_a_old[mol_i] = self.dict_a_new[mol_i]
                self.dict_t_old[mol_i] = self.dict_t_new[mol_i]
                self.count_molecules += 1

            self.t_full = GFL.ab_to_tr(self.a_full)
            self.temperature_old = self.temperature_new
            self.pressure_old = self.pressure_new
            self.pathlength_old = self.pathlength_new
            self.wlmin_old = self.wlmin_new
            self.wlmax_old = self.wlmax_new
            self.slitsize_old = self.slitsize_new


            if not self.bool_remake:
                self.bool_remake = True
                self.signal_ftir_simulation_create_plot.emit(["Simulated Data", self.w_full, self.t_full])
                #self.parent.worker_plotting.single_transmission_plot(["Simulated Data", self.w_full, self.t_full])
                self.tab.layout.label_simulation.setText("Simulation done")
                self.tab.layout.label_simulation.setStyleSheet("color: green")
            else:
                #self.parent.worker_plotting.update_single_transmission_plot(
                #    ["Simulated Data", self.w_full, self.t_full])
                self.signal_ftir_simulation_update_plot.emit(["Simulated Data", self.w_full, self.t_full])
                self.tab.layout.label_simulation.setText("Simulation done")
                self.tab.layout.label_simulation.setStyleSheet("color: green")

        except:
            self.bool_remake = False
            #self.parent.worker_plotting.update_single_transmission_plot(
            #    ["Simulated Data", range(100), np.zeros(100)])
            self.signal_ftir_simulation_update_plot.emit(["Simulated Data", range(100), np.zeros(100)])
            self.tab.layout.label_simulation.setText("Can't do simulation, not enough conditions given")
            self.tab.layout.label_simulation.setStyleSheet("color: red")

class WorkerSaver(QObject):
    """
    Worker used for saving data
    """
    def __init__(self, parent):
        super().__init__()
        self.directory = ""
        self.name = ""
        self.data = []
        self.parent = parent

        mpl.use('Qt5Agg')
        mpl.rcParams['axes.edgecolor'] = "white"
        mpl.rcParams['axes.labelcolor'] = "white"
        mpl.rcParams['axes.titlecolor'] = "white"
        mpl.rcParams['axes.facecolor'] = "black"
        mpl.rcParams['xtick.color'] = "white"
        mpl.rcParams['ytick.color'] = "white"
        mpl.rcParams['legend.labelcolor'] = "white"
        mpl.rcParams['figure.figsize'] = (22, 11.5)

    def save_data_in_txt(self, data):
        self.directory = data[0]
        self.name = data[1]

        if data[3] == "Fitted FTIR Data":
            folder = "Data (in text files)"
            if not os.path.exists(self.directory + "/" + folder):
                os.makedirs(self.directory + "/" + folder)
            self.data = data[2]
            np.savetxt(self.directory + "/" + folder + "/" + self.name + ".txt", self.data, delimiter='\t',
                       header="Wavelength \t\t\tExperimental Transmission \t   Fitted Transmission \t\t      "
                              "Residual")
            print(self.directory)
        elif data[3] == "Simulated Data":
            self.data = data[2]
            np.savetxt(self.directory + self.name + ".txt", self.data, delimiter='\t',
                       header="Wavelength \t\t\tSimulated data")

    def save_data_as_png(self, data):
        self.directory = data[0]
        self.name = data[1]
        self.data = data[2]


        if data[3] == "Fitted FTIR Data":
            folder = "Data (in png's)"
            if not os.path.exists(self.directory + "/" + folder):
                os.makedirs(self.directory + "/" + folder)

            if not self.parent.tab_ftir_fitting.absorbance_bool:
                self.state = "Transmission"
            elif self.parent.tab_ftir_fitting.absorbance_bool:
                self.state = "Absorbance"
            figure = Figure(constrained_layout=True, facecolor='black')
            axes1 = figure.add_subplot(211)
            axes2 = figure.add_subplot(212)
            axes1.set_title("Experimental Spectra vs Fitted Spectra")
            axes1.title.set_fontsize(20)
            axes1.set_xlabel("Wavenumber / (cm$^{-1}$)")
            axes1.xaxis.label.set_fontsize(16)
            axes1.set_ylabel("Intensity / (A.U.)")
            axes1.yaxis.label.set_fontsize(16)
            axes1.get_xaxis().set_visible(False)
            axes2.set_title("Experimental Spectra vs Fitted Spectra")
            axes2.title.set_fontsize(20)
            axes2.set_xlabel("Wavenumber / (cm$^{-1}$)")
            axes2.xaxis.label.set_fontsize(16)
            axes2.set_ylabel("Intensity / (A.U.)")
            axes2.yaxis.label.set_fontsize(16)

            x_temp = data[2][0]
            x_temp = list(reversed(x_temp))
            [wlmin_view, wlmax_view] = \
                self.parent.worker_plotting.dict_plots["ftir_fitting_" + self.name].p1.getViewBox().state["viewRange"][
                    0]
            wlmin = x_temp[0]
            wlmax = x_temp[-1]
            wlmin_i = 0
            wlmax_i = 0

            if wlmin_view > wlmin:
                wlmin = wlmin_view
                wlmin_i = x_temp.index(GFL.take_closest(x_temp, wlmin))
            if wlmax_view < wlmax:
                wlmax = wlmax_view
                wlmax_i = x_temp.index(GFL.take_closest(x_temp, wlmax))

            if wlmax_i == 0:
                x = x_temp[wlmin_i:]
                if self.parent.tab_ftir_simulator.absorbance_bool:
                    y_exp = list(reversed(GFL.tr_to_ab(data[2][1])))[wlmin_i:]
                    y_fit = list(reversed(GFL.tr_to_ab(data[2][2])))[wlmin_i:]
                    y_res = list(reversed(GFL.tr_to_ab(data[2][2]) - GFL.tr_to_ab(data[2][1])))[wlmin_i:]
                else:
                    y_exp = list(reversed(data[2][1]))[wlmin_i:]
                    y_fit = list(reversed(data[2][2]))[wlmin_i:]
                    y_res = list(reversed(data[2][3]))[wlmin_i:]
            else:
                x = x_temp[wlmin_i:wlmax_i]
                if self.parent.tab_ftir_simulator.absorbance_bool:
                    y_exp = list(reversed(GFL.tr_to_ab(data[2][1])))[wlmin_i:wlmax_i]
                    y_fit = list(reversed(GFL.tr_to_ab(data[2][2])))[wlmin_i:wlmax_i]
                    y_res = list(reversed(GFL.tr_to_ab(data[2][2]) - GFL.tr_to_ab(data[2][1])))[wlmin_i:wlmax_i]
                else:
                    y_exp = list(reversed(data[2][1]))[wlmin_i:wlmax_i]
                    y_fit = list(reversed(data[2][2]))[wlmin_i:wlmax_i]
                    y_res = list(reversed(data[2][3]))[wlmin_i:wlmax_i]

            plot1 = axes1.plot(x, y_exp, label="experimental data", linewidth=1, color="red")
            plot2 = axes1.plot(x, y_fit, label="fit", linewidth=1, color="blue")
            plot3 = axes2.plot(x, y_res, label="fit", linewidth=1, color="green")
            axes1.tick_params(axis='both', which='major', labelsize=14)
            axes1.tick_params(axis='both', which='minor', labelsize=12)
            axes2.tick_params(axis='both', which='major', labelsize=14)
            axes2.tick_params(axis='both', which='minor', labelsize=12)
            figure.savefig(self.directory + "/" + folder + "/" + self.name + "_" + self.state + ".png", dpi=100)

        elif data[3] == "Simulated Data":
            if not self.parent.tab_ftir_simulator.absorbance_bool:
                self.state = "Transmission"
            elif self.parent.tab_ftir_simulator.absorbance_bool:
                self.state = "Absorbance"

            figure = Figure(constrained_layout=True, facecolor='black')
            axes = figure.add_subplot(111)
            axes.set_title("Simulated Spectra")
            axes.title.set_fontsize(20)
            axes.set_xlabel("Wavenumber / (cm$^{-1}$)")
            axes.xaxis.label.set_fontsize(16)
            axes.set_ylabel("Intensity / (A.U.)")
            axes.yaxis.label.set_fontsize(16)

            x_temp = data[2][0].tolist()
            [wlmin_view, wlmax_view] = \
                self.parent.worker_plotting.dict_plots["Simulated Data"].p1.getViewBox().state["viewRange"][0]
            wlmin = x_temp[0]
            wlmax = x_temp[-1]
            wlmin_i = 0
            wlmax_i = 0

            if wlmin_view > wlmin:
                wlmin = wlmin_view
                wlmin_i = x_temp.index(GFL.take_closest(x_temp, wlmin))
            if wlmax_view < wlmax:
                wlmax = wlmax_view
                wlmax_i = x_temp.index(GFL.take_closest(x_temp, wlmax))
            if wlmax_i == 0:
                x = x_temp[wlmin_i:]
                y = data[2][1][wlmin_i:]
            else:
                x = x_temp[wlmin_i:wlmax_i]
                y = data[2][1][wlmin_i:wlmax_i]

            plot = axes.plot(x, y, label="simulated data", linewidth=1, color="red")
            axes.tick_params(axis='both', which='major', labelsize=14)
            axes.tick_params(axis='both', which='minor', labelsize=12)
            figure.savefig(self.directory + self.name + ".png", dpi=100)

class WorkerHitran(QObject):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
    
    
    @Slot()
    def fetch(self, molecule_list):
        gotten_string = ""
        for molecule in molecule_list:
            self.parent.tab_get_database.layout.label_info_to_get.setText("Getting Hitran Data for molecules: " + molecule)
            
            self.parent.tab_get_database.layout.label_info_to_get.setStyleSheet("font-weight: bold ; font-size:12pt; background-color: white")
            
            
            fetch_hitran(molecule, cache="regen", extra_params="all")

                    
            gotten_string += molecule + ", " 
            self.parent.tab_get_database.layout.label_info_gotten.setText("Got Data for molecules: " + gotten_string)
            self.parent.tab_get_database.layout.label_info_gotten.setStyleSheet("font-weight: bold ; font-size:12pt; background-color: blue")
        self.parent.tab_get_database.layout.label_info_gotten.setStyleSheet("font-weight: bold ; font-size:12pt; background-color: green")

