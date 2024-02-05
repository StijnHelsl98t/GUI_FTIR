"""
In this file, one can find the used workers.
"""

import sys

sys.path.append("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\Code Gui")

from PySide6.QtGui import QFont
from PySide6.QtCore import QObject, Signal
import pyqtgraph as pg
from radis import Spectrum
import numpy as np
import Code_Gui.Gui_General_Code.General_Functions_Library as GFL
import Code_Gui.Gui_General_Code.Gas_Mixtures_Spectra_Library as GMSL
import h5py
import os
import os.path as op
import matplotlib as mpl
from matplotlib.figure import Figure
from lmfit import Model


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

    def create_empty_plot(self, name):
        if name not in self.dict_plots.keys():
            self.dict_plots[name] = pg.GraphicsLayoutWidget()
        else:
            del self.dict_plots[name]
            self.dict_plots[name] = pg.GraphicsLayoutWidget()

        return self.dict_plots[name]

    def single_transmission_plot(self, data):
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

        self.dict_plots[name].p1.x = data[1]
        if self.parent.tab_ftir_simulator.absorbance_bool:
            self.dict_plots[name].p1.y = GFL.tr_to_ab(data[2])
        else:
            self.dict_plots[name].p1.y = data[2]
        self.dict_plots[name].p1.enableAutoRange('y', enable=True)
        self.dict_plots[name].p1.dataline = \
            self.dict_plots[name].p1.plot(self.dict_plots[name].p1.x, self.dict_plots[name].p1.y,
                                          pen=self.dict_plots[name].p1.pen)

    def update_single_transmission_plot(self, data):
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

    def update_fitting_and_residual_plot(self, data):
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

    def fitting_and_residual_plot(self, data):
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
        self.fitting_functions = []
        self.fitting_functions.append(GMSL.spectra_fit_1_molecule)
        self.fitting_functions.append(GMSL.spectra_fit_2_molecules)
        self.fitting_functions.append(GMSL.spectra_fit_3_molecules)
        self.fitting_functions.append(GMSL.spectra_fit_4_molecules)
        self.fitting_functions.append(GMSL.spectra_fit_5_molecules)
        self.fitting_functions.append(GMSL.spectra_fit_6_molecules)
        self.fitting_functions.append(GMSL.spectra_fit_7_molecules)
        self.fitting_functions.append(GMSL.spectra_fit_8_molecules)
        self.parent = parent
        self.dict_molecules = {}
        self.dict_molecules_errors = {}


    def ftir_fit_one(self, s):
        self.tab = self.parent.tab_ftir_fitting
        self.inner_tab = self.parent.tab_ftir_fitting.inner_tab
        files_in_directory = self.inner_tab.files_in_directory
        len_files_in_dir = len(files_in_directory)
        index = self.inner_tab.currentIndex()
        current_file = files_in_directory[index]
        data_text_file_selected = self.tab.directory_save_invenioR_processed + "\\\\Data (in text files)" + "\\\\" + \
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

            pathlength = self.tab.list_ftir_pathlength[0]

            if bool_molecules and bool_temperature and bool_pressure:
                if not op.exists(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5"):
                    hf = h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                                   "w")
                self.tab.layout.label_info_fit.setText("Fitting for "
                                                       + current_file
                                                       + " ....")

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

                with h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5", "r+") as hf:
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
                        elif list_molecules[mol_i] + "_concentration_errors" not in keys:
                            hf.create_dataset(list_molecules[mol_i] + "_concentration_errors",
                                              data=np.zeros(len(files_in_directory)))
                        elif not len(hf[list_molecules[mol_i] + "_concentration"]) == len_files_in_dir:
                            del hf[list_molecules[mol_i] + "_concentration"]
                            hf.create_dataset(list_molecules[mol_i] + "_concentration",
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
                    dict_molecules = {}
                    for mol_i in range(len(list_molecules)):
                        dict_molecules[list_molecules[mol_i]] = 0.001

                    spectrum_dictionary = {}
                    dict_w = {}
                    dict_t = {}
                    for mol_i in range(len(dict_molecules)):
                        molecule = list(dict_molecules.keys())[mol_i]

                        spectrum_dictionary[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule],
                                                                                    pres=pressure, temp=temperature,
                                                                                    path_l=pathlength,
                                                                                    wl_min=self.wlmin, wl_max=self.wlmax, step=0.001)
                        dict_w[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[0]
                        dict_t[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[1]

                    GMSL.storage_for_dict(dict_w, dict_t, pathlength, self.tab.slit_size, self.tab.offset_left,
                                          self.tab.offset_right)
                    model = Model(self.fitting_functions[len(list_molecules) - 1])
                    if len(list_molecules) == 1:
                        params = model.make_params(c1=1)
                    elif len(list_molecules) == 2:
                        params = model.make_params(c1=1,c2=1)
                    elif len(list_molecules) == 3:
                        params = model.make_params(c1=1,c2=1,c3=1)
                    elif len(list_molecules) == 4:
                        params = model.make_params(c1=1,c2=1,c3=1,c4=1)
                    elif len(list_molecules) == 5:
                        params = model.make_params(c1=1,c2=1,c3=1,c4=1,c5=1)
                    elif len(list_molecules) == 6:
                        params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1)
                    elif len(list_molecules) == 7:
                        params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1)
                    elif len(list_molecules) == 8:
                        params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1, c8=1)

                    for i in range(len(list_molecules)):
                        params['c' + str(i+1)].set(min=-0.001, max=1000)

                    result = model.fit(t_exp_baseline_corrected, params, w=w_exp)

                    list_c_old = {}
                    c_list_minus = {}
                    list_c_errors_old = {}
                    for i in range(len(list_molecules)):
                        if result.best_values['c' + str(i + 1)] * dict_molecules[list_molecules[i]] < 0:
                            list_c_old[list_molecules[i]] = 0
                            list_c_errors_old[list_molecules[i]] = 0
                            c_list_minus[list_molecules[i]] = result.best_values['c' + str(i + 1)] * dict_molecules[
                                list_molecules[i]]
                        else:
                            list_c_old[list_molecules[i]] = result.best_values['c' + str(i + 1)] * dict_molecules[
                                list_molecules[i]]

                            list_c_errors_old[list_molecules[i]] = (
                                    result.params['c' + str(i + 1)].stderr/result.best_values['c' + str(i + 1)])


                    # Re-fit if needed
                    refit_bool = True
                    for fit_i in range(5):
                        print(fit_i)
                        if refit_bool:
                            spectrum_dictionary_refit = {}
                            dict_w_refit = {}
                            dict_t_refit = {}
                            for mol_i in range(len(dict_molecules)):
                                molecule = list(dict_molecules.keys())[mol_i]
                                spectrum_dictionary_refit[molecule] = (
                                    GFL.spectrum_in_air_creator(mol=molecule,
                                                                mf=list_c_old[molecule],
                                                                pres=pressure, temp=temperature,
                                                                path_l=pathlength,
                                                                wl_min=self.wlmin,
                                                                wl_max=self.wlmax,
                                                                step=0.001))
                                dict_w_refit[molecule] = \
                                    spectrum_dictionary_refit[molecule].get("transmittance_noslit")[0]
                                dict_t_refit[molecule] = \
                                    spectrum_dictionary_refit[molecule].get("transmittance_noslit")[1]

                            GMSL.storage_for_dict(dict_w_refit, dict_t_refit, pathlength,self.tab.slit_size, self.tab.offset_left,
                                          self.tab.offset_right)

                            model_refit = Model(self.fitting_functions[len(list_molecules) - 1])
                            if len(list_molecules) == 1:
                                params_refit = model_refit.make_params(c1=1)
                            elif len(list_molecules) == 2:
                                params_refit = model_refit.make_params(c1=1, c2=1)
                            elif len(list_molecules) == 3:
                                params_refit = model_refit.make_params(c1=1, c2=1, c3=1)
                            elif len(list_molecules) == 4:
                                params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1)
                            elif len(list_molecules) == 5:
                                params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1)
                            elif len(list_molecules) == 6:
                                params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1)
                            elif len(list_molecules) == 7:
                                params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1)
                            elif len(list_molecules) == 8:
                                params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1,
                                                                       c8=1)

                            for i in range(len(list_molecules)):
                                params_refit['c' + str(i + 1)].set(min=0.8, max=1.2)

                            result_refit = model_refit.fit(t_exp, params_refit, w=w_exp)

                            list_c_new = {}
                            list_c_errors_new={}
                            for i in range(len(list_molecules)):
                                if result_refit.best_values['c' + str(i + 1)] * list_c_old[list_molecules[i]] < 0:
                                    list_c_new[list_molecules[i]] = 0
                                    list_c_errors_new[list_molecules[i]] = 0
                                else:
                                    list_c_new[list_molecules[i]] = result_refit.best_values['c' + str(i + 1)] * \
                                                                           list_c_old[list_molecules[i]]

                                    if result_refit.params['c' + str(i + 1)].stderr is not None:
                                        list_c_errors_new[list_molecules[i]] = (
                                            result_refit.params['c' + str(i + 1)].stderr/result_refit.best_values['c' + str(i + 1)])
                                    else:
                                        list_c_errors_new[list_molecules[i]] = 0
                            correctness_count = 0
                            for i in range(len(list_molecules)):
                                if (np.abs(
                                        list_c_new[list_molecules[i]] - list_c_old[list_molecules[i]]) /
                                    list_c_old[
                                        list_molecules[i]] > 0.005 or
                                    np.abs(list_c_new[list_molecules[i]] - list_c_old[
                                        list_molecules[i]])) < 1 * 10 ** -6:
                                    correctness_count += 1

                                list_c_old[list_molecules[i]] = list_c_new[list_molecules[i]]
                                print(list_c_errors_old, list_c_errors_new)
                                list_c_errors_old[list_molecules[i]] += list_c_errors_new[list_molecules[i]]
                            if correctness_count == 0:
                                refit_bool = False

                    spectrum_dictionary_final = {}
                    molecule_final = list(dict_molecules.keys())[0]
                    dict_w_final = {}
                    dict_t_final = {}
                    dict_a_final = {}
                    a_full = []
                    for mol_i in range(len(dict_molecules)):
                        molecule = list(dict_molecules.keys())[mol_i]
                        spectrum_dictionary_final[molecule] = GFL.spectrum_in_air_creator(mol=molecule,
                                                                                          mf=list_c_new[molecule],
                                                                                          pres=pressure, temp=temperature,
                                                                                          path_l=pathlength,
                                                                                          wl_min=self.wlmin, wl_max=self.wlmax,
                                                                                          step=0.001)
                        dict_w_final[molecule] = spectrum_dictionary_final[molecule].get("transmittance_noslit")[0]
                        dict_t_final[molecule] = spectrum_dictionary_final[molecule].get("transmittance_noslit")[1]
                        dict_a_final[molecule] = GFL.tr_to_ab(dict_t_final[molecule])
                        if mol_i == 0:
                            a_full = dict_a_final[molecule]
                        else:
                            a_full += dict_a_final[molecule]
                    t_full = GFL.ab_to_tr(a_full)

                    spec_new2 = Spectrum({"wavenumber": dict_w_final[molecule_final], "transmittance_noslit": t_full},
                                         wunit='cm-1',
                                         units={"transmittance_noslit": ""})

                    spec_new2.apply_slit(self.tab.slit_size, unit="cm-1", norm_by="area", inplace=True)

                    w_temp, t_temp = spec_new2.get("transmittance")
                    offset_a = (self.tab.offset_right - self.tab.offset_left) / (4000 - 1000)
                    offset_b = self.tab.offset_right - offset_a * 4000

                    list_offset_temp = offset_a * w_temp + offset_b

                    w_new_temp = np.zeros(len(w_temp))
                    for i in range(len(w_temp)):
                        w_new_temp[i] = w_temp[i] + list_offset_temp[i]
                    spec_new_final = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                                         units={"transmittance": ""}, conditions={"path_length": pathlength})
                    spec_new_final.resample(w_exp, inplace=True)
                    w_new2, t_new2 = spec_new_final.get("transmittance")
                    residual = t_exp_baseline_corrected - result_refit.best_fit
                    for mol_i in range(len(dict_molecules)):
                        molecule = list(dict_molecules.keys())[mol_i]
                        dict_molecules[molecule] = list_c_new[molecule]
                        dict_molecules_errors[molecule] = list_c_errors_old[molecule]*list_c_new[molecule]

                    for mol_i in range(len(list_molecules)):
                        if dict_molecules[list_molecules[mol_i]] < 1 * 10 ** -6:
                            dict_molecules[list_molecules[mol_i]] = 0
                            print("Mole fraction " + str(list_molecules[mol_i]) + ": " + str(0))
                        else:
                            dict_molecules[list_molecules[mol_i]] = dict_molecules[list_molecules[mol_i]]
                            print(
                                "Mole fraction " + str(list_molecules[mol_i]) + ": " + str(dict_molecules[list_molecules[mol_i]]))

                    self.signal_fitting_plot.emit(
                        ["ftir_fitting_" + current_file, w_new2, t_exp, result_refit.best_fit, residual])
                    with (h5py.File(
                            self.parent.tab_ftir_fitting.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                            "r+") as hf):
                        if hf["Temperature"][index] != temperature:
                            hf["Temperature"][index] = temperature
                        if hf["Pressure"][index] != pressure:
                            hf["Pressure"][index] = pressure

                        for mol_i in range(len(list_molecules)):
                            hf[list_molecules[mol_i] + "_concentration"][index] = dict_molecules[list_molecules[mol_i]]
                            hf[list_molecules[mol_i] + "_concentration_errors"][index] = \
                                dict_molecules_errors[list_molecules[mol_i]]

                    self.tab.layout.label_info_fit.setText("Fitting for " + current_file + " made")
                else:
                    self.tab.layout.label_info_fit.setText("No fit, background-file")
            else:
                self.tab.layout.label_info_fit.setText("Can't fit, "
                                                       "incorrect/no molecules, pressure, "
                                                       "temperature or pathlength")

    def ftir_fit_all(self, s):
        self.tab = self.parent.tab_ftir_fitting
        self.inner_tab = self.tab.inner_tab
        files_in_directory = self.inner_tab.files_in_directory
        len_files_in_dir = len(files_in_directory)
        list_molecules = self.tab.molecule_storage_for_fitting
        dict_molecules = {}
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
                data_text_file_selected = self.tab.directory_save_invenioR_processed + \
                                          "\\\\Data (in text files)" + "\\\\" + current_file + ".txt"
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
                        if not op.exists(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5"):
                            hf = h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                                           "w")

                        with h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
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

                            for mol_i in range(len(list_molecules)):
                                if list_molecules[mol_i] not in list(dict_molecules.keys()):
                                    dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                elif not len(dict_molecules[list_molecules[mol_i]]) == len_files_in_dir:
                                    dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                if list_molecules[mol_i] + "_concentration" not in keys:
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                      data=np.zeros(len(files_in_directory)))
                                elif not len(hf[list_molecules[mol_i] + "_concentration"]) == len_files_in_dir:
                                    del hf[list_molecules[mol_i] + "_concentration"]
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                      data=np.zeros(len(files_in_directory)))

                            for mol_i in range(len(list_molecules)):
                                if list_molecules[mol_i] not in list(dict_molecules.keys()):
                                    dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                elif not len(dict_molecules[list_molecules[mol_i]]) == len_files_in_dir:
                                    dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                if list_molecules[mol_i] + "_concentration" not in keys:
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                      data=np.zeros(len(files_in_directory)))
                                elif not len(hf[list_molecules[mol_i] + "_concentration"]) == len_files_in_dir:
                                    del hf[list_molecules[mol_i] + "_concentration"]
                                    hf.create_dataset(list_molecules[mol_i] + "_concentration",
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

                        spectrum_dictionary = {}
                        dict_w = {}
                        dict_t = {}
                        for mol_i in range(len(dict_molecules)):
                            molecule = list(dict_molecules.keys())[mol_i]

                            spectrum_dictionary[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule][index],
                                                                                        pres=pressure, temp=temperature,
                                                                                        path_l=pathlength,
                                                                                        wl_min=self.wlmin, wl_max=self.wlmax,
                                                                                        step=0.001)
                            dict_w[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[0]
                            dict_t[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[1]

                        GMSL.storage_for_dict(dict_w, dict_t, pathlength, self.tab.slit_size, self.tab.offset_left,
                                              self.tab.offset_right)

                        model = Model(self.fitting_functions[len(list_molecules) - 1])
                        if len(list_molecules) == 1:
                            params = model.make_params(c1=1)
                        elif len(list_molecules) == 2:
                            params = model.make_params(c1=1, c2=1)
                        elif len(list_molecules) == 3:
                            params = model.make_params(c1=1, c2=1, c3=1)
                        elif len(list_molecules) == 4:
                            params = model.make_params(c1=1, c2=1, c3=1, c4=1)
                        elif len(list_molecules) == 5:
                            params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1)
                        elif len(list_molecules) == 6:
                            params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1)
                        elif len(list_molecules) == 7:
                            params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1)
                        elif len(list_molecules) == 8:
                            params = model.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1, c8=1)

                        for i in range(len(list_molecules)):
                            params['c' + str(i + 1)].set(min=-0.001, max=1000)
                        result = model.fit(t_exp_baseline_corrected, params, w=w_exp)

                        c_list_old = {}
                        for i in range(len(list_molecules)):
                            if result.best_values['c' + str(i + 1)] * dict_molecules[list_molecules[i]][index] < 0:
                                c_list_old[list_molecules[i]] = 0
                            else:
                                c_list_old[list_molecules[i]] = result.best_values['c' + str(i + 1)] * dict_molecules[
                                    list_molecules[i]][index]

                        # Re-fit if needed
                        refit_bool = True
                        for fit_i in range(5):
                            if refit_bool:
                                spectrum_dictionary_refit = {}
                                dict_w_refit = {}
                                dict_t_refit = {}
                                for mol_i in range(len(dict_molecules)):
                                    molecule = list(dict_molecules.keys())[mol_i]
                                    spectrum_dictionary_refit[molecule] = (
                                        GFL.spectrum_in_air_creator(mol=molecule,
                                                                    mf=c_list_old[molecule],
                                                                    pres=pressure, temp=temperature,
                                                                    path_l=pathlength,
                                                                    wl_min=self.wlmin,
                                                                    wl_max=self.wlmax,
                                                                    step=0.001))
                                    dict_w_refit[molecule] = \
                                    spectrum_dictionary_refit[molecule].get("transmittance_noslit")[0]
                                    dict_t_refit[molecule] = \
                                    spectrum_dictionary_refit[molecule].get("transmittance_noslit")[1]

                                GMSL.storage_for_dict(dict_w_refit, dict_t_refit, pathlength,self.tab.slit_size, self.tab.offset_left,
                                          self.tab.offset_right)

                                model_refit = Model(self.fitting_functions[len(list_molecules) - 1])
                                if len(list_molecules) == 1:
                                    params_refit = model_refit.make_params(c1=1)
                                elif len(list_molecules) == 2:
                                    params_refit = model_refit.make_params(c1=1, c2=1)
                                elif len(list_molecules) == 3:
                                    params_refit = model_refit.make_params(c1=1, c2=1, c3=1)
                                elif len(list_molecules) == 4:
                                    params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1)
                                elif len(list_molecules) == 5:
                                    params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1)
                                elif len(list_molecules) == 6:
                                    params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1)
                                elif len(list_molecules) == 7:
                                    params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1)
                                elif len(list_molecules) == 8:
                                    params_refit = model_refit.make_params(c1=1, c2=1, c3=1, c4=1, c5=1, c6=1, c7=1,
                                                                           c8=1)

                                for i in range(len(list_molecules)):
                                    params_refit['c' + str(i + 1)].set(min=0.8, max=1.2)

                                result_refit = model_refit.fit(t_exp, params_refit, w=w_exp)

                                c_list_new = {}
                                for i in range(len(list_molecules)):
                                    if result_refit.best_values['c' + str(i + 1)] * c_list_old[list_molecules[i]] < 0:
                                        c_list_new[list_molecules[i]] = 0
                                    else:
                                        c_list_new[list_molecules[i]] = result_refit.best_values['c' + str(i + 1)] * \
                                                                    c_list_old[
                                                                        list_molecules[i]]
                                correctness_count = 0
                                for i in range(len(list_molecules)):
                                    if (np.abs(c_list_new[list_molecules[i]] - c_list_old[list_molecules[i]]) /
                                        c_list_old[
                                            list_molecules[i]] > 0.005 or
                                        np.abs(c_list_new[list_molecules[i]] - c_list_old[
                                            list_molecules[i]])) < 1 * 10 ** -6:
                                        correctness_count += 1

                                    c_list_old[list_molecules[i]] = c_list_new[list_molecules[i]]
                                if correctness_count == 0:
                                    refit_bool = False

                        spectrum_dictionary_final = {}
                        molecule_final = list(dict_molecules.keys())[0]
                        dict_w_final = {}
                        dict_t_final = {}
                        dict_a_final = {}
                        a_full = []
                        for mol_i in range(len(dict_molecules)):
                            molecule = list(dict_molecules.keys())[mol_i]
                            spectrum_dictionary_final[molecule] = GFL.spectrum_in_air_creator(
                                mol=molecule,
                                mf=c_list_new[molecule],
                                pres=pressure,
                                temp=temperature,
                                path_l=pathlength,
                                wl_min=self.wlmin, wl_max=self.wlmax,
                                step=0.001)
                            dict_w_final[molecule] = spectrum_dictionary_final[molecule].get("transmittance_noslit")[0]
                            dict_t_final[molecule] = spectrum_dictionary_final[molecule].get("transmittance_noslit")[1]
                            dict_a_final[molecule] = GFL.tr_to_ab(dict_t_final[molecule])
                            if mol_i == 0:
                                a_full = dict_a_final[molecule]
                            else:
                                a_full += dict_a_final[molecule]
                        t_full = GFL.ab_to_tr(a_full)
                        spec_new2 = Spectrum({"wavenumber": dict_w_final[molecule_final], "transmittance_noslit": t_full},
                                             wunit='cm-1',
                                             units={"transmittance_noslit": ""})

                        spec_new2.apply_slit(self.tab.slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")

                        w_temp, t_temp = spec_new2.get("transmittance")

                        offset_a = (self.tab.offset_right - self.tab.offset_left) / (4000 - 1000)
                        offset_b = self.tab.offset_right - offset_a * 4000

                        list_offset_temp = offset_a * w_temp + offset_b

                        w_new_temp = np.zeros(len(w_temp))
                        for i in range(len(w_temp)):
                            w_new_temp[i] = w_temp[i] + list_offset_temp[i]
                        spec_new2 = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                                             units={"transmittance": ""}, conditions={"path_length": pathlength})
                        spec_new2.resample(w_exp, inplace=True)
                        w_new2, t_new2 = spec_new2.get("transmittance")
                        residual = t_exp_baseline_corrected - result_refit.best_fit
                        values_now = []
                        for mol_i in range(len(dict_molecules)):
                            molecule = list(dict_molecules.keys())[mol_i]
                            dict_molecules[molecule][index] = c_list_new[molecule]
                            values_now.append(dict_molecules[molecule][index])

                        for mol_i in range(len(list_molecules)):
                            if dict_molecules[list_molecules[mol_i]][index] < 1 * 10 ** -6:
                                dict_molecules[list_molecules[mol_i]][index] = 0
                                print("Mole fraction " + str(list_molecules[mol_i]) + ": " + str(0))
                            else:
                                dict_molecules[list_molecules[mol_i]][index] = dict_molecules[list_molecules[mol_i]][index]
                                print("Mole fraction " + str(list_molecules[mol_i]) + ": " + str(
                                    dict_molecules[list_molecules[mol_i]][index]))

                        self.signal_fitting_plot.emit(
                            ["ftir_fitting_" + current_file, w_new2, t_exp, result_refit.best_fit, residual])
                        with h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                                       "r+") as hf:
                            if hf["Temperature"][index] != temperature:
                                hf["Temperature"][index] = temperature
                            if hf["Pressure"][index] != pressure:
                                hf["Pressure"][index] = pressure
                            for mol_i in range(len(list_molecules)):
                                hf[list_molecules[mol_i] + "_concentration"][index] = dict_molecules[list_molecules[mol_i]][index]

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
    signal_ftir_simulation_plot = Signal(list)

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

    def start_simulation(self, s):
        try:
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
                self.parent.worker_plotting.single_transmission_plot(["Simulated Data", self.w_full, self.t_full])
                self.tab.layout.label_simulation.setText("Simulation done")
            else:
                self.parent.worker_plotting.update_single_transmission_plot(
                    ["Simulated Data", self.w_full, self.t_full])
                self.tab.layout.label_simulation.setText("Simulation done")

        except:
            self.bool_remake = False
            self.parent.worker_plotting.update_single_transmission_plot(
                ["Simulated Data", range(100), np.zeros(100)])
            self.tab.layout.label_simulation.setText("Can't do simulation, not enough conditions given")


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
            if not os.path.exists(self.directory + "\\\\" + folder):
                os.makedirs(self.directory + "\\\\" + folder)
            self.data = data[2]
            np.savetxt(self.directory + "\\\\" + folder + "\\\\" + self.name + ".txt", self.data, delimiter='\t',
                       header="Wavelength \t\t\tExperimental Transmission \t   Fitted Transmission \t\t      "
                              "Residual")

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
            if not os.path.exists(self.directory + "\\\\" + folder):
                os.makedirs(self.directory + "\\\\" + folder)

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
            figure.savefig(self.directory + "\\\\" + folder + "\\\\" + self.name + "_" + self.state + ".png", dpi=100)

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
