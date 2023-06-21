"""
In this file, one can find the used workers.
"""

from PySide6.QtGui import QFont
from PySide6.QtCore import QObject, Signal
import pyqtgraph as pg
from radis import Spectrum
import numpy as np
import scipy.optimize as sc
import General_Functions_Library as GFL
import Gas_Mixtures_Spectra_Library as GMSL
import h5py
import os
import os.path as op
import matplotlib as mpl
from matplotlib.figure import Figure


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
        self.fitting_functions.append(GMSL.spectra_linear_wavenumber_dependent_offset_1_molecule)
        self.fitting_functions.append(GMSL.spectra_linear_wavenumber_dependent_offset_2_molecules)
        self.fitting_functions.append(GMSL.spectra_linear_wavenumber_dependent_offset_3_molecules)
        self.fitting_functions.append(GMSL.spectra_linear_wavenumber_dependent_offset_4_molecules)
        self.fitting_functions.append(GMSL.spectra_linear_wavenumber_dependent_offset_5_molecules)
        self.parent = parent
        self.dict_molecules = {}

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
            bool_pathlength = False

            if len(list_molecules) != 0:
                bool_molecules = True

            list_pressure = self.tab.list_ftir_pressure
            if len(list_pressure) == 1:
                pressure = list_pressure[0]
                list_pressure = np.ones(len_files_in_dir) * pressure
                bool_pressure = True
            elif len(list_pressure) == len(list_molecules) and len(list_molecules) != 0 \
                    and list_pressure.any() != 0:
                bool_pressure = True

            list_temperature = self.tab.list_ftir_temperature
            if len(list_temperature) == 1:
                temperature = list_temperature[0]
                list_temperature = np.ones(len_files_in_dir) * temperature
                bool_temperature = True
            elif len(list_temperature) == len(list_molecules) and len(list_molecules) != 0 \
                    and list_temperature.any() != 0:
                bool_temperature = True

            list_pathlength = self.tab.list_ftir_pathlength

            if len(list_pathlength) == 1:
                pathlength = list_pathlength[0]
                list_pathlength = np.ones(len_files_in_dir) * pathlength
                bool_pathlength = True
            elif len(list_pathlength) == len(list_molecules) and len(list_molecules) != 0 \
                    and list_pathlength.any() != 0:
                bool_pathlength = True

            if bool_molecules and bool_pathlength and bool_temperature and bool_pressure:
                if not op.exists(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5"):
                    hf = h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                                   "w")
                self.tab.layout.label_info_fit.setText("Fitting for "
                                                       + current_file
                                                       + " ....")

                fit_params_p0 = []
                fit_params_bound_min = []
                fit_params_bound_max = []
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
                print(list_time)
                list_offset_min = np.zeros(len_files_in_dir)
                list_offset_max = np.zeros(len_files_in_dir)
                list_slit = np.zeros(len_files_in_dir)

                print(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5")
                with h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5", "r+") as hf:
                    keys = list(hf.keys())
                    if "Time" not in keys:
                        hf.create_dataset("Time", data=list_time)
                    elif not len(hf["Time"][:]) == len_files_in_dir:
                        del hf["Time"]
                        hf.create_dataset("Time", data=list_time)

                    if "Offset_minimum" not in keys:
                        hf.create_dataset("Offset_minimum", data=list_offset_min)
                    elif not len(hf["Offset_minimum"][:]) == len_files_in_dir:
                        del hf["Offset_minimum"]
                        hf.create_dataset("Offset_minimum", data=list_offset_min)

                    if "Offset_maximum" not in keys:
                        hf.create_dataset("Offset_maximum", data=list_offset_max)
                    elif not len(hf["Offset_maximum"][:]) == len_files_in_dir:
                        del hf["Offset_maximum"]
                        hf.create_dataset("Offset_maximum", data=list_offset_max)

                    if "Slit_size" not in keys:
                        hf.create_dataset("Slit_size", data=list_slit)
                    elif not len(hf["Slit_size"][:]) == len_files_in_dir:
                        del hf["Slit_size"]
                        hf.create_dataset("Slit_size", data=list_slit)

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

                for i in range(len(list_molecules)):
                    fit_params_p0.append(1)
                    fit_params_bound_min.append(0)
                    fit_params_bound_max.append(2000)

                fit_params_p0.append(0)
                fit_params_p0.append(0.4)
                fit_params_p0.append(0.3)

                fit_params_bound_min.append(-0.4)
                fit_params_bound_min.append(-0.2)
                fit_params_bound_min.append(0.01)

                fit_params_bound_max.append(0.4)
                fit_params_bound_max.append(1)
                fit_params_bound_max.append(2)

                fit_params_tot = len(fit_params_p0)
                try:
                    data_wavelength = np.array(data.get_range("AB")[0:-1])
                    data_transmission = np.array(data["AB"][0:-1])
                except:
                    data_wavelength = data[0]
                    data_transmission = data[1]


                s_exp = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                 units={"transmittance": ""})
                w_exp, t_exp = s_exp.get("transmittance")

                wlmin = min(data_wavelength)
                wlmax = max(data_wavelength)
                dict_molecules = {}
                for mol_i in range(len(list_molecules)):
                    dict_molecules[list_molecules[mol_i]] = 0.000010

                spectrum_dictionary = {}
                dict_w = {}
                dict_t = {}
                for mol_i in range(len(dict_molecules)):
                    molecule = list(dict_molecules.keys())[mol_i]

                    spectrum_dictionary[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule],
                                                                                pres=pressure, temp=temperature,
                                                                                path_l=pathlength,
                                                                                wl_min=wlmin, wl_max=wlmax, step=0.002)
                    dict_w[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[0]
                    dict_t[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[1]

                GMSL.storage_for_dict(dict_w, dict_t, pathlength)

                popt, pcov = sc.curve_fit(f=self.fitting_functions[len(list_molecules) - 1], xdata=w_exp,
                                          ydata=t_exp,
                                          p0=fit_params_p0,
                                          bounds=(fit_params_bound_min, fit_params_bound_max))
                spectrum_dictionary_final = {}
                molecule_final = list(dict_molecules.keys())[0]
                dict_w_final = {}
                dict_t_final = {}
                dict_a_final = {}
                a_full = []
                for mol_i in range(len(dict_molecules)):
                    molecule = list(dict_molecules.keys())[mol_i]
                    spectrum_dictionary_final[molecule] = GFL.spectrum_in_air_creator(mol=molecule,
                                                                                      mf=dict_molecules[molecule] * popt[
                                                                                          mol_i],
                                                                                      pres=pressure, temp=temperature,
                                                                                      path_l=pathlength,
                                                                                      wl_min=wlmin, wl_max=wlmax,
                                                                                      step=0.002)
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

                spec_new2.apply_slit(popt[-1], unit="cm-1", norm_by="area", inplace=True)

                w_new2_temp, t_new2_temp = spec_new2.get("transmittance")
                w_new_temp = np.zeros(len(w_new2_temp))

                offset_min = popt[-3]
                offset_max = popt[-2]

                offset_max_temp = offset_max
                if offset_max < offset_min:
                    offset_max = offset_min
                    offset_min = offset_max_temp
                offset_list_lin = np.linspace(offset_min, offset_max, len(w_new2_temp))
                for i in range(len(w_new2_temp)):
                    w_new_temp[i] = w_new2_temp[i] + offset_list_lin[i]
                spec_new2 = Spectrum({"wavenumber": w_new_temp, "transmittance": t_new2_temp}, wunit='cm-1',
                                     units={"transmittance": ""}, conditions={"path_length": 20.5})
                spec_new2.resample(w_exp, inplace=True)
                w_new2, t_new2 = spec_new2.get("transmittance")
                residual = t_exp - t_new2
                for mol_i in range(len(dict_molecules)):
                    molecule = list(dict_molecules.keys())[mol_i]
                    dict_molecules[molecule] = dict_molecules[molecule] * popt[mol_i]

                if all(i < 1 * 10 ** -6 for i in list(dict_molecules.values())):
                    print("Offset_min: " + str(0))
                    print("Offset_max: " + str(0))
                    print("Slit Size: " + str(0))
                    list_offset_min[index] = 0
                    list_offset_max[index] = 0
                    list_slit[index] = 0
                else:
                    print("Offset_min: " + str(popt[-3]))
                    print("Offset_max: " + str(popt[-2]))
                    print("Slit Size: " + str(popt[-1]))
                    list_offset_min[index] = offset_min
                    list_offset_max[index] = offset_max
                    list_slit[index] = popt[-1]

                for mol_i in range(len(list_molecules)):
                    if dict_molecules[list_molecules[mol_i]] < 1 * 10 ** -6:
                        dict_molecules[list_molecules[mol_i]] = 0
                        print("Mole fraction " + str(list_molecules[mol_i]) + ": " + str(0))
                    else:
                        dict_molecules[list_molecules[mol_i]] = dict_molecules[list_molecules[mol_i]]
                        print(
                            "Mole fraction " + str(list_molecules[mol_i]) + ": " + str(dict_molecules[list_molecules[mol_i]]))

                self.signal_fitting_plot.emit(
                    ["ftir_fitting_" + current_file, w_new2, t_exp, t_new2, residual])
                with h5py.File(
                        self.parent.tab_ftir_fitting.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                        "r+") as hf:
                    hf["Offset_minimum"][index] = list_offset_min[index]
                    hf["Offset_maximum"][index] = list_offset_max[index]
                    hf["Slit_size"][index] = list_slit[index]
                    for mol_i in range(len(list_molecules)):
                        hf[list_molecules[mol_i] + "_concentration"][index] = dict_molecules[list_molecules[mol_i]]

                self.tab.layout.label_info_fit.setText("Fitting for " + current_file + " made")
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

        time_list = []
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
                    time_list.append(time_spend)
                else:
                    time_spend = 60 * 60 * (int(hour) - hour_0) + 60 * (int(minute) - minute_0) + (
                            int(second) - second_0)
                    time_list.append(time_spend)
            except:
                time_list.append(0)
        offset_min_list = np.zeros(len_files_in_dir)
        offset_max_list = np.zeros(len_files_in_dir)
        slit_list = np.zeros(len_files_in_dir)

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
        elif len(list_pressure) == len(list_molecules) and len(list_molecules) != 0 \
                and list_pressure.any() != 0:
            bool_pressure = True

        list_temperature = self.tab.list_ftir_temperature
        if len(list_temperature) == 1:
            temperature = list_temperature[0]
            list_temperature = np.ones(len_files_in_dir) * temperature
            bool_temperature = True
        elif len(list_temperature) == len(list_molecules) and len(list_molecules) != 0 \
                and list_temperature.any() != 0:
            bool_temperature = True

        list_pathlength = self.tab.list_ftir_pathlength

        if len(list_pathlength) == 1:
            pathlength = list_pathlength[0]
            list_pathlength = np.ones(len_files_in_dir) * pathlength
            bool_pathlength = True
        elif len(list_pathlength) == len(list_molecules) and len(list_molecules) != 0 \
                and list_pathlength.any() != 0:
            bool_pathlength = True

        for index in range(len_files_in_dir):
            current_file = files_in_directory[index]
            data_text_file_selected = self.tab.directory_save_invenioR_processed + \
                                      "\\\\Data (in text files)" + "\\\\" + current_file + ".txt"
            data = self.inner_tab.array_data[index]

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

                    fit_params_p0 = []
                    fit_params_bound_min = []
                    fit_params_bound_max = []

                    with h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                                   "r+") as hf:
                        keys = list(hf.keys())
                        if "Time" not in keys:
                            hf.create_dataset("Time", data=time_list)
                        elif not len(hf["Time"][:]) == len_files_in_dir:
                            del hf["Time"]
                            hf.create_dataset("Time", data=time_list)

                        if "Offset_minimum" not in keys:
                            hf.create_dataset("Offset_minimum", data=offset_min_list)
                        elif not len(hf["Offset_minimum"][:]) == len_files_in_dir:
                            del hf["Offset_minimum"]
                            hf.create_dataset("Offset_minimum", data=offset_min_list)

                        if "Offset_maximum" not in keys:
                            hf.create_dataset("Offset_maximum", data=offset_max_list)
                        elif not len(hf["Offset_maximum"][:]) == len_files_in_dir:
                            del hf["Offset_maximum"]
                            hf.create_dataset("Offset_maximum", data=offset_max_list)

                        if "Slit_size" not in keys:
                            hf.create_dataset("Slit_size", data=slit_list)
                        elif not len(hf["Slit_size"][:]) == len_files_in_dir:
                            del hf["Slit_size"]
                            hf.create_dataset("Slit_size", data=slit_list)

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
                            try:
                                print(dict_molecules[list_molecules[mol_i]])
                            except:
                                print("Nope")

                            if list_molecules[mol_i] not in list(dict_molecules.keys()):
                                dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                                print(dict_molecules)
                            elif not len(dict_molecules[list_molecules[mol_i]]) == len_files_in_dir:
                                dict_molecules[list_molecules[mol_i]] = np.zeros(len(files_in_directory))
                            if list_molecules[mol_i] + "_concentration" not in keys:
                                hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                  data=np.zeros(len(files_in_directory)))
                            elif not len(hf[list_molecules[mol_i] + "_concentration"]) == len_files_in_dir:
                                del hf[list_molecules[mol_i] + "_concentration"]
                                hf.create_dataset(list_molecules[mol_i] + "_concentration",
                                                  data=np.zeros(len(files_in_directory)))

                    for i in range(len(list_molecules)):
                        fit_params_p0.append(1)
                        fit_params_bound_min.append(0)
                        fit_params_bound_max.append(20000)

                    fit_params_p0.append(0.1)
                    fit_params_p0.append(0.4)
                    fit_params_p0.append(0.3)

                    fit_params_bound_min.append(-0.4)
                    fit_params_bound_min.append(-0.2)
                    fit_params_bound_min.append(0.01)

                    fit_params_bound_max.append(0.4)
                    fit_params_bound_max.append(1)
                    fit_params_bound_max.append(1)

                    fit_params_tot = len(fit_params_p0)
                    try:
                        data_wavelength = np.array(data.get_range("AB")[0:-1])
                        data_transmission = np.array(data["AB"][0:-1])
                    except:
                        data_wavelength = data[0]
                        data_transmission = data[1]

                    s_exp = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                     units={"transmittance": ""})
                    w_exp, t_exp = s_exp.get("transmittance")

                    wlmin = min(data_wavelength)
                    wlmax = max(data_wavelength)
                    for mol_i in range(len(list_molecules)):
                        dict_molecules[list_molecules[mol_i]][index] = 0.000010
                    spectrum_dictionary = {}
                    dict_w = {}
                    dict_t = {}
                    for mol_i in range(len(dict_molecules)):
                        molecule = list(dict_molecules.keys())[mol_i]

                        spectrum_dictionary[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule][index],
                                                                                    pres=pressure, temp=temperature,
                                                                                    path_l=pathlength,
                                                                                    wl_min=wlmin, wl_max=wlmax,
                                                                                    step=0.002)
                        dict_w[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[0]
                        dict_t[molecule] = spectrum_dictionary[molecule].get("transmittance_noslit")[1]

                    GMSL.storage_for_dict(dict_w, dict_t, pathlength)

                    popt, pcov = sc.curve_fit(f=self.fitting_functions[len(list_molecules) - 1], xdata=w_exp,
                                              ydata=t_exp,
                                              p0=fit_params_p0,
                                              bounds=(fit_params_bound_min, fit_params_bound_max))
                    spectrum_dictionary_final = {}
                    molecule_final = list(dict_molecules.keys())[0]
                    dict_w_final = {}
                    dict_t_final = {}
                    dict_a_final = {}
                    a_full = []
                    for mol_i in range(len(dict_molecules)):
                        molecule = list(dict_molecules.keys())[mol_i]
                        spectrum_dictionary_final[molecule] = GFL.spectrum_in_air_creator(mol=molecule,
                                                                                          mf=dict_molecules[molecule][index] * popt[
                                                                                              mol_i],
                                                                                          pres=pressure,
                                                                                          temp=temperature,
                                                                                          path_l=pathlength,
                                                                                          wl_min=wlmin, wl_max=wlmax,
                                                                                          step=0.002)
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

                    spec_new2.apply_slit(popt[-1], unit="cm-1", norm_by="area", inplace=True)

                    w_new2_temp, t_new2_temp = spec_new2.get("transmittance")
                    w_new_temp = np.zeros(len(w_new2_temp))
                    offset_list_lin = np.linspace(popt[-3], popt[-2], len(w_new2_temp))
                    for i in range(len(w_new2_temp)):
                        w_new_temp[i] = w_new2_temp[i] + offset_list_lin[i]
                    spec_new2 = Spectrum({"wavenumber": w_new_temp, "transmittance": t_new2_temp}, wunit='cm-1',
                                         units={"transmittance": ""}, conditions={"path_length": 20.5})
                    spec_new2.resample(w_exp, inplace=True)
                    w_new2, t_new2 = spec_new2.get("transmittance")
                    residual = t_exp - t_new2
                    values_now = []
                    for mol_i in range(len(dict_molecules)):
                        molecule = list(dict_molecules.keys())[mol_i]
                        dict_molecules[molecule][index] = dict_molecules[molecule][index] * popt[mol_i]
                        values_now.append(dict_molecules[molecule][index] * popt[mol_i])

                    if all(i < 1 * 10 ** -6 for i in values_now):
                        print("Offset_min: " + str(0))
                        print("Offset_max: " + str(0))
                        print("Slit Size: " + str(0))
                        offset_min_list[index] = 0
                        offset_max_list[index] = 0
                        slit_list[index] = 0
                    else:
                        print("Offset_min: " + str(popt[-3]))
                        print("Offset_max: " + str(popt[-2]))
                        print("Slit Size: " + str(popt[-1]))
                        offset_min_list[index] = popt[-3]
                        offset_max_list[index] = popt[-2]
                        slit_list[index] = popt[-1]

                    for mol_i in range(len(list_molecules)):
                        if dict_molecules[list_molecules[mol_i]][index] < 1 * 10 ** -6:
                            dict_molecules[list_molecules[mol_i]][index] = 0
                            print("Mole fraction " + str(list_molecules[mol_i]) + ": " + str(0))
                        else:
                            dict_molecules[list_molecules[mol_i]][index] = dict_molecules[list_molecules[mol_i]][index]
                            print("Mole fraction " + str(list_molecules[mol_i]) + ": " + str(
                                dict_molecules[list_molecules[mol_i]][index]))

                    self.signal_fitting_plot.emit(
                        ["ftir_fitting_" + current_file, w_new2, t_exp, t_new2, residual])
                    with h5py.File(self.tab.directory_save_invenioR_processed + "\\" + "Measurement_file.h5",
                                   "r+") as hf:
                        hf["Offset_minimum"][index] = offset_min_list[index]
                        hf["Offset_maximum"][index] = offset_max_list[index]
                        hf["Slit_size"][index] = slit_list[index]
                        for mol_i in range(len(list_molecules)):
                            hf[list_molecules[mol_i] + "_concentration"][index] = dict_molecules[list_molecules[mol_i]][index]
                            print(hf[list_molecules[mol_i] + "_concentration"][index])

                    self.tab.layout.label_info_fit.setText("Fitting for " + current_file + " made")
                    all_fits_made_bool = True
                else:
                    self.tab.layout.label_info_fit.setText("Can't fit for "
                                                           + current_file + ", ""incorrect/no molecules, pressure, "
                                                                            "temperature or pathlength")
                    all_fits_made_bool = False
                    break

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
                        or self.wlmin_new != self.wlmin_old or self.wlmax_new != self.wlmax_old:
                    self.dict_spectrum[mol_i] = GFL.spectrum_in_air_creator(mol=mol_i,
                                                                                  mf=self.dict_molecules_new[mol_i],
                                                                                  pres=self.pressure_new,
                                                                                  temp=self.temperature_new,
                                                                                  path_l=self.pathlength_new,
                                                                                  wl_min=self.wlmin_new,
                                                                                  wl_max=self.wlmax_new, step=0.002)

                    self.dict_w_new[mol_i] = self.dict_spectrum[mol_i].get("transmittance_noslit")[0]
                    self.dict_t_new[mol_i] = self.dict_spectrum[mol_i].get("transmittance_noslit")[1]
                    self.dict_a_new[mol_i] = GFL.tr_to_ab(self.dict_t_new[mol_i])

                    self.tab.layout.label_simulation.setText(
                        "Simulated data for: " + str(self.list_simulated_molecules))

                if mol_i in self.dict_molecules_old.keys() and self.temperature_new == self.temperature_old \
                        and self.pressure_new == self.pressure_old and self.pathlength_new == self.pathlength_old \
                        and self.wlmin_new == self.wlmin_old and self.wlmax_new == self.wlmax_old:
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
