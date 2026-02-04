"""
File containing all of the code to create the Tabs shown within the GUI
"""

import Code_Gui.Gui_General_Code.General_Functions_Library as GFL
from PySide6.QtCore import QObject, Signal, Qt
from PySide6.QtWidgets import QWidget, QLabel, QTabWidget, QVBoxLayout, QHBoxLayout, QGridLayout, \
    QLineEdit, QPushButton, QCheckBox, QSizePolicy
import os
import numpy as np
import os.path as op
import radis as rd
import sys
import platform

class create_tab_ftir_simulator(QWidget):
    signal_ftir_simulation_start = Signal(list)
    signal_ftir_simulation_set_parent = Signal(list)
    signal_ftir_simulation_tr_vs_ab = Signal(list)
    signal_ftir_simulation_save_data_in_txt = Signal(list)
    signal_ftir_simulation_save_data_as_png = Signal(list)

    def __init__(self):
        super().__init__()
        self.absorbance_bool = False
        self.molecule_list = ["H2O", "CO2", "CO", "CH4", "NO", "NO2", "N2O", "HNO3", "C2H2", "C2H4", "O3", "NH3", "OH",
                              "C2H6", "CH3OH", "C4H2", "H2CO", "HCN", "H2O2", "O", "HCOOH", "HO2", "C2N2", "CH3OH"]
        self.len_table = len(self.molecule_list) / 6
        self.rows = len(self.molecule_list) / 6
        self.molecule_storage_for_simulation = {}
        self.molecule_storage_for_errors = []
        self.molecule_storage_for_simulation_final = {}
        self.temperature_for_simulation = 0
        self.pressure_for_simulation = 0
        self.pathlength_for_simulation = 0

        self.directory_for_saving = ""
        self.name_for_saving = ""

    def create_tab(self, parent):
        self.parent = parent
        self.layout = QGridLayout(self)

        self.layout.label_directory = QLabel("Enter directory for saving data here (not essential): ")
        self.layout.label_directory.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_directory = QLineEdit()
        self.layout.line_edit_directory.textChanged.connect(lambda: self.text_changed("Directory"))
        self.layout.line_edit_directory.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_directory = QPushButton("Save directory")
        self.layout.button_directory.clicked.connect(self.button_pressed)
        self.layout.label_directory2 = QLabel("No directory given yet")
        self.layout.label_directory2.setAlignment(Qt.AlignCenter)


        self.layout.addWidget(self.layout.label_directory, 0,0,1,2)
        self.layout.addWidget(self.layout.line_edit_directory, 0,2,1,3)
        self.layout.addWidget(self.layout.button_directory, 0,5,1,1)
        self.layout.addWidget(self.layout.label_directory2, 0,6,1,2)

        self.layout.label_temperature = QLabel("Enter gas temperature here (in Kelvin): ")
        self.layout.label_temperature.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_temperature = QLineEdit("296.35")
        self.layout.line_edit_temperature.textChanged.connect(lambda: self.text_changed("Temperature"))
        self.layout.line_edit_temperature.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_save_temperature = QPushButton("Save temperature")
        self.layout.button_save_temperature.clicked.connect(self.button_pressed)
        self.layout.label_temperature2 = QLabel("No temperature saved yet")
        self.layout.label_temperature2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_temperature, 1,0,1,2)
        self.layout.addWidget(self.layout.line_edit_temperature, 1,2,1,3)
        self.layout.addWidget(self.layout.button_save_temperature, 1,5, 1, 1)
        self.layout.addWidget(self.layout.label_temperature2, 1, 6, 1,2)

        self.layout.label_pressure = QLabel("Enter gas pressure here (in bar): ")
        self.layout.label_pressure.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_pressure = QLineEdit("1.005")
        self.layout.line_edit_pressure.textChanged.connect(
            lambda: self.text_changed("Pressure"))
        self.layout.line_edit_pressure.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_save_pressure = QPushButton("Save pressure")
        self.layout.button_save_pressure.clicked.connect(self.button_pressed)
        self.layout.label_pressure2 = QLabel("No pressure saved yet")
        self.layout.label_pressure2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_pressure, 2, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_pressure, 2, 2, 1, 3)
        self.layout.addWidget(self.layout.button_save_pressure, 2, 5, 1, 1)
        self.layout.addWidget(self.layout.label_pressure2, 2, 6, 1, 2)

        self.layout.label_pathlength = QLabel("Enter pathlength of your gascel here (in cm): ")
        self.layout.label_pathlength.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_pathlength = QLineEdit("20.062")
        self.layout.line_edit_pathlength.textChanged.connect(lambda: self.text_changed("Pathlength"))
        self.layout.line_edit_pathlength.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_save_pathlength = QPushButton("Save pathlength")
        self.layout.button_save_pathlength.clicked.connect(self.button_pressed)
        self.layout.label_pathlength2 = QLabel("No pathlength saved yet")
        self.layout.label_pathlength2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_pathlength, 3, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_pathlength, 3, 2, 1,3)
        self.layout.addWidget(self.layout.button_save_pathlength, 3, 5, 1,1)
        self.layout.addWidget(self.layout.label_pathlength2, 3, 6, 1, 2)

        self.layout.label_wavenumber_min = QLabel("Enter minimum wavenumber(in cm^-1)")
        self.layout.label_wavenumber_min.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_wavenumber_min = QLineEdit("400")
        self.layout.line_edit_wavenumber_min.textChanged.connect(lambda: self.text_changed("Minimum wavenumber"))
        self.layout.line_edit_wavenumber_min.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.label_wavenumber_max = QLabel(" and maximum wavenumber(in cm^-1)")
        self.layout.label_wavenumber_max.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_wavenumber_max = QLineEdit("4000")
        self.layout.line_edit_wavenumber_max.textChanged.connect(lambda: self.text_changed("Maximum wavenumber"))
        self.layout.line_edit_wavenumber_max.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_wavenumber = QPushButton("Save wavenumbers")
        self.layout.button_wavenumber.clicked.connect(self.button_pressed)
        self.layout.label_wavenumber = QLabel("No wavenumbers saved yet")
        self.layout.label_wavenumber.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_wavenumber_min, 4, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_wavenumber_min, 4, 2, 1, 1)
        self.layout.addWidget(self.layout.label_wavenumber_max, 4, 3, 1, 1)
        self.layout.addWidget(self.layout.line_edit_wavenumber_max, 4, 4, 1, 1)
        self.layout.addWidget(self.layout.button_wavenumber, 4, 5, 1, 1)
        self.layout.addWidget(self.layout.label_wavenumber, 4, 6, 1, 2)

        self.layout.label_slitsize = QLabel("Extra Parameters // \t\t\t\t Slit size: ")
        self.layout.label_slitsize.setStyleSheet("font: bold ; font-size:11pt")
        self.layout.line_edit_slitsize = QLineEdit("")
        self.layout.line_edit_slitsize.textChanged.connect(lambda: self.text_changed("Slit size"))
        self.layout.line_edit_slitsize.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_slitsize = QPushButton("Save slit size")
        self.layout.button_slitsize.clicked.connect(self.button_pressed)
        self.layout.label_slitsize_saved = QLabel("No slit size saved yet")
        self.layout.label_slitsize_saved.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_slitsize, 5, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_slitsize, 5, 2, 1, 3)
        self.layout.addWidget(self.layout.button_slitsize, 5, 5, 1, 1)
        self.layout.addWidget(self.layout.label_slitsize_saved, 5, 6, 1, 2)


        self.layout.label_molecule = \
            QLabel("Choose molecules for simulation (give their mole fraction (1000 ppm = 0.001 mole fraction))")
        self.layout.label_molecule.setStyleSheet("font-weight: bold ; font-size:14pt")
        self.layout.label_molecule.setAlignment(Qt.AlignCenter)
        self.layout.label_molecule.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.addWidget(self.layout.label_molecule, 6,0,1,8)

        # Create Molecule Grid
        self.layout.grid_molecule = QGridLayout()
        self.dict_molecule_labels={}
        self.dict_molecule_text={}
        for mol_i in range(len(self.molecule_list)):
            temp_layout = QHBoxLayout()
            column = int(mol_i/self.rows)
            mol = self.molecule_list[mol_i]
            mol_label_temp = QLabel(mol)
            text_space_temp = QLineEdit("0")
            self.dict_molecule_labels[mol] = mol_label_temp
            self.dict_molecule_text[mol] = text_space_temp


            temp_layout.addWidget(self.dict_molecule_labels[mol])
            temp_layout.addWidget(self.dict_molecule_text[mol])
            self.dict_molecule_text[mol].textChanged.connect(self.molecules_concentration_chosen)
            self.layout.grid_molecule.addLayout(temp_layout, mol_i % 4, column)
        self.layout.addLayout(self.layout.grid_molecule,7,0,2,8)

        self.layout.label_molecules = QLabel("No molecules chosen/No concentrations given yet")
        self.layout.label_molecules.setAlignment(Qt.AlignCenter)
        self.layout.button_save_molecules = QPushButton("Save chosen molecules and concentrations")
        self.layout.button_save_molecules.clicked.connect(self.button_pressed)

        self.layout.addWidget(self.layout.button_save_molecules,9,0,1,1)
        self.layout.addWidget(self.layout.label_molecules, 9,1,1,7)

        self.layout.button_start_simulation = QPushButton("Start simulation")
        self.layout.button_start_simulation.clicked.connect(self.button_pressed)

        self.signal_ftir_simulation_start.connect(self.parent.worker_ftir_simulation.start_simulation)
        self.layout.addWidget(self.layout.button_start_simulation, 10,0,1,4)

        self.layout.label_simulation = QLabel("No simulation done yet")
        self.layout.addWidget(self.layout.label_simulation,10,4,1,2)

        self.layout.button_transmission_vs_absorbance = QPushButton("Showing transmission")
        self.layout.button_transmission_vs_absorbance.clicked.connect(self.button_pressed)
        self.signal_ftir_simulation_tr_vs_ab.connect(self.parent.worker_plotting.update_single_transmission_plot)
        self.layout.addWidget(self.layout.button_transmission_vs_absorbance,10,6,1,2)

        plot_temp = self.parent.worker_plotting.create_empty_plot("Simulated Data")

        self.layout.addWidget(plot_temp, 11,0,12,8)

        self.layout.label_name = QLabel("Enter name for data files: ")
        self.layout.label_name.setStyleSheet("font-weight: bold ; font-size:11pt")
        self.layout.line_edit_name = QLineEdit()
        self.layout.line_edit_name.textChanged.connect(lambda: self.text_changed("Name"))

        self.layout.button_save_simulated_data_in_txt = QPushButton("Save simulated data in text file")
        self.layout.button_save_simulated_data_in_txt.clicked.connect(self.button_pressed)
        self.signal_ftir_simulation_save_data_in_txt.connect(self.parent.worker_saving.save_data_in_txt)
        self.layout.button_save_simulated_data_as_png = QPushButton("Save simulated data as png")
        self.layout.button_save_simulated_data_as_png.clicked.connect(self.button_pressed)
        self.signal_ftir_simulation_save_data_as_png.connect(self.parent.worker_saving.save_data_as_png)
        self.layout.label_save = QLabel("No data saved yet")
        self.layout.label_save.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_name, 23, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_name, 23, 2, 1, 2)
        self.layout.addWidget(self.layout.button_save_simulated_data_in_txt, 23, 4, 1, 1)
        self.layout.addWidget(self.layout.button_save_simulated_data_as_png, 23, 5, 1, 1)
        self.layout.addWidget(self.layout.label_save, 23, 6, 1, 2)

        self.parent.tab_ftir_simulator.setLayout(self.layout)

    def molecules_concentration_chosen(self,s):
        for mol_i in self.molecule_list:
            if self.sender() == self.dict_molecule_text[mol_i]:
                try:
                    if float(s) < 1 and float(s) != 0:
                        if not mol_i in self.molecule_storage_for_simulation.keys():
                            self.molecule_storage_for_simulation[mol_i] = float(s)
                        else:
                            del self.molecule_storage_for_simulation[mol_i]
                            self.molecule_storage_for_simulation[mol_i] = float(s)
                        self.layout.label_molecules.setText(
                            mol_i + " concentration has been updated")

                        if mol_i in self.molecule_storage_for_errors:
                            index = self.molecule_storage_for_errors.index(mol_i)
                            del self.molecule_storage_for_errors[index]
                    elif float(s) == 0:
                        if mol_i in self.molecule_storage_for_simulation.keys():
                            del self.molecule_storage_for_simulation[mol_i]

                    else:
                        if not mol_i in self.molecule_storage_for_errors:
                            self.molecule_storage_for_errors.append(mol_i)
                            self.layout.label_molecules.setText(
                                mol_i + " concentration couldn't be updated, incorrect concentration given")

                except:
                    if s == "":
                        if mol_i in self.molecule_storage_for_simulation.keys():
                            del self.molecule_storage_for_simulation[mol_i]

                        self.layout.label_molecules.setText(
                            mol_i + " concentration has been updated")

                        if mol_i in self.molecule_storage_for_errors:
                            index = self.molecule_storage_for_errors.index(mol_i)
                            del self.molecule_storage_for_errors[index]

                    else:
                        if not mol_i in self.molecule_storage_for_errors:
                            self.molecule_storage_for_errors.append(mol_i)
                        self.layout.label_molecules.setText(
                            mol_i + " concentration couldn't be updated, incorrect concentration given")

    def text_changed(self, s):
        if s == "Temperature":
            self.layout.label_temperature2.setText("Temperature changed, not yet stored")
        elif s == "Pressure":
            self.layout.label_pressure2.setText("Pressure changed, not yet stored")
        elif s == "Pathlength":
            self.layout.label_pathlength2.setText("Pathlength changed, not yet stored")
        elif s== "Minimum wavenumber":
            self.layout.label_wavenumber.setText("Minimum wavenumber changed, no wavenumbers stored yet")
        elif s=="Maximum wavenumber":
            self.layout.label_wavenumber.setText("Maximum wavenumber changed, no wavenumbers stored yet")
        elif s =="Directory":
            self.layout.label_directory2.setText("Directory changed, not stored yet")
        elif s == "Name":
            self.name_for_saving = self.layout.line_edit_name.text()
        elif s=="Slit size":
            self.layout.label_slitsize_saved.setText("Slit size changed, not yet stored")

    def button_pressed(self,s):
        if self.sender().text() == "Save directory":
            self.directory_for_saving = self.layout.line_edit_directory.text().replace("\\", "\\\\") + "\\\\"
            print(self.directory_for_saving)
            if self.directory_for_saving == "":
                self.layout.label_directory2.setText("No directory given.")
            elif platform.system() == "Windows":
                if ("C:\\" not in self.directory_for_saving and "P:\\" not in self.directory_for_saving and
                      "D:\\" not in self.directory_for_saving and "F:\\" not in self.directory_for_saving and
                      "K:\\" not in self.directory_for_saving and "L:\\" not in self.directory_for_saving):
                    self.layout.label_directory2.setText("Incorrect directory given.")
                else:
                    if not os.path.exists(self.directory_for_saving):
                        os.makedirs(self.directory_for_saving)
                        self.layout.label_directory2.setText("Directory exists.")

                    if os.path.exists(self.directory_for_saving):
                        self.layout.label_directory2.setText("Directory exists.")
            elif platform.system() == "Darwin":
                if not os.path.exists(self.directory_for_saving):
                    os.makedirs(self.directory_for_saving)
                    self.layout.label_directory2.setText("Directory exists.")

                if os.path.exists(self.directory_for_saving):
                    self.layout.label_directory2.setText("Directory exists.")
        elif self.sender().text() == "Save temperature":
            try:
                self.temperature_for_simulation = float(self.layout.line_edit_temperature.text())
                self.layout.label_temperature2.setText("Temperature stored")
            except:
                self.layout.label_temperature2.setText("Uncorrect temperature given")
        elif self.sender().text() == "Save pressure":
            try:
                self.pressure_for_simulation = float(self.layout.line_edit_pressure.text())
                self.layout.label_pressure2.setText("Pressure stored")
            except:
                self.layout.label_pressure2.setText("Uncorrect pressure given")
        elif self.sender().text() == "Save pathlength":
            try:
                self.pathlength_for_simulation = float(self.layout.line_edit_pathlength.text())
                self.layout.label_pathlength2.setText("Pathlength stored")
            except:
                self.layout.label_pathlength2.setText("Uncorrect pathlength given")
        elif self.sender().text() == "Save slit size":
            try:
                self.slitsize_for_simulation = float(self.layout.line_edit_slitsize.text())
                self.layout.label_slitsize_saved.setText("Slit size stored")
            except:
                if self.layout.line_edit_slitsize.text() == "":
                    self.slitsize_for_simulation = 0
                    self.layout.label_slitsize_saved.setText("Slit size stored (no slit used)")
                else:
                    self.layout.label_slitsize_saved.setText("Uncorrect slitsize given")
        elif self.sender().text() == "Save chosen molecules and concentrations":
            if len(self.molecule_storage_for_errors) == 0:
                if len(self.molecule_storage_for_simulation) > 0:
                    self.molecule_storage_for_simulation_final = self.molecule_storage_for_simulation
                    self.layout.label_molecules.setText("All molecules and their concentrations have been stored")
                else:
                    self.layout.label_molecules.setText("No molecules chosen yet")
            else:
                self.layout.label_molecules.setText("Still some errors remain in given concentrations, please check")
        elif self.sender().text() == "Start simulation":
            self.layout.label_simulation.setText("Simulation started...")
            self.signal_ftir_simulation_start.emit([True])
        elif self.sender().text() == "Save wavenumbers":
                try:
                    self.wavenumber_min_for_simulation = float(self.layout.line_edit_wavenumber_min.text())
                    self.wavenumber_max_for_simulation = float(self.layout.line_edit_wavenumber_max.text())
                    self.layout.label_wavenumber.setText("Both wavenumbers stored")
                except:
                    self.layout.label_wavenumber.setText("Uncorrect minimum or maximum wavenumber given")
        elif self.sender().text() == "Showing transmission":
            try:
                self.absorbance_bool = True
                self.signal_ftir_simulation_tr_vs_ab.emit(["Simulated Data", self.absorbance_bool])
                self.layout.button_transmission_vs_absorbance.setText("Showing absorbance")
            except:
                self.absorbance_bool = True
                self.layout.button_transmission_vs_absorbance.setText("Showing absorbance")
        elif self.sender().text() == "Showing absorbance":
            try:
                self.absorbance_bool = False
                self.signal_ftir_simulation_tr_vs_ab.emit(["Simulated Data",self.absorbance_bool])
                self.layout.button_transmission_vs_absorbance.setText("Showing transmission")
            except:
                self.absorbance_bool = False
                self.layout.button_transmission_vs_absorbance.setText("Showing transmission")
        elif self.sender().text() == "Save simulated data in text file":
            if self.name_for_saving != "" and self.directory_for_saving != "":
                try:
                    self.signal_ftir_simulation_save_data_in_txt.emit([
                        self.directory_for_saving, self.name_for_saving,
                        np.transpose([self.parent.worker_plotting.dict_plots["Simulated Data"].p1.x,
                                      self.parent.worker_plotting.dict_plots["Simulated Data"].p1.y]),
                        "Simulated Data"])
                    self.layout.label_save.setText("Simulated data saved in text file")
                except:
                    self.layout.label_save.setText("Can't save, no simulated data yet")
            else:
                self.layout.label_save.setText("No directory or name given")
        elif self.sender().text() == "Save simulated data as png":
            if self.name_for_saving != "" and self.directory_for_saving != "":
                try:
                    self.signal_ftir_simulation_save_data_as_png.emit([
                        self.directory_for_saving, self.name_for_saving,
                        [self.parent.worker_plotting.dict_plots["Simulated Data"].p1.x,
                         self.parent.worker_plotting.dict_plots["Simulated Data"].p1.y],
                        "Simulated Data"])
                    self.layout.label_save.setText("Simulated data saved in png")
                except:
                    self.layout.label_save.setText("Can't save, no simulated data yet")
            else:
                self.layout.label_save.setText("No directory or name given")

class create_tab_ftir_fitting(QWidget):
    signal_create_inner_tab = Signal(list)
    def __init__(self):
        super().__init__()
        self.absorbance_bool = False

    def create_tab(self, parent):
        self.parent = parent
        self.molecule_list = ["H2O", "CO2", "CO", "CH4", "NO", "NO2", "N2O", "HNO3", "C2H2", "C2H4", "O3", "NH3", "OH",
                              "C2H6", "CH3OH", "C4H2", "H2CO", "HCN", "H2O2", "O", "HCOOH", "HO2", "C2N2", "CH3OH" ]
        self.len_table = len(self.molecule_list) / 6
        self.rows = len(self.molecule_list) / 6
        self.molecule_storage_for_fitting = []
        self.ftir_temperature_list = []
        self.ftir_pressure_list = []
        self.ftir_pathlength_list = []
        self.directory_data = ""
        self.directory_data_InvenioR = ""
        self.directory_save_invenioR_processed = ""
        self.new_directory_bool = False
        self.background_bool = False
        self.lnc_mct_detector_bool = False
        self.layout = QGridLayout()

        # Create row for directory choosing/creation
        self.layout.label_directory = QLabel("Choose folder with FTIR data for analysis: ")
        self.layout.label_directory.setStyleSheet("font-weight: bold ; font-size:10pt")
        self.layout.line_edit_directory = QLineEdit()
        self.layout.line_edit_directory.textChanged.connect(lambda:self.text_changed("Directory"))
        self.layout.line_edit_directory.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_directory = QPushButton("Create directory")
        self.layout.button_directory.clicked.connect(self.button_pressed)
        self.layout.label_directory2 = QLabel("Save directory Unknown")
        self.layout.label_directory2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_directory, 0, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_directory, 0, 2, 1, 5)
        self.layout.addWidget(self.layout.button_directory, 0, 7, 1, 1)
        self.layout.addWidget(self.layout.label_directory2, 0, 8, 1, 2)

        # Create row to add temperature data
        self.layout.label_temperature = QLabel("Add temperature data (either single value or list) in Kelvin: ")
        self.layout.label_temperature.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_temperature = QLineEdit()
        self.layout.line_edit_temperature.textChanged.connect(lambda:self.text_changed("Temperature"))
        self.layout.line_edit_temperature.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_temperature = QPushButton("Save temperature data")
        self.layout.button_temperature.clicked.connect(self.button_pressed)
        self.layout.label_temperature2 = QLabel("No temperature data given yet.")
        self.layout.label_temperature2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_temperature, 1, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_temperature, 1, 2, 1, 5)
        self.layout.addWidget(self.layout.button_temperature, 1, 7, 1, 1)
        self.layout.addWidget(self.layout.label_temperature2,1, 8, 1, 2)

        # Create row to add pressure data
        self.layout.label_pressure = QLabel("Add pressure data (either single value or list) in bar: ")
        self.layout.label_pressure.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_pressure = QLineEdit()
        self.layout.line_edit_pressure.textChanged.connect(lambda:self.text_changed("Pressure"))
        self.layout.line_edit_pressure.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_pressure = QPushButton("Save pressure data")
        self.layout.button_pressure.clicked.connect(self.button_pressed)
        self.layout.label_pressure2 = QLabel("No pressure data given yet.")
        self.layout.label_pressure2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_pressure, 2, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_pressure, 2, 2, 1, 5)
        self.layout.addWidget(self.layout.button_pressure, 2, 7, 1, 1)
        self.layout.addWidget(self.layout.label_pressure2, 2, 8, 1, 2)

        # Create Row for pathlength
        self.layout.label_pathlength = QLabel(
            "Add path length used during experiment (only of actual gas-cel) in cm: ")
        self.layout.label_pathlength.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_pathlength = QLineEdit("20.062")
        self.layout.line_edit_pathlength.textChanged.connect(lambda:self.text_changed("Pathlength"))
        self.layout.line_edit_pathlength.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_pathlength = QPushButton("Save path length")
        self.layout.button_pathlength.clicked.connect(self.button_pressed)
        self.layout.label_pathlength2 = QLabel("No path length given yet.")
        self.layout.label_pathlength2.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_pathlength, 3, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_pathlength, 3, 2, 1, 5)
        self.layout.addWidget(self.layout.button_pathlength, 3, 7, 1, 1)
        self.layout.addWidget(self.layout.label_pathlength2, 3, 8, 1, 2)

        # Create row for wavenumber choosing
        self.layout.label_wavenumber_min = QLabel("Enter minimum wavenumber(in cm^-1): ")
        self.layout.label_wavenumber_min.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_wavenumber_min = QLineEdit("1000")
        self.layout.line_edit_wavenumber_min.textChanged.connect(lambda: self.text_changed("Minimum wavenumber"))
        self.layout.line_edit_wavenumber_min.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.label_wavenumber_max = QLabel(" and maximum wavenumber(in cm^-1) : ")
        self.layout.label_wavenumber_max.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_wavenumber_max = QLineEdit("4150")
        self.layout.line_edit_wavenumber_max.textChanged.connect(lambda: self.text_changed("Maximum wavenumber"))
        self.layout.line_edit_wavenumber_max.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.button_wavenumber = QPushButton("Save wavenumbers")
        self.layout.button_wavenumber.clicked.connect(self.button_pressed)
        self.layout.label_wavenumber = QLabel("No wavenumbers saved yet")
        self.layout.label_wavenumber.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_wavenumber_min, 4, 0, 1, 2)
        self.layout.addWidget(self.layout.line_edit_wavenumber_min, 4, 2, 1, 2)
        self.layout.addWidget(self.layout.label_wavenumber_max, 4, 4, 1, 1)
        self.layout.addWidget(self.layout.line_edit_wavenumber_max, 4, 5, 1, 2)
        self.layout.addWidget(self.layout.button_wavenumber, 4, 7, 1, 1)
        self.layout.addWidget(self.layout.label_wavenumber, 4, 8, 1, 2)

        # Create fitting parameters row
        self.layout.label_slitsize = QLabel("Fitting parameters: \t\t\t\t\t Slit size: ")
        self.layout.label_slitsize.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_slitsize = QLineEdit("0.2768186233821143")
        self.layout.line_edit_slitsize.textChanged.connect(lambda: self.text_changed("Slit size"))
        self.layout.line_edit_slitsize.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        self.layout.label_k0 = QLabel("k0 (if known): ")
        self.layout.label_k0.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_k0 = QLineEdit("-0.001120285207139915")
        self.layout.line_edit_k0.textChanged.connect(lambda: self.text_changed("k0"))
        self.layout.label_k1 = QLabel("k1 (if known): ")
        self.layout.label_k1.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.line_edit_k1 = QLineEdit("1.0001065125867772")
        self.layout.line_edit_k1.textChanged.connect(lambda: self.text_changed("k1"))
        self.layout.button_fit_parameters = QPushButton("Save fitting parameters")
        self.layout.button_fit_parameters.clicked.connect(self.button_pressed)
        self.layout.label_fit_parameters_saved = QLabel("No fitting parameters saved yet")
        self.layout.label_fit_parameters_saved.setAlignment(Qt.AlignCenter)

        self.layout.addWidget(self.layout.label_slitsize, 5, 0,1,2)
        self.layout.addWidget(self.layout.line_edit_slitsize, 5, 2,1,1)
        self.layout.addWidget(self.layout.label_k0, 5, 3,1,1)
        self.layout.addWidget(self.layout.line_edit_k0, 5, 4,1,1)
        self.layout.addWidget(self.layout.label_k1, 5, 5,1,1)
        self.layout.addWidget(self.layout.line_edit_k1, 5, 6, 1, 1)
        self.layout.addWidget(self.layout.button_fit_parameters, 5, 7, 1, 1)
        self.layout.addWidget(self.layout.label_fit_parameters_saved, 5, 8,1,2)


        # Create necessary stuff for the fitting Row
        self.layout.label_fitting_stuff = QLabel("Fitting extra's: ")
        self.layout.label_fitting_stuff.setStyleSheet("font: bold ; font-size:10pt")
        self.layout.check_box_background = QCheckBox("Is the first loaded file a background? ")
        self.layout.check_box_background.toggled.connect(self.check_box_background)
        self.layout.check_box_lnc_mct_detector = \
            QCheckBox("Is a MCT (lnc) detector used (apply heating effects in the fitting procedure)? ")
        self.layout.check_box_lnc_mct_detector.toggled.connect(self.check_box_lnc_mct_detector)

        self.layout.addWidget(self.layout.label_fitting_stuff, 6, 0, 1, 2)
        self.layout.addWidget(self.layout.check_box_background, 6, 2, 1, 4)
        self.layout.addWidget(self.layout.check_box_lnc_mct_detector, 6, 6, 1, 4)

        self.layout.layout_grid_title = QHBoxLayout()
        self.layout.layout_grid_title.grid_label = QLabel("Fittable Molecules")
        self.layout.layout_grid_title.grid_label.setAlignment(Qt.AlignCenter)
        self.layout.layout_grid_title.addWidget(
            self.layout.layout_grid_title.grid_label)
        self.layout.layout_grid_title.grid_label.setStyleSheet("font: bold ; font-size:12pt")
        self.layout.layout_grid_title.grid_label.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)

        self.layout.addLayout(self.layout.layout_grid_title,7,0,1,10)

        self.layout.grid_molecules = QGridLayout()
        self.dict_molecules={}
        for mol_i in range(len(self.molecule_list)):
            column = int(mol_i/self.rows)
            mol = self.molecule_list[mol_i]
            checkbox_temp = QCheckBox(mol)
            self.dict_molecules[mol] = checkbox_temp
            self.layout.grid_molecules.addWidget(self.dict_molecules[mol],mol_i%4,column)
            self.dict_molecules[mol].toggled.connect(self.check_box_molecules)
        self.layout.addLayout(self.layout.grid_molecules,8,0,2,8)


        self.layout.button_ftir = QPushButton("Show data in directory")
        self.layout.button_ftir.clicked.connect(self.button_pressed)
        self.signal_create_inner_tab.connect(self.parent.tab_ftir_fitting.inner_tab.create_inner_tabs)
        self.layout.addWidget(self.layout.button_ftir, 10, 0, 1, 4)

        self.layout.empty_space = QLabel("")
        self.layout.addWidget(self.layout.empty_space, 11, 0, 2, 10)

        self.layout.addWidget(self.parent.tab_ftir_fitting.inner_tab, 13, 0, 11, 10)

        self.layout.empty_space2 = QLabel("")
        self.layout.addWidget(self.layout.empty_space2, 24, 0, 1, 10)

        self.parent.tab_ftir_fitting.setLayout(self.layout)

    def text_changed(self,s):
        if s == "Directory":
            self.layout.label_directory2.setText("Directory has been chosen, not yet new folder created")
        elif s == "Temperature":
            self.layout.label_temperature2.setText("Temperature changed, not yet stored")
        elif s == "Pressure":
            self.layout.label_pressure2.setText("Pressure changed, not yet stored")
        elif s == "Pathlength":
            self.layout.label_pathlength2.setText("Path length changed, not yet stored")
        elif s== "Minimum wavenumber":
            self.layout.label_wavenumber.setText("Minimum wavenumber changed, no wavenumbers stored yet")
        elif s=="Maximum wavenumber":
            self.layout.label_wavenumber.setText("Maximum wavenumber changed, no wavenumbers stored yet")
        elif s == "Slit size":
            self.layout.label_fit_parameters_saved.setText("Slit size was changed")
        elif s == "k0":
            self.layout.label_fit_parameters_saved.setText("k0 was changed")
        elif s == "k1":
            self.layout.label_fit_parameters_saved.setText("k1 was changed")

    def button_pressed(self,s):
        if self.sender().text() == "Create directory":
            self.directory_data_InvenioR = self.layout.line_edit_directory.text().replace("\\", "\\\\") + "\\\\"
            self.directory_save_invenioR_processed = self.layout.line_edit_directory.text() + "\\Processed_Data"

            if self.directory_save_invenioR_processed == "":
                self.layout.label_directory2.setText("No directory given.")
            elif platform.system() == "Windows":
                if ("C:\\" not in self.directory_data_InvenioR and "P:\\" not in self.directory_data_InvenioR and
                      "D:\\" not in self.directory_data_InvenioR and "F:\\" not in self.directory_data_InvenioR and
                      "K:\\" not in self.directory_data_InvenioR and "L:\\" not in self.directory_data_InvenioR):
                    self.layout.label_directory2.setText("Incorrect directory given.")
                else:
                    if not os.path.exists(self.directory_save_invenioR_processed):
                        os.makedirs(self.directory_save_invenioR_processed)
                        self.layout.label_directory2.setText("Dirctory exists.")

                    if os.path.exists(self.directory_save_invenioR_processed):
                        self.layout.label_directory2.setText("Directory exists.")
                    self.new_directory_bool = True
            elif platform.system() == "Darwin":
                if not os.path.exists(self.directory_save_invenioR_processed):
                    os.makedirs(self.directory_save_invenioR_processed)
                    self.layout.label_directory2.setText("Dirctory exists.")

                if os.path.exists(self.directory_save_invenioR_processed):
                    self.layout.label_directory2.setText("Directory exists.")
                self.new_directory_bool = True
        if self.sender().text() == "Save temperature data":
            try:
                self.list_ftir_temperature = []
                s_list = self.layout.line_edit_temperature.text().split(sep=",")
                for s_i in range(len(s_list)):
                    self.list_ftir_temperature.append(float(s_list[s_i]))
                self.layout.label_temperature2.setText("Temperature data has been stored.")
            except:
                self.layout.label_temperature2.setText("Incorrect temperature data given.")
        if self.sender().text() == "Save pressure data":
            try:
                self.list_ftir_pressure = []
                s_list = self.layout.line_edit_pressure.text().split(sep=",")
                for s_i in range(len(s_list)):
                    self.list_ftir_pressure.append(float(s_list[s_i]))
                self.layout.label_pressure2.setText("Pressure data has been stored.")
            except:
                self.layout.label_pressure2.setText("Incorrect pressure data given.")
        if self.sender().text() == "Save path length":
            try:
                self.list_ftir_pathlength = []
                s_list = self.layout.line_edit_pathlength.text().split(sep=",")
                for s_i in range(len(s_list)):
                    self.list_ftir_pathlength.append(float(s_list[s_i]))
                self.layout.label_pathlength2.setText("Pathlength has been stored.")
            except:
                self.layout.label_pathlength2.setText("Invalid pathlength given.")
        if self.sender().text() == "Show data in directory" and self.new_directory_bool==True:
            self.signal_create_inner_tab.emit([self.parent])
            self.new_directory_bool=False
        if self.sender().text() == "Save wavenumbers":
            try:
                self.wavenumber_min= float(self.layout.line_edit_wavenumber_min.text())
                self.wavenumber_max= float(self.layout.line_edit_wavenumber_max.text())
                self.layout.label_wavenumber.setText("Both wavenumbers stored")
            except:
                self.layout.label_wavenumber.setText("Uncorrect minimum or maximum wavenumber given")
        if self.sender().text() == "Save fitting parameters":
            try:
                self.slit_size = float(self.layout.line_edit_slitsize.text())
                self.k0 = float(self.layout.line_edit_k0.text())
                self.k1 = float(self.layout.line_edit_k1.text())
                self.layout.label_fit_parameters_saved.setText('Fit parameters stored')
            except:
                self.layout.label_fit_parameters_saved.setText('Uncorrect fit parameters given')


    def send_button_signal(self,s):
        if s=="Show data in directory" and self.new_directory_bool==True:
            self.signal_create_inner_tab.emit([self.parent])
            self.new_directory_bool=False

    def check_box_molecules(self,s):
        if s:
            self.molecule_storage_for_fitting.append(self.sender().text())
        if not s:
            self.molecule_storage_for_fitting.remove(self.sender().text())

    def check_box_background(self,s):
        if s:
            self.background_bool = True
        if not s:
            self.background_bool = False

    def check_box_lnc_mct_detector(self,s):
        if s:
            self.lnc_mct_detector_bool = True
        if not s:
            self.lnc_mct_detector_bool = False


class create_inner_tab_ftir_fitting(QTabWidget):
    signal_plot_fitting_and_residual = Signal(list)
    signal_button_fit_single = Signal(list)
    signal_button_fit_single_redo = Signal(list)
    signal_button_fit_all = Signal(list)
    signal_button_fit_all_redo = Signal(list)
    signal_save_data_in_text_files = Signal(list)
    signal_save_data_as_png = Signal(list)
    signal_ftir_fitting_tr_vs_ab = Signal(list)

    def __init__(self):
        super().__init__()
        self.list_plots = []
    def create_inner_tabs(self, parent):
        self.parent = parent[0]
        self.signal_plot_fitting_and_residual.connect(self.parent.worker_plotting.fitting_and_residual_plot)
        self.directory = self.parent.tab_ftir_fitting.directory_data_InvenioR
        self.directory_for_saving = self.parent.tab_ftir_fitting.directory_save_invenioR_processed
        self.files_in_directory_temp = [f for f in os.listdir(self.directory) if op.isfile(op.join(self.directory, f))]
        self.files_in_directory = []
        for file in self.files_in_directory_temp:
            if file[-2:] == ".0":
                self.files_in_directory.append(file)
            elif file[-4:] == ".csv":
                self.files_in_directory.append(file)
            elif file[-4:] == ".dat":
                self.files_in_directory.append(file)
            elif file[-4:] == ".tsv":
                self.files_in_directory.append(file)
        self.array_data = GFL.read_opus_data_from_folder_into_array_for_gui(self.directory)
        self.list_tabs = []
        self.w_dict = {}
        self.t_dict = {}
        self.background = ""

        try:
            for f in range(self.files_in_directory_previous):
                self.removeTab(0)
            self.files_in_directory_previous = len(self.files_in_directory)

            for i in range(len(self.array_data)):
                if i == 0 and self.parent.tab_ftir_fitting.background_bool:
                    self.background = self.files_in_directory_temp[i]
                try:
                    data_wavelength = np.array(self.array_data[i].get_range("AB")[0:-1])
                    if np.round(np.mean(self.array_data[i]["AB"][0:-1])) == 0:
                        self.array_data[i]["AB"][0:-1] = GFL.ab_to_tr(self.array_data[i]["AB"][0:-1])

                    data_transmission = np.array(self.array_data[i]["AB"][0:-1])
                except:
                    data_wavelength = self.array_data[i][0][0:-1]
                    if np.round(np.mean(self.array_data[i][1][0:-1])) == 0:
                        self.array_data[i][1][0:-1] = GFL.ab_to_tr(self.array_data[i][1][0:-1])
                    data_transmission = self.array_data[i][1][0:-1]


                s_temp = rd.Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                     units={"transmittance": ""})

                self.w_dict[self.files_in_directory_temp[i]], self.t_dict[self.files_in_directory_temp[i]] = s_temp.get("transmittance")

                tab_temp = QWidget()
                tab_temp_layout = QVBoxLayout()
                plot_temp = self.parent.worker_plotting.create_empty_plot("ftir_fitting_" + str(self.files_in_directory[i]))
                self.parent.worker_plotting.fitting_and_residual_plot(
                    ["ftir_fitting_" + str(self.files_in_directory[i]), self.w_dict[self.files_in_directory_temp[i]],
                     self.t_dict[self.files_in_directory_temp[i]]])

                tab_temp_layout.addWidget(plot_temp)
                tab_temp.setLayout(tab_temp_layout)
                self.addTab(tab_temp, self.files_in_directory[i])
        except:
            for i in range(len(self.array_data)):
                if i == 0 and self.parent.tab_ftir_fitting.background_bool:
                    self.background = self.files_in_directory_temp[i]
                try:
                    data_wavelength = np.array(self.array_data[i].get_range("AB")[0:-1])
                    if np.round(np.mean(self.array_data[i]["AB"][0:-1])) == 0:
                        self.array_data[i]["AB"][0:-1] = GFL.ab_to_tr(self.array_data[i]["AB"][0:-1])

                    data_transmission = np.array(self.array_data[i]["AB"][0:-1])
                except:
                    data_wavelength = self.array_data[i][0][0:-1]
                    if np.round(np.mean(self.array_data[i][1][0:-1])) == 0:
                        self.array_data[i][1][0:-1] = GFL.ab_to_tr(self.array_data[i][1][0:-1])
                    data_transmission = self.array_data[i][1][0:-1]


                s_exp = rd.Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                 units={"transmittance": ""})
                self.w_dict[self.files_in_directory_temp[i]], self.t_dict[self.files_in_directory_temp[i]] = s_exp.get(
                    "transmittance")

                tab_temp = QWidget()
                tab_temp_layout = QVBoxLayout()
                plot_temp = self.parent.worker_plotting.create_empty_plot(
                    "ftir_fitting_" + str(self.files_in_directory[i]))
                self.parent.worker_plotting.fitting_and_residual_plot(
                    ["ftir_fitting_" + str(self.files_in_directory[i]), self.w_dict[self.files_in_directory_temp[i]],
                     self.t_dict[self.files_in_directory_temp[i]]])

                tab_temp_layout.addWidget(plot_temp)
                tab_temp.setLayout(tab_temp_layout)
                self.addTab(tab_temp, self.files_in_directory[i])

            self.parent.tab_ftir_fitting.layout.button_fit_single = QPushButton("Start single fit")
            self.parent.tab_ftir_fitting.layout.button_fit_single.clicked.connect(
                lambda: self.button_pressed("Start single fit"))
            self.signal_button_fit_single.connect(self.parent.worker_ftir_fitting.ftir_fit_one)
            self.parent.tab_ftir_fitting.layout.addWidget(self.parent.tab_ftir_fitting.layout.button_fit_single,
                                                          10, 4, 1, 3)

            self.parent.tab_ftir_fitting.layout.button_fit_single_redo = QPushButton("Redo single fit")
            self.parent.tab_ftir_fitting.layout.button_fit_single_redo.clicked.connect(
                lambda: self.button_pressed("Redo single fit"))
            self.signal_button_fit_single_redo.connect(self.parent.worker_ftir_fitting.ftir_fit_one)
            self.parent.tab_ftir_fitting.layout.addWidget(self.parent.tab_ftir_fitting.layout.button_fit_single_redo,
                                                          10, 7, 1, 3)

            self.parent.worker_ftir_fitting.signal_fitting_plot.connect(
                self.parent.worker_plotting.update_fitting_and_residual_plot)

            self.parent.tab_ftir_fitting.layout.label_info_fit = QLabel("No fitting done yet.")
            self.parent.tab_ftir_fitting.layout.label_info_fit.setAlignment(Qt.AlignCenter)
            self.parent.tab_ftir_fitting.layout.addWidget(self.parent.tab_ftir_fitting.layout.label_info_fit,
                                                          11, 0, 1, 4)

            self.parent.tab_ftir_fitting.layout.button_fit_all = QPushButton("Start fitting all")
            self.parent.tab_ftir_fitting.layout.button_fit_all.clicked.connect(
                lambda: self.button_pressed("Start fitting all"))
            self.signal_button_fit_all.connect(self.parent.worker_ftir_fitting.ftir_fit_all)
            self.parent.tab_ftir_fitting.layout.addWidget(self.parent.tab_ftir_fitting.layout.button_fit_all,
                                                          11, 4, 1, 3)

            self.parent.tab_ftir_fitting.layout.button_fit_all_redo = QPushButton("Redo fitting all")
            self.parent.tab_ftir_fitting.layout.button_fit_all_redo.clicked.connect(lambda: self.button_pressed(
                "Redo fitting all"))
            self.signal_button_fit_all_redo.connect(self.parent.worker_ftir_fitting.ftir_fit_all)
            self.parent.tab_ftir_fitting.layout.addWidget(self.parent.tab_ftir_fitting.layout.button_fit_all_redo,
                                                          11, 7, 1, 3)

            self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance = QPushButton("Showing transmission")
            self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.clicked.connect(
                lambda: self.button_pressed("Change transmission/absorption"))
            self.signal_ftir_fitting_tr_vs_ab.connect(self.parent.worker_plotting.update_fitting_and_residual_plot)
            self.parent.tab_ftir_fitting.layout.addWidget(
                self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance,
                12, 7, 1, 3)

            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_in_txt = \
                QPushButton("Save fitted data in text file")
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_in_txt.clicked.connect(
                lambda: self.button_pressed("Save fitted data in text file"))
            self.signal_save_data_in_text_files.connect(self.parent.worker_saving.save_data_in_txt)
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_as_png = QPushButton("Save fitted data as png")
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_as_png.clicked.connect(
                lambda: self.button_pressed("Save fitted data as png"))
            self.signal_save_data_as_png.connect(self.parent.worker_saving.save_data_as_png)
            self.parent.tab_ftir_fitting.layout.label_save = QLabel("No data saved yet")
            self.parent.tab_ftir_fitting.layout.label_save.setAlignment(Qt.AlignCenter)
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_in_txt_all = \
                QPushButton("Save all fitted data in text files")
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_in_txt_all.clicked.connect(
                lambda: self.button_pressed("Save all fitted data in text files"))
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_as_png_all = \
                QPushButton("Save all fitted data as png's")
            self.parent.tab_ftir_fitting.layout.button_save_simulated_data_as_png_all.clicked.connect(
                lambda: self.button_pressed("Save all fitted data as png's"))

            self.parent.tab_ftir_fitting.layout.addWidget(
                self.parent.tab_ftir_fitting.layout.button_save_simulated_data_in_txt, 24, 0, 1, 2)
            self.parent.tab_ftir_fitting.layout.addWidget(
                self.parent.tab_ftir_fitting.layout.button_save_simulated_data_as_png, 24, 2, 1, 2)
            self.parent.tab_ftir_fitting.layout.addWidget(
                self.parent.tab_ftir_fitting.layout.button_save_simulated_data_in_txt_all, 24, 4, 1, 2)
            self.parent.tab_ftir_fitting.layout.addWidget(
                self.parent.tab_ftir_fitting.layout.button_save_simulated_data_as_png_all, 24, 6, 1, 2)
            self.parent.tab_ftir_fitting.layout.addWidget(self.parent.tab_ftir_fitting.layout.label_save, 24, 8, 1, 2)
            self.files_in_directory_previous = len(self.files_in_directory)

    def button_pressed(self, s):
        if s == "Start single fit":
            self.signal_button_fit_single.emit([False])
        elif s == "Redo single fit":
            self.signal_button_fit_single_redo.emit([True])
        elif s == "Start fitting all":
            self.signal_button_fit_all.emit([False])
        elif s == "Redo fitting all":
            self.signal_button_fit_all_redo.emit([True])
        elif s == "Save fitted data in text file":
            files_in_directory = self.parent.tab_ftir_fitting.inner_tab.files_in_directory
            index = self.parent.tab_ftir_fitting.inner_tab.currentIndex()
            current_file = files_in_directory[index]
            bool_fit = all(i == 0 for i in self.parent.worker_plotting.dict_plots[
                "ftir_fitting_" + current_file].p2.y_res)
            if not bool_fit:
                self.signal_save_data_in_text_files.emit([self.parent.tab_ftir_fitting.directory_save_invenioR_processed,
                                                          current_file,
                                                          np.transpose([self.parent.worker_plotting.dict_plots[
                                                                            "ftir_fitting_" + current_file].p1.x_fit,
                                                                        self.parent.worker_plotting.dict_plots[
                                                                            "ftir_fitting_" + current_file].p1.y_exp,
                                                                        self.parent.worker_plotting.dict_plots[
                                                                            "ftir_fitting_" + current_file].p1.y_fit,
                                                                        self.parent.worker_plotting.dict_plots[
                                                                            "ftir_fitting_" + current_file].p2.y_res]),
                                                          "Fitted FTIR Data"])
                self.parent.tab_ftir_fitting.layout.label_save.setText("Saved fitted data for '" + current_file +
                                                                       "' in text file")
            else:
                self.parent.tab_ftir_fitting.layout.label_save.setText(
                    "Nothing is saved for '" + current_file + "', no fitted data yet")
        elif s == "Change transmission/absorption":
            if self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.text() == "Showing transmission":
                files_in_directory = self.parent.tab_ftir_fitting.inner_tab.files_in_directory
                for index in range(len(files_in_directory)):
                    current_file = files_in_directory[index]

                    try:
                        self.parent.tab_ftir_fitting.absorbance_bool = True
                        self.signal_ftir_fitting_tr_vs_ab.emit(["ftir_fitting_" + current_file,
                                                                self.parent.tab_ftir_fitting.absorbance_bool])
                        self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.setText("Showing absorbance")
                    except:
                        self.parent.tab_ftir_fitting.absorbance_bool = True
                        self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.setText("Showing absorbance")
            elif self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.text() == "Showing absorbance":
                files_in_directory = self.parent.tab_ftir_fitting.inner_tab.files_in_directory
                for index in range(len(files_in_directory)):
                    current_file = files_in_directory[index]
                    try:
                        self.parent.tab_ftir_fitting.absorbance_bool = False
                        self.signal_ftir_fitting_tr_vs_ab.emit(["ftir_fitting_" + current_file,
                                                                self.parent.tab_ftir_fitting.absorbance_bool])
                        self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.setText("Showing transmission")
                    except:
                        self.parent.tab_ftir_fitting.absorbance_bool = False
                        self.parent.tab_ftir_fitting.layout.button_transmission_vs_absorbance.setText("Showing transmission")
        elif s == "Save all fitted data in text files":
            files_in_directory = self.parent.tab_ftir_fitting.inner_tab.files_in_directory
            failed_fit_list = []
            for f in range(len(files_in_directory)):
                current_file = files_in_directory[f]
                bool_fit = all(i == 0 for i in self.parent.worker_plotting.dict_plots[
                    "ftir_fitting_" + current_file].p2.y_res)
                if not bool_fit:
                    self.signal_save_data_in_text_files.emit(
                        [self.parent.tab_ftir_fitting.directory_save_invenioR_processed,
                         current_file,
                         np.transpose([self.parent.worker_plotting.dict_plots[
                                           "ftir_fitting_" + current_file].p1.x_fit,
                                       self.parent.worker_plotting.dict_plots[
                                           "ftir_fitting_" + current_file].p1.y_exp,
                                       self.parent.worker_plotting.dict_plots[
                                           "ftir_fitting_" + current_file].p1.y_fit,
                                       self.parent.worker_plotting.dict_plots[
                                           "ftir_fitting_" + current_file].p2.y_res]),
                         "Fitted FTIR Data"])
                    self.parent.tab_ftir_fitting.layout.label_save.setText("Saved all fitted data for in text files")
                else:
                    failed_fit_list.append(current_file)
                    self.parent.tab_ftir_fitting.layout.label_save.setText(
                        "Nothing is saved for '" + str(failed_fit_list) + "', no fitted data yet")
        elif s == "Save fitted data as png":
            files_in_directory = self.parent.tab_ftir_fitting.inner_tab.files_in_directory
            index = self.parent.tab_ftir_fitting.inner_tab.currentIndex()
            current_file = files_in_directory[index]
            bool_fit = all(i == 0 for i in self.parent.worker_plotting.dict_plots[
                "ftir_fitting_" + current_file].p2.y_res)

            if not bool_fit:
                self.signal_save_data_as_png.emit([
                    self.parent.tab_ftir_fitting.directory_save_invenioR_processed, current_file,
                    [self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p1.x_fit,
                     self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p1.y_exp,
                     self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p1.y_fit,
                     self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p2.y_res],
                    "Fitted FTIR Data"])
                self.parent.tab_ftir_fitting.layout.label_save.setText("Saved fitted data for '" + current_file +
                                                                       "' as png")
            else:
                self.parent.tab_ftir_fitting.layout.label_save.setText(
                    "Nothing is saved for '" + current_file + "', no fitted data yet")
        elif s == "Save all fitted data as png's":
            failed_fit_list = []
            files_in_directory = self.parent.tab_ftir_fitting.inner_tab.files_in_directory
            for f in range(len(files_in_directory)):
                current_file = files_in_directory[f]
                bool_fit = all(i == 0 for i in self.parent.worker_plotting.dict_plots[
                    "ftir_fitting_" + current_file].p2.y_res)
                if not bool_fit:
                    self.signal_save_data_as_png.emit([
                        self.parent.tab_ftir_fitting.directory_save_invenioR_processed,current_file,
                        [self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p1.x_fit,
                         self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p1.y_exp,
                         self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p1.y_fit,
                         self.parent.worker_plotting.dict_plots["ftir_fitting_" + current_file].p2.y_res],
                        "Fitted FTIR Data"])
                    self.parent.tab_ftir_fitting.layout.label_save.setText("Saved all fitted data as png's")
                else:
                    failed_fit_list.append(current_file)
                    self.parent.tab_ftir_fitting.layout.label_save.setText(
                        "Nothing is saved for '" + str(failed_fit_list) + "', no fitted data yet")


class create_tab_ftir_fit_analysis(QObject):
    def __init__(self):
        super().__init__()

    def create(self):
        self.tab_ftir_fit_analysis = QWidget()
        self.tab_ftir_fit_analysis.layout_full = QGridLayout()


