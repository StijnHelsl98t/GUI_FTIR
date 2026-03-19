"""
Main file. In this file, a single class is used to create and start the needed QThreads (to make the GUI faster and more
responsive), after which the GUI is created (see GUI_Tab_Library.py to see what Tabs are created within the GUI).

Author: Stijn Helsloot (stijn.helsloot@maastrichtuniversity.nl)
"""

from PySide6.QtWidgets import QWidget, QApplication, QTabWidget, QVBoxLayout
from PySide6.QtGui import QIcon
from PySide6.QtCore import QThread
import sys
sys.path.append("..")
import Code_Gui.Gui_Library.Gui_Workers as GW
import Code_Gui.Gui_Library.Gui_Tab_Library as GTL
import pyqtgraph as pg

import matplotlib
matplotlib.use("Agg")

class Window(QWidget):
    def __init__(self, parent=None):
        super().__init__()
        self.setWindowTitle("FTIR GUI")
        self.UI()
        self.showMaximized()

    def UI(self):
        # Create worker for simulation
        self.worker_ftir_simulation = GW.WorkerFTIRSimulation(self)
        self.thread_ftir_simulation = QThread()
        self.thread_ftir_simulation.start()
        self.worker_ftir_simulation.moveToThread(self.thread_ftir_simulation)

        # Create worker to plot data
        self.worker_plotting = GW.WorkerPlotter(self)
        
        # Create worker for hitran
        self.worker_hitran = GW.WorkerHitran(self)
        self.thread_hitran = QThread()
        self.thread_hitran.start()
        self.worker_hitran.moveToThread(self.thread_hitran)

        # Create worker to fit simulated data to experimental data
        self.worker_ftir_fitting = GW.WorkerFTIRFitter(self)
        self.thread_ftir_fitting = QThread()
        self.thread_ftir_fitting.start()
        self.worker_ftir_fitting.moveToThread(self.thread_ftir_fitting)
        

        # Create worker to save data
        self.worker_saving = GW.WorkerSaver(self)
        self.thread_saver = QThread()
        self.thread_saver.start()
        self.worker_saving.moveToThread(self.thread_saver)

        # Create main GUI
        mainLayout = QVBoxLayout()
        self.tabs = QTabWidget()

        
        # Create tab for getting / updating database
        self.tab_get_database = GTL.create_tab_database_updater()
        self.tab_get_database.create_tab(self)

        # Create tab for simulations
        self.tab_ftir_simulator  = GTL.create_tab_ftir_simulator()
        self.tab_ftir_simulator.create_tab(self)

        # Create tab for fitting procedure
        self.tab_ftir_fitting = GTL.create_tab_ftir_fitting()
        self.tab_ftir_fitting.inner_tab = GTL.create_inner_tab_ftir_fitting()
        self.tab_ftir_fitting.create_tab(self)

        self.tabs.addTab(self.tab_get_database, "Get/Update Hitran Database")
        self.tabs.addTab(self.tab_ftir_simulator, "Simulate FTIR Data")
        self.tabs.addTab(self.tab_ftir_fitting, "Fit FTIR Data")

        mainLayout.addWidget(self.tabs)
        self.setLayout(mainLayout)

        self.setWindowIcon(QIcon('Spherical_Cow_NoBackground.ico'))



if __name__== "__main__":
    App = QApplication(sys.argv)
    window = Window()
    App.exec()