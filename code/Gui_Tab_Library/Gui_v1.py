"""
Main file. In this file, a single class is used to create and start the needed QThreads (to make the GUI faster), after
which the GUI is created (see GUI_Tab_Library.py to see what Tabs are created within the GUI).

Author: Stijn Helsloot (stijn.helsloot@maastrichtuniversity.nl)
"""

import sys
from PySide6.QtWidgets import QWidget, QApplication, QTabWidget, QVBoxLayout
from PySide6.QtCore import QThread
import Gui_Workers as GW
import Gui_Tab_Library as GTB


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
        self.worker_ftir_simulation.moveToThread(self.thread_ftir_simulation)
        self.thread_ftir_simulation.start()

        # Create worker to plot data
        self.worker_plotting = GW.WorkerPlotter(self)
        self.thread_plotting = QThread()
        self.worker_plotting.moveToThread(self.thread_plotting)
        self.thread_plotting.start()

        # Create worker to fit simulated data to experimental data
        self.worker_ftir_fitting = GW.WorkerFTIRFitter(self)
        self.thread_ftir_fitting = QThread()
        self.worker_ftir_fitting.moveToThread(self.thread_ftir_fitting)
        self.thread_ftir_fitting.start()

        # Create worker to save data
        self.worker_saving = GW.WorkerSaver(self)
        self.thread_saver = QThread()
        self.worker_saving.moveToThread(self.thread_saver)
        self.thread_saver.start()

        # Create main GUI
        mainLayout = QVBoxLayout()
        self.tabs = QTabWidget()

        # Create tab for simulations
        self.tab_ftir_simulator  = GTB.create_tab_ftir_simulator()
        self.tab_ftir_simulator.create_tab(self)

        # Create tab for fitting procedure
        self.tab_ftir_fitting = GTB.create_tab_ftir_fitting()
        self.tab_ftir_fitting.inner_tab = GTB.create_inner_tab_ftir_fitting()
        self.tab_ftir_fitting.create_tab(self)

        self.tabs.addTab(self.tab_ftir_simulator, "Simulate FTIR Data")
        self.tabs.addTab(self.tab_ftir_fitting, "Fit FTIR Data")

        mainLayout.addWidget(self.tabs)
        self.setLayout(mainLayout)



if __name__== "__main__":
    App = QApplication(sys.argv)
    window = Window()
    App.exec()