"""
To use this file and create an Gui, run the following line in your python environment:

- python Setup.py build (to directly create the gui)
- python Setup.py bdist_msi (to create an installer)

"""

import sys
from cx_Freeze import setup, Executable

import pathlib

directory_current = str(pathlib.Path().resolve())
directory_packages = directory_current + "\\venv\\Lib\\site-packages"
directory_GUI = directory_current + "\\Code_Gui"

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {
    "excludes": ["tkinter"],
    "zip_include_packages": ["encodings", "PySide6"],
    "include_files":[
        ( directory_packages + "\\tables.libs",
         "lib\\tables.libs"), #Adding tables.libs, cause otherwise not imported
       ( directory_packages + "\\astropy\\",
         "lib\\astropy\\"), #Adding astropy, since otherwise not imported correctly
       ( directory_GUI + "\\Gui_Calibration_Code\\",
         "lib\\Gui_Calibration_Code\\"),
        ( directory_GUI + "\\Gui_General_Code\\",
         "lib\\Gui_General_Code\\"),
        ( directory_GUI + "\\Gui_Library\\",
         "lib\\Gui_Library\\"),
        (directory_packages + "\\pyqtgraph\\",
         "lib\\pyqtgraph\\"),
        (directory_packages+ "\\radis\\",
         "lib\\radis\\"),
        (directory_packages + "\\lmfit\\",
         "lib\\lmfit\\"),
        (directory_packages + "\\brukeropusreader\\",
         "lib\\brukeropusreader\\"),
        (directory_packages + "\\hjson\\",
         "lib\\hjson\\")
    ]
}

# base="Win32GUI" should be used only for Windows GUI app
base = "Win32GUI" if sys.platform == "win32" else None

setup(
    name="FTIR_GUI",
    version="6.0",
    optimize=2,
    description="A GUI that can be used to both simulate FTIR spectra or fit concentrations to experimental FTIR data.",
    options={"build_exe":build_exe_options},
    executables=[Executable(directory_GUI + "/FTIR_GUI.py", base=base, icon="Spherical_Cow_NoBackground.ico")],)