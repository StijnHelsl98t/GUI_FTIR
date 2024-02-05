"""
To use this file and create an Gui, run the following line in your python environment:

- python Setup.py build (to directly create the gui)
- python Setup.py bdist_msi (to create an installer)

"""

import sys
from cx_Freeze import setup, Executable


# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {
    "excludes": ["tkinter"],
    "zip_include_packages": ["encodings", "PySide6"],
    "include_files":[
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\tables.libs",
         "lib\\tables.libs"), #Adding tables.libs, cause otherwise not imported
       ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\astropy\\",
         "lib\\astropy\\"), #Adding astropy, since otherwise not imported correctly
       ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\Code_Gui\\Gui_Calibration_Code\\",
         "lib\\Gui_Calibration_Code\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\Code_Gui\\Gui_General_Code\\",
         "lib\\Gui_General_Code\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\Code_Gui\\Gui_Library\\",
         "lib\\Gui_Library\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\pyqtgraph\\",
         "lib\\pyqtgraph\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\radis\\",
         "lib\\radis\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\lmfit\\",
         "lib\\lmfit\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\brukeropusreader\\",
         "lib\\brukeropusreader\\"),
        ("C:\\Users\\P70085588\\Documents\\General_CCE\\PycharmProjects\\Project_Gui\\venv\\Lib\\site-packages\\hjson\\",
         "lib\\hjson\\")
    ]
}

# base="Win32GUI" should be used only for Windows GUI app
base = "Win32GUI" if sys.platform == "win32" else None

setup(
    name="FTIR_GUI",
    version="3.0",
    optimize=2,
    description="A GUI that can be used to both simulate FTIR spectra or fit concentrations to experimental FTIR data.",
    options={"build_exe":build_exe_options},
    executables=[Executable("C:/Users/P70085588/Documents/General_CCE/PycharmProjects/Project_Gui/Code_Gui/FTIR_GUI.py", base=base, icon="Spherical_Cow_NoBackground.ico")],)