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
         "lib\\astropy\\") #Adding astropy, since otherwise not imported correctly
    ]
}

# base="Win32GUI" should be used only for Windows GUI app
base = "Win32GUI" if sys.platform == "win32" else None

setup(
    name="Gui_v1",
    version="0.1",
    optimize=2,
    description="My GUI application!",
    options={"build_exe":build_exe_options},
    executables=[Executable("Gui_v1.py", base=base)],)