"""
To use this file and create an Gui, run the following line in your python environment:

- python Setup.py build (to directly create the gui)
- python Setup.py bdist_msi (to create an installer for Windows)
- python Setup.py bdist_mac (to create a .app bundle for macOS)
- python Setup.py bdist_dmg (to create a disk image for macOS)

"""

import sys
import os
import sysconfig
from cx_Freeze import setup, Executable
from setuptools import find_packages


# Dependencies are automatically detected, but it might need fine tuning.
# We use relative paths for local code and dynamic site-packages for libraries.
include_files = [
    ("Code_Gui/Gui_Calibration_Code/", "lib/Gui_Calibration_Code/"),
    ("Code_Gui/Gui_General_Code/", "lib/Gui_General_Code/"),
    ("Code_Gui/Gui_Library/", "lib/Gui_Library/"),
]

# List of libraries that were manually included in the Windows version
libs_to_include = [
    "tables.libs",
    "astropy",
    "pyqtgraph",
    "radis",
    "lmfit",
    "brukeropusreader",
    "hjson",
]

for lib in libs_to_include:
    # Check both purelib and platlib for the library in the current environment
    for path_type in ["purelib", "platlib"]:
        lib_path = os.path.join(sysconfig.get_paths()[path_type], lib)
        if os.path.exists(lib_path):
            include_files.append((lib_path, f"lib/{lib}"))
            break

build_exe_options = {
    "excludes": ["tkinter"],
    "zip_include_packages": ["encodings", "PySide6"],
    "include_files": include_files,
}

# base="Win32GUI" should be used only for Windows GUI app
base = "Win32GUI" if sys.platform == "win32" else None

# macOS specific options for creating a .app bundle
bdist_mac_options = {
    "bundle_name": "FTIR_GUI",
    "iconfile": "Spherical_Cow_NoBackground.ico",  # Note: .icns is preferred for Mac
}

setup(
    name="FTIR_GUI",
    version="1",
    optimize=2,
    description="A GUI that can be used to both simulate FTIR spectra or fit concentrations to experimental FTIR data.",
    packages=find_packages(),
    options={
        "build_exe": build_exe_options,
        "bdist_mac": bdist_mac_options,
    },
    executables=[
        Executable(
            "Code_Gui/FTIR_GUI.py", base=base, icon="Spherical_Cow_NoBackground.ico"
        )
    ],
)
