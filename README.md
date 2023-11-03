# GUI_FTIR
Fully open source Gui made for FTIR simulation/data analysis (using python RADIS)

# Purpose
- Within my own project, we are focussing on the injection of trace methane behind a microwave plasma. The source gas for this microwave plasma will be either hydrogen or oxygen. The objective is to destroy as much of the methane as possible, while using as little energy as necesarry. 
- A colleague of mine is interested in using a mix of nitrogen and oxygen as an source gas for the microwave plasma.

Because of our focus, a GUI has been created that can both simulate FTIR spectra as well as fit calculated spectra to the gained experimental spectra. Within this GUI, in total 24 molecules have been chosen to be listed, which include all of the Hitran molecule-options which contain a combination of H, N, O and C. 

If one wants to run this GUI, there are multiple options to do this:
- The first option is the simplest, if one doesn't want to use python at all. Simply download the .msi file (installer) and use this to install the gui.
- The second option does use python. In this case, one should install the "Code Gui" folder and throw this into it's own environment. After this, install the necessary modules and run "FTIR-GUI.py". This will open the GUI. Important to note: everytime one wants to open this GUI now, one will need to do this via python. No stand-alone GUI is now created.
- The third possibility is to create this stand-alone GUI themselves via the "setup,py" file (again, python is needed).
