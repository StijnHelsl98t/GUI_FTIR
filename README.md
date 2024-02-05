# FTIR GUI
Fully open source Gui made for FTIR simulation/data analysis. Both are done via using the open-source python module RADIS.

# Purpose
This GUI can be used for both for the simulating of theoretical FTIR spectra or for the fitting to experimental FTIR spectra. For both, a tab has been created within the GUI. 

![Tab1](https://github.com/StijnHelsl98t/GUI_FTIR/assets/133780753/ae29160f-8db5-4146-b377-2d4f056d9a5e)


Within my own project, we are focussing on the injection of trace methane behind a microwave plasma. The source gas for this microwave plasma will be either hydrogen or oxygen. The objective is to destroy as much of the methane as possible, while using as little energy as necesarry. 
- A colleague of mine is interested in using a mix of nitrogen and oxygen as an source gas for the microwave plasma.

Because of our focus, a GUI has been created that can both simulate FTIR spectra as well as fit calculated spectra to the gained experimental spectra. Within this GUI, in total 24 molecules have been chosen to be listed, which include all of the Hitran molecule-options which contain a combination of H, N, O and C. 

# Usage
There are three possible options for someone to use this GUI:
- The first option includes e-mailing the creator and asking if you can get the .msi file (installer) and use this to install the gui, to which he will happily comply. This email is the following: stijn.helsloot@maastrichtuniversity.nl
- The second option would be to actually create this installer themselves, using the "setup,py"-file and python. To do this, one would need to install all the files within the folder and copy-paste the files inside into a single python environment. Once this is done, one can run the "setup,py"-file to get the wanted installer. The creation of this installer might take a couple of hours, so please keep this in mind.
-  The last option would be to us python to directly create the GUI. For this, one needs to install al the files within the folder and copy-paste the files inside into a single python environment. Once this is done, once can run the "FTIR-GUI.py"-file. After a few second, this will open GUI. Important to note: everytime one wants to open this GUI now, one will need to do this via python. No stand-alone GUI is now created.
