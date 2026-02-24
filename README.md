# FTIR GUI
Fully open source Gui made for FTIR simulation/data analysis. Both are done via using the open-source python module RADIS.

# Installation
## Windows
- Install the requirements: `pip install -r requirements.txt`
- Install the package: `pip install -e .`
- Run the GUI: `python -m Code_Gui.FTIR_GUI`

## MacOS
- Install Homebrew if you don't have it already.
- Install c-blosc and HDF5: `brew install c-blosc hdf5`
- Install a virtual environment: `python3 -m venv .venv` (Python 3.14.2 tested on MacOS Sequoia) 
- Activate the virtual environment: `source .venv/bin/activate`
- Install the requirements: `pip install -r requirements.txt`
- Install the package: `pip install -e .`
- Run the GUI: `python -m Code_Gui.FTIR_GUI`

# Purpose
This GUI can be used for both for the simulating of theoretical FTIR spectra or for the fitting to experimental FTIR spectra. For both, a tab has been created within the GUI. 

## Simulation Tab
The simulation tab can be seen in the figure below.

![Tab1](https://github.com/StijnHelsl98t/GUI_FTIR/assets/133780753/ae29160f-8db5-4146-b377-2d4f056d9a5e)

As one can see, a quick simulation has been done here for a FTIR spectra with water, carbon-dioxide and carbon-monoxide. The parameters that one can change (as of now) are:
- **Temperature:** The temperature of the gas decides the final spectra. This effect can be simulated using the temperature-parameter.
- **Pressure:** The pressure of the gas decides the final spectra, as a lower pressure will mean that there are less particles that can absorb any of the light. This effect can be simulated using the pressure-parameter.
- **Pathlength:** The distance that light has to travel through gas matters for how much light is absorbed, deciding the final spectra. This effect can be simulated using the pathlength-parameter. 
- **Slitsize:** Within actual experiments, the shape of your spectra are decided by your chosen aperature. This effect can be simulated using the slitsize-paramter. If left empty, no slitsize is used, giving the perfect theoretical spectrum.
- **Concentration:** The actual gas composition of course decides the final spectra. In total, 24 molecules have been added to this GUI. As we are using Hitran, we are capable of adding a total of 55 molecules. However, as most of these are not useful for us, we decided to leave them out. Those that want to use this GUI but need some of these other molecules, are completely allowed to change my python code and add/remove the molecules they want to add/remove.

One can also save the created spectra, both in txt and picture format. For this, one needs to first add a directory where one wants to save this files.

## Fitting Tab
The fitting tab can be seen in the figure below.

![Tab2](https://github.com/StijnHelsl98t/GUI_FTIR/assets/133780753/2adb9b76-8dec-47ff-af3f-f9f60636fd79)

This tab is used, when one wants to know the concentration of the species within their gas. In order to fit, the following information/parameters are required:
- **Directory:** One needs to add the directory where your data has been saved. As of this moment, only .0- (bruker) and .csv-files are seen as a possible data-file. Within later additions of this GUI, we can add more option (when required).
- **Temperature:** One needs to add the temperature at which the spectrum has been taken.
- **Pressure:** One needs to add the pressure at which the spectrum has been taken
- **Pathlength:** One needs to add the pathlength of the gas-cell used within their experiments. The standard value is 20.062, which is the value of our own gas-cell. If one choses to re-create this GUI, they can of course change this as well.
- **Slitsize:** One needs to add the slitsize as well. This is a way of simulating the apperature used within our experiments. The standard value of this slitsize is 0.2768..., which is the slitsize for our own FTIR spectrometer (when we use an apperature of 1.5 mm). Important to note: this value needs to be gained via calibration of your FTIR, so the properness at which this is done, will decide all further data.
- **Offset:** No machine (FTIR spectrometer) is perfect, so sometimes they have a slight offset compared to the perfect theoretical spectra. For this reason, one can add this offset as well (can be either a constant value (choose offset left and right as the same value) or a linear function (draws a straight line between offset left and right)). Again, this offset has to be calibrated for your machine, so this needs to be done with certain care.
- **LNC MCT detector:** It's possible to use an liquid nitrogen cooled (LNC) MCT detector. Using one, one will note that a certain peak seems to pop up around 3135 cm$^{-1}$. This peak exists due to ice forming on the detector due to the rapid cooling, which will slowly melt away during your experiments. This ice-peak can be very annoying if one has to fit molecules which lie inside of this range (for example CH$_{4}$). For this, we created the "ice-peak fitting"-function, which can be used to fit the ice-peak and remove this from the spectrum during the fitting of the species within our gas.
- **Background:** If one notes their background has a lot of carbon-dioxide or water, one can use the "remove-background"-function to remove this from the fitting procedure.

If one wants to fit all their spectra in one sitting, this is also possible. One only needs to add all their temperature/pressure data at once (by the following method 1,2,3,4,5,....). Make sure the created pressure/temperature list is of equal length as the amount of data-files one has.

After fitting, the graph shown on the bottom is created, which includes the original spectrum, the created fit and the residual.

Once done, one can choose to save this graph in either txt or png format. If one chooses to save this spectrum as a txt-file, the next time one presses the "start single fit"-button, it will not fit but instead again show the already created fit (even after closing and re-opening). If one is not satisfied with this fit, click the "Redo single fit"-button. This will ignore the previously created fit and redo it. The fitted concentrations are saved inside of an HDF5 file. If one wants to open this file, one needs to install a HDF5 viewer. Please see (download HDF)[https://www.hdfgroup.org/downloads/].

# Usage
There are three possible options for someone to use this GUI:
- The first option includes e-mailing the creator and asking if you can get an installer and use this to install the gui, to which he will happily comply. This email is the following: stijn.helsloot@maastrichtuniversity.nl.
- The second option would be to actually create this installer/gui.exe themselves, using the "setup,py"-file and python. To do this, one would need to install all the files within the folder and copy-paste the files inside into a single python environment. Once this is done, one can run the "setup,py"-file to get the wanted installer. The creation of this installer might take a couple of hours, so please keep this in mind. Also important to keep in mind, this will need one needs to download all the necessary modules.
-  The last option would be to us python to directly create the GUI. For this, one needs to install al the files within the folder and copy-paste the files inside into a single python environment. Once this is done, once can run the "FTIR-GUI.py"-file. After a few second, this will open the GUI. Important to note: everytime one wants to open this GUI now, one will need to do this via python. No stand-alone GUI is now created. Also, again, one needs to make sure the necessary modules are downloaded.
