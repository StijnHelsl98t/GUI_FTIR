import numpy as np
import brukeropusreader as bor
import radis as rd
import os
import os.path as op
import pandas as pd
from brukeropusreader import read_file
from bisect import bisect_left
import csv

def read_in_opus_file_with_keys(dir, file_name, read_keys=False):
    """
    Function created to read in opus-files. If the correct boolean is chosen, all possible keys within the opus-file
    will also be printed.

    :param dir: Directory to folder from which one wants to read in their opus-file.
    :param file_name: Name of their opus-file.
    :param read_keys: Boolean which decides if a table with all possible keys within the opus-file is printed or not.
    Default is false.
    :return: The data read in from the given opus-file.
    """
    opus_data = bor.read_file(dir + file_name)

    if read_keys:
        opus_data_keys_list = list(opus_data.keys())
        opus_data_keys_table = [[]]
        for key_i in range(len(opus_data_keys_list)):
            opus_data_keys_table[-1].append(opus_data_keys_list[key_i])
            if len(opus_data_keys_table[-1]) == 3:
                opus_data_keys_table.append([])
        opus_data_keys_table_printable = pd.DataFrame(opus_data_keys_table)
        print(f"The chosen dataset contains the following usable keys: ")
        print(opus_data_keys_table_printable)
    return opus_data

def read_opus_data_from_folder_into_array_for_gui(dir_inv):
    """
    The function gets the usable data-files from the given directory. The usable files are either OPUS-files, csv-files
    or dat-files (not yet implemented as not yet used).

    :param dir_inv: given directory where to find the data-files
    :return: array filled with the gained data-files from the given directory
    """
    directory_data = dir_inv
    files_in_directory_temp = [f for f in os.listdir(directory_data) if op.isfile(op.join(directory_data,f))]
    files_in_directory = []
    for file in files_in_directory_temp:
        if file[-2:] == ".0":
            files_in_directory.append(file)
        elif file[-4:] == ".csv":
            files_in_directory.append(file)
        elif file[-4:] == ".dat":
            files_in_directory.append(file)
    dat_array = []
    for file in files_in_directory:
        if file[-2:] == ".0":
            dat_array.append(read_file(directory_data + file))
        elif file[-4:] == ".csv":
            w_data = []
            t_data = []
            with open(directory_data + "\\" + file) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=",")
                line_count = 0
                for row in csv_reader:
                    w_data.append(float(row[0]))
                    t_data.append(float(row[1]))
            full_data = [np.array(w_data), np.array(t_data)]
            dat_array.append(full_data)
        #elif file[-4:] == ".dat":
        #    print("hello")
    return dat_array

def change_wavenumber_range(w_data, new_wmin = -1, new_wmax = 0):
    """
    When one would want to change the minimum or maximum wavenumber, one would have to know this exactly, which one
    often does not. This function is able to find the wavenumbers within the data closest to the given wavenumbers and
    gives their value and position.

    :param w_data: The wavenumber data
    :param new_wmin: New preferred minimum wavenumber (-1 if no new minimum needs to be found)
    :param new_wmax: New preferred maximum wavenumber (0 if no new maximum needs to be found)
    :return: Gives back the new minimum wavenumber, new maximum wavenumber, position of the new minimum wavenumber and
    the position of the new maximum wavenumber.
    """
    wmin, wmax, wmin_pos, wmax_pos = 0, 0, 0, 0
    wmin_bool, wmax_bool = True, True
    for w in w_data:
        if w <= new_wmin and wmin_bool:
            wmin_bool = False
            wmin_pos = np.where(w_data == w)[0][0]
            if abs(w_data[wmin_pos] - new_wmin) < abs(w_data[wmin_pos-1] - new_wmin):
                wmin = w
            else:
                wmin_pos -= 1
                wmin = w_data[wmin_pos]

        if w <= new_wmax and wmax_bool:
            wmax_bool = False
            wmax_pos = np.where(w_data == w)[0][0]
            if abs(w_data[wmax_pos] - new_wmax) < abs(w_data[wmax_pos-1] - new_wmax):
                wmax = w
            else:
                wmax_pos -= 1

    if new_wmin == -1:
        wmin = w_data[-1]
        wmin_pos = -1

    elif new_wmax == 0:
        wmax = w_data[0]
        wmax_pos = 0

    return wmin, wmax, wmin_pos, wmax_pos

def spectrum_in_air_creator(mol, mf, pres, temp, path_l, wl_min, wl_max, step):
    """
    Function that can be used to simulate a spectrum for a gas. It's not much changed from the original code of radis, but
    it's made like this, so that in the original file, it's slightly cleaner, since some fo the parameters that are
    constant, don't have to be filled in (like the data coming from Hitran for example).

    :param mol: Type of molecule one wants to gain the spectra from.
    :param mf: The mole fraction one wants to simulate.
    :param pres: The pressure at which one wants their gas to be.
    :param temp: The temperature at which one wants their gas to be.
    :param path_l: The path length the light goes through the gas one wants to look at.
    :param wl_min: The minimum wavelength one wants to simulate their spectrum for.
    :param wl_max: The maximum wavelength one wants to simulate their spectrum for.
    :param step: The step size (for the wavelength) one wants to simulate their spectrum for
    :return: Returns the created simulated spectrum.
    """
    spectrum = rd.calc_spectrum(wl_min, wl_max, molecule=mol, isotope="1", pressure=pres, Tgas=temp,
                                    wstep=step, path_length=path_l, databank="hitran", mole_fraction=mf,
                                    medium="air", warnings={"AccuracyError": "ignore"})
    return spectrum

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before

def trold_to_trnew(Told, c_factor):
    """
    Function that transforms an old transmission spectra to a new transmission spectra using a factor c. This function
    is used to fasten the fitting process.

    :param Told: The old transmission spectra
    :param c_factor: The factor with which one wants to change their transmission spectra
    :return: The new transmission spectra
    """
    Tnew = 10 ** (-c_factor * np.log10(1 / Told))
    return Tnew

def tr_to_ab(transmission):
    """
    Function which transforms a transmission spectra to an absorption spectra.

    :param transmission: The transmission spectra that needs to be transformed
    :return: The new absorption spectra
    """
    tr = np.array(transmission)
    return np.log10(1/tr)

def ab_to_tr(absorption):
    """
    Function that transforms an absorption spectra to a transmission spectra.

    :param absorption: The absorption spectra that needs to be transformed
    :return: The new transmission spectra
    """
    ab = np.array(absorption)
    return 10**(-ab)


