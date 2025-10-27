import radis as rd
import os
import os.path as op
import scipy
from brukeropusreader import read_file
from bisect import bisect_left
import csv
from scipy.sparse.linalg import spsolve
from scipy import sparse
import numpy as np
from lmfit.models import ExponentialModel, SkewedGaussianModel


def read_opus_data_from_folder_into_array_for_gui(dir_inv):
    """
    The function gets the usable data-files from the given directory. The usable files are either OPUS-files, csv-files
    or dat-files (not yet implemented as they are not yet used).

    :param dir_inv: given directory where to find the data-files
    :return: array filled with the gained data-files from the given directory
    """
    directory_data = dir_inv
    # Make a list will all files within the folder (for the given directory)
    files_in_directory_temp = [f for f in os.listdir(directory_data) if op.isfile(op.join(directory_data, f))]
    files_in_directory = []
    # Recreate this list with all files within the folder, but now only the files containing the ftir-data
    for file in files_in_directory_temp:
        if file[-2:] == ".0":
            files_in_directory.append(file)
        elif file[-4:] == ".csv":
            files_in_directory.append(file)
        elif file[-4:] == ".dat":
            files_in_directory.append(file)
        elif file[-4:] == ".tsv":
            files_in_directory.append(file)


    # Create array and get the experimental data from the files and add this to the array
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
        elif file[-4:] == ".tsv":
            w_data = []
            t_data = []
            with open(directory_data + "\\" + file) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter="\t")
                line_count = 0
                for row in csv_reader:
                    if row[0] != "X":
                        w_data.append(float(row[0]))
                        t_data.append(float(row[1]))
            full_data = [np.array(w_data), np.array(t_data)]
            dat_array.append(full_data)
    return dat_array


def spectrum_in_air_creator(mol, mf, pres, temp, path_l, wl_min, wl_max, step, dil='air'):
    """
    Function that can be used to simulate a spectrum for a gas. It's not much changed from the original code of radis,
    but it's made like this, so that in the original file, it's slightly cleaner, since some fo the parameters that are
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
    iso = "all"

    # Certain molecules might be special in their choice of isotopes, not choosing all of them:
    # - NO2 is such a special molecule, as the secondary isotope only has known data between 1500 and 1700, so we ignore
    # this isotope. Maybe later, we will try and see if we can also incorporate this isotope.
    if mol == "NO2":
        iso = 1

    # Create needed spectrum
    spectrum = rd.calc_spectrum(wl_min, wl_max, molecule=mol, isotope=iso, pressure=pres, Tgas=temp,
                                wstep=step, path_length=path_l, databank="hitran", mole_fraction=mf, medium="air",
                                warnings={"AccuracyError": "ignore"}, diluent=dil)
    return spectrum


def take_closest(array, value):
    """
    Function that calculates the closest value within an array (or list), doing this for the given value
    :param array: Array or list in which the function needs to search
    :param value: The given value one wants to find the closest value for
    :return: The closest value
    """
    array = np.asarray(array)
    idx = np.nanargmin((np.abs(array - value)))
    return array[idx]


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
    return np.log10(1 / tr)


def ab_to_tr(absorption):
    """
    Function that transforms an absorption spectra to a transmission spectra.

    :param absorption: The absorption spectra that needs to be transformed
    :return: The new transmission spectra
    """
    ab = np.array(absorption)
    return 10 ** (-ab)


def baseline_correction_calculator(y, baseline_diff, lam, p, niter=10):
    """
    Function that calculates the baseline of any given array

    :param y: Array from which one wants to calculate the baseline
    :param baseline_diff: If one wants to know an "average baseline", one can choose to do this for all values. However,
    sometimes we want to ignore peaks or anything like this, which will change the average value a lot. For this
    purpose, we added the option that of giving a "baseline_diff", which is the maximum that a value can differ from the
    baseline, for it to be taken into account in calculating the average
    :param lam: The lam-factor decides how "flat" the baseline is made. If this vaue is low, averaging happens over only
    a few values at a time, so the data is followed quite a lot and thus noise and peaks/valleys can still be slightly
    visible within the baseline. The higher this value, averaging happens over more values and the stronger noise is
    flattened
    :param p: The p-factor decides which part of the data is more "important". If this value is 0.5, the baseline is
    calculated at the exact half of the data. If this value becomes lower, peaks are ignored more and more and only the
    lower values "count". If this value becomes higher, valleys are ignored more and more and only the higher values
    "count".
    :param niter: Amount of iterations before we accept the found baseline. Normal value is 10, can be made higher if
    found solution isn't acceptable
    :return: 1) The baseline calculated for y. 2) The average of the baseline found.
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    D = lam * D.dot(D.transpose())  # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    baseline = np.array([])

    for i_bas in range(niter):
        W.setdiag(w)  # Update diagonal values of the matrix
        Z = W + D
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
        baseline = z

    list_baseline = []
    for i in range(len(baseline)):
        if np.abs(y[i] - baseline[i]) <= baseline_diff:
            list_baseline.append(y[i])
    baseline_average = np.mean(list_baseline)
    return baseline, baseline_average


def fit_baseline(w, a, a_b, bool_mct):
    """
    Function to fit baseline. Important to note, However, due to the problem arising within our machine (see IceBand-Fitting), we also
    added the option to fit an extra skewed gaussian and exponential to fight this problem (depending on the amount of
    time the mct-detector has gotten to warm up, it can lead to a misfitting of around 10 to 1000 ppm of methane).

    :param w: The list with wavenumbers for which the data was gained
    :param a: The list containing the experimental absorption spectra for which a baseline needs to be fitted
    :param a_b: The list containing the experimental absorption spectra which is used as a background
    :param bool_mct: The boolean that decides if the fitted baseline also has to include the mct_detector fit or not.
    :return: baseline: The complete baseline used to take away the background and iceband
    """

    def mct_baseline_creator(x, A1, center, sigma, gamma, A2, decay):
        # Function for the mct-detector fitting.
        if A1 > 0:
            skewed = (A1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-1 * (x - center) ** 2 / 2 / sigma ** 2) * \
                     (1 + scipy.special.erf(gamma * (x - center) / (sigma * np.sqrt(2))))
        else:
            skewed = np.zeros(len(x))

        exp_now = A2 * np.exp(-1 * x / decay)

        full_baseline = skewed + exp_now
        return full_baseline

    # Create the baseline for the background file
    bas_b, b_b = baseline_correction_calculator(a_b, 0.0001, 1000000, 0.5)

    # Remove background baseline from experimental absorption
    a_background_removed = a - bas_b

    # Create baseline for background_removed (used for fitting of mct-detector)
    bas_background_removed, b_background_removed = baseline_correction_calculator(
        a_background_removed, 0.0001, 1000000, 0.5)

    # We only use values which are close enough to the baseline, so that no peaks/valleys are taken into account during
    # the fitting of the mct-detector
    bas_line_temp = []
    bas_w = []
    for k in range(len(w)):
        if (np.abs(a_background_removed[k] - bas_background_removed[k])) <= 0.00005:
            bas_w.append(w[k])
            bas_line_temp.append(a_background_removed[k])

    bas_w = np.array(bas_w)
    bas_line_temp = np.array(bas_line_temp)

    # Using lmfit to fit a skewed gaussian to the iceband and an exponential to the heating of the mct-detector
    if bool_mct:
        gaus = SkewedGaussianModel(prefix="g_")
        g_center_min = 3134
        g_center_max = 3139
        pars = gaus.guess(bas_line_temp, x=bas_w)
        pars["g_center"].set(value=3136, min=g_center_min, max=g_center_max)
        pars["g_sigma"].set(value=155.3, min=162, max=166.5)
        pars["g_gamma"].set(value=2.33, min=1.95, max=2.8)

        exp = ExponentialModel(prefix="e_")
        pars += exp.guess(bas_line_temp, x=bas_w)
        model = gaus + exp


        out = model.fit(bas_line_temp, pars, x=bas_w)

        g_center = out.best_values["g_center"]
        g_amplitude = out.best_values["g_amplitude"]
        g_sigma = out.best_values["g_sigma"]
        g_gamma = out.best_values['g_gamma']
        e_amplitude = out.best_values["e_amplitude"]
        e_decay = out.best_values["e_decay"]

        # Create final baseline
        baseline = mct_baseline_creator(x=w, A1=g_amplitude, center=g_center, sigma=g_sigma, gamma=g_gamma, A2=e_amplitude,
                                    decay=e_decay)
    else:
        baseline = bas_b

    return baseline
