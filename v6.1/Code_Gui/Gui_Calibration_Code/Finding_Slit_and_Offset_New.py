"""
This document is used to find the following parameters:
- slit size = 0.30926986906880655
- offset left (at 2200 cm-1) = -0.07197559987285236
- offset right (at 4150 cm-1) = -0.14074011578087747

Two things to note:
- Technically, one can also use this to find the correct pathlength. However, for this, perfect fitting data is needed,
without any air dillution. In a later stage, I might get back to this.
- Though the offset is computed, one can choose to have this remain a fitting parameter. This is necessary if one wants
to fit in a different range than from 2200 to 4150 (as the offset boundaries are dependent on the wavenumber-boundaries)
"""



import os
import os.path as op

from scipy.signal import *
import scipy.integrate as inte
from Code_Gui.Gui_General_Code import General_Functions_Library as GFL
import numpy as np
from radis import Spectrum
import matplotlib.pyplot as plt
from lmfit import Model
from lmfit import minimize, Parameters, fit_report
from ttictoc import tic,toc
from matplotlib.gridspec import GridSpec
import wedme
from scipy.special import voigt_profile
wedme.dev()

def spectral_fit_residuals(pars: Parameters, w_exp_short, t_exp_short, test1=False, test2 = False):
    """
       This fuction is used in order to fit a calculated spectra to the gained experimental data. This formula can only be
       used when one wants to fit a single molucele. Because of the way scipy.curve_fit() works, this function is repeated
       below, but with more c's, so more molecules can be fitted at the same time.

       :param w: wavenumber range needed to match simulated with experimental data
       :param offset_min: Minimimal offset-error compared with simulated data
       :param offset_max: Maximal offset-error compared with simulated data
       :param slit_size: Size of the needed slit to match simulated with experimental data
       :return: return fitted transmittance
       """
    tic()
    c_list = [pars['c_ch4'].value, pars['c_h2o'].value, pars['c_co2'].value]
    w_g = pars['w_g'].value
    w_l = pars['w_l']
    k0 = pars['k0'].value
    k1 = pars['k1'].value

    w_short = w_exp_short
    t_short = t_exp_short


    t_dict_new = {}
    a_dict_new = {}
    a_full = np.array([])

    for i in range(len(dict_w3)):
        w_full = dict_w3[list(dict_w3.keys())[i]]
        t_dict_new[list(dict_w3.keys())[i]] = GFL.trold_to_trnew(dict_t3[list(dict_w3.keys())[i]], c_list[i])
        a_dict_new[list(dict_w3.keys())[i]] = GFL.tr_to_ab(t_dict_new[list(dict_w3.keys())[i]])
        if i == 0:
            a_full = a_dict_new[list(dict_w3.keys())[i]]
        else:
            a_full += a_dict_new[list(dict_w3.keys())[i]]

    #t_full = GFL.ab_to_tr(a_full)

    #spec_temp = Spectrum({"wavenumber": w_full, "transmittance_noslit": t_full}, wunit='cm-1',
    #                    units={"transmittance_noslit": ""},
    #                    conditions={"path_length": pathlength_0, "self_absorption": None})

    av_w = np.mean(w_full)
    w_full_voigt = w_full - av_w
    voigt_full = voigt_profile(w_full_voigt, w_g, w_l)

    a_temp = convolve(a_full, voigt_full, mode="same")
    i_temp = inte.simpson(a_temp, x=w_full)
    i_full = inte.simpson(a_full, x=w_full)
    a_temp *= (i_full/i_temp)
    t_temp = GFL.ab_to_tr(a_temp)
    w_full_temp = k0 + k1*w_full

    #a_exp = GFL.tr_to_ab(t_exp_short)
    #plt.plot(w_exp_short, a_exp, alpha=0.7)
    #plt.plot(w_full, a_full, alpha=0.7)
    #plt.plot(w_full, a_temp, alpha=0.7)
    #plt.show()

    spec_new = Spectrum({"wavenumber": w_full_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength_0})

    spec_new.resample(w_short, inplace=True)
    w_new, t_new = spec_new.get("transmittance")


    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1

    if test1 is True:
        plt.plot(w_exp_short, t_short, alpha=0.7)
        plt.plot(w_new, t_new, alpha=0.7)
        plt.show()
        print(f"Execution time: {toc() * 1e3:.0f}ms")
        return (t_new - t_short) ** 2
    elif test2 is True:
        return t_new

    else:
        print(f"Execution time: {toc() * 1e3:.0f}ms")
        return (t_new - t_short) ** 2

directory = "C:\\Users\\P70085588\\Data\\FTIR\\Invenio-r\\2025\\8\\6_2\\"

files_in_directory_temp = [f for f in os.listdir(directory) if op.isfile(op.join(directory, f))]
print(files_in_directory_temp)
array_data = GFL.read_opus_data_from_folder_into_array_for_gui(directory)
c_ch4_list = []
c_co2_list = []
c_h2o_list = []
pathlength_list = []
list_k0 = []
list_k1 = []
list_w_g = []
list_w_l = []

for i in range(0, len(files_in_directory_temp)):
    data_text_file_selected = directory + files_in_directory_temp[i]
    print(files_in_directory_temp[i])
    data = array_data[i]

    if os.path.exists(data_text_file_selected) and os.stat(data_text_file_selected).st_size > 0:
        t_list = [0,298.15,298.25,298.25,298.25,298.25,298.25]
        p_list = [0,0.0975,0.0977,0.0984,0.0984,0.0975,0.0975]
        c_list = [0,0.002,0.002,0.0012,0.0012,0.0002,0.0002]
        pathlength_0 = 20.062

        w_exp = []
        t_exp = []
        t_fit = []
        t_residual = []
        dict_w2 = {}
        dict_t2 = {}
        dict_w3 = {}
        dict_t3 = {}
        try:
            data_wavelength = np.array(data.get_range("AB")[0:-1])
            data_transmission = np.array(data["AB"][0:-1])
        except:
            data_wavelength = data[0]
            data_transmission = data[1]

        if i==0:
            s_background = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                             units={"transmittance": ""})
            w_background, t_background = s_background.get("transmittance")
            wlmin = min(data_wavelength)
            wlmax= max(data_wavelength)
            wlmax_j = 0
            wlmin_j= len(w_background)
            a_background = GFL.tr_to_ab(t_background)


        else:
            s_exp = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                 units={"transmittance": ""})
            w_exp, t_exp = s_exp.get("transmittance")

            a_exp = GFL.tr_to_ab(t_exp)

            bas = GFL.fit_baseline(w_exp, a_exp, a_background,True)
            a_exp_baseline_corrected = a_exp - bas
            t_exp_baseline_corrected = GFL.ab_to_tr(a_exp_baseline_corrected)

            for j in range(len(w_background)):
                if w_background[j] <= 1100:
                    wlmin2=w_background[j]
                    wlmin_j2 = j
                if w_background[j] <= 4000:
                    wlmax2=w_background[j]
                    wlmax_j2 = j

            w_exp = w_exp[wlmin_j2:wlmax_j2]
            t_exp_baseline_corrected= t_exp_baseline_corrected[wlmin_j2:wlmax_j2]

            dict_molecules3 = {}
            dict_molecules3["CO2"] = 0.0002
            dict_molecules3["H2O"] = 0.0003
            dict_molecules3["CH4"] = c_list[i]
            list_molecules3 = ["CH4", "H2O", "CO2"]

            for j in range(len(list_molecules3)):
                molecule = list_molecules3[j]
                spec = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules3[molecule], pres=p_list[i], temp=t_list[i],
                                                       path_l=pathlength_0, wl_min=w_exp[0], wl_max=w_exp[-1], step=0.002)
                dict_w3[molecule], dict_t3[molecule] = spec.get("transmittance_noslit")

            fitParameters = Parameters()
            fitParameters.add("c_ch4", value=1, min=0.95, max=1.05, vary=False)
            fitParameters.add("c_h2o", value=1, min=0, max=1000, vary=True)
            fitParameters.add("c_co2", value=1, min=0.00, max=1000, vary=True)
            fitParameters.add("k0", value=0, min=-10, max=10, vary=True)
            fitParameters.add("k1", value=1, min=0.95, max=1.05, vary=True)
            fitParameters.add("w_g", value=0.02, min=0.0001, max=1, vary=True)
            fitParameters.add("w_l", value=0.0, min=0.000, max=1, vary=True)


            spectral_fit_residuals(fitParameters, w_exp, t_exp_baseline_corrected,test1=True)
            out = minimize(spectral_fit_residuals, fitParameters,
                           args=(w_exp, t_exp_baseline_corrected), method='leastq')
            fit = spectral_fit_residuals(out.params, w_exp,t_exp_baseline_corrected, test2=True)
            res = t_exp_baseline_corrected - fit
            print(fit_report(out))

            c_list_old = {}
            c_list_minus = {}

            c_list_old["CH4"] = out.params['c_ch4'] * dict_molecules3["CH4"]
            c_list_old["H2O"] = out.params['c_h2o'] * dict_molecules3["H2O"]
            c_list_old["CO2"] = out.params['c_co2'] * dict_molecules3["CO2"]


            c_h2o_list.append(c_list_old["H2O"]*1000000)
            c_ch4_list.append(c_list_old["CH4"]*1000000)
            c_co2_list.append( c_list_old["CO2"]*1000000)
            list_k0.append(out.params["k0"].value)
            list_k1.append(out.params["k1"].value)
            list_w_g.append(out.params["w_g"].value)
            list_w_l.append(out.params["w_l"].value)

            print(c_ch4_list)
            print(c_co2_list)
            print(c_h2o_list)
            print(list_k0)
            print(list_k1)
            print(list_w_g)
            print(list_w_l)

            fig = plt.figure()
            gs = GridSpec(nrows=2, ncols=1, height_ratios=[4, 1])
            ax0 = fig.add_subplot(gs[0, 0])
            ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)


            ax0.plot(w_exp, t_exp_baseline_corrected, alpha=0.7)
            ax0.plot(w_exp, fit, alpha=0.7)
            #ax0.plot(w_exp, result_refit.best_fit, alpha=0.7)
            plt.tick_params('x', labelbottom=False)
            ax0.legend(loc='upper right', fontsize=12)
            ax0.set_ylabel('Intensity (a.u.)')
            ax0.grid(True)

            ax1.plot(w_exp, res)
            ax1.set_xlabel("wavenumber ($cm^{-1}$)")
            ax1.set_ylabel('Residual')
            ax1.grid(True)
            fig.tight_layout()

            plt.show()
