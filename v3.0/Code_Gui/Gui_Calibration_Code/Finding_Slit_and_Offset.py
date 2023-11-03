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
from Code_Gui.Gui_General_Code import General_Functions_Library as GFL
import numpy as np
from radis import Spectrum
import matplotlib.pyplot as plt
import Fitting_Ice_Band as FIB
from lmfit import Model

def spectra_linear_wavenumber_dependent_offset_3_molecules(w, c_ch4, c_h2o, c_co2, offset_left, offset_right, slit_size):
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
    c_list = [c_ch4, c_h2o, c_co2]
    t_dict_new = {}
    a_dict_new = {}
    a_full = np.array([])

    for i in range(len(dict_w3)):
        t_dict_new[list(dict_w3.keys())[i]] = GFL.trold_to_trnew(dict_t3[list(dict_w3.keys())[i]], c_list[i])
        a_dict_new[list(dict_w3.keys())[i]] = GFL.tr_to_ab(t_dict_new[list(dict_w3.keys())[i]])
        if i == 0:
            a_full = a_dict_new[list(dict_w3.keys())[i]]
        else:
            a_full += a_dict_new[list(dict_w3.keys())[i]]

    t_full = GFL.ab_to_tr(a_full)
    molecule = list(dict_w3.keys())[0]

    spec_new = Spectrum({"wavenumber": dict_w3[molecule], "transmittance_noslit": t_full}, wunit='cm-1',
                        units={"transmittance_noslit": ""},
                        conditions={"path_length": pathlength_0, "self_absorption": None})
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")

    w_temp, t_temp = spec_new.get("transmittance")

    offset_a = (offset_right - offset_left)/(4150-2200)
    offset_b = offset_right-offset_a*4150

    list_offset_temp = offset_a*w_temp + offset_b

    w_new_temp = np.zeros(len(w_temp))
    for i in range(len(w_temp)):
        w_new_temp[i] = w_temp[i] + list_offset_temp[i]

    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength_0})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")


    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new

directory = "C:\\Users\\P70085588\\Data\\Invenio-r\\2023\\07\\17\\Fit_Params_Test\\"

files_in_directory_temp = [f for f in os.listdir(directory) if op.isfile(op.join(directory, f))]
array_data = GFL.read_opus_data_from_folder_into_array_for_gui(directory)
c_ch4_list = []
c_co2_list = []
c_h2o_list = []
pathlength_list = []
list_offset_left = []
list_offset_right = []
list_slitsize = []

for i in range(0, len(files_in_directory_temp)-3):
    data_text_file_selected = directory + files_in_directory_temp[i]
    print(files_in_directory_temp[i])
    data = array_data[i]

    if os.path.exists(data_text_file_selected) and os.stat(data_text_file_selected).st_size > 0:
        t_list = [295.55, 295.55, 296.65,296.75,296.85,296.85,295.95,295.85,295.85,295.85,295.85,295.85,295.85,295.85]
        p_list = [1.041, 1.0357, 1.0357, 1.0357, 1.0357, 1.0286, 1.0286, 1.0392, 1.0392, 1.0353, 1.0353, 1.041, 1.041, 1.041]
        c_list = [0,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,0,0,0]
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

            for j in range(len(w_background)):
                if w_background[j] >= 1000:
                    wlmin=w_background[j]
                    wlmin_j = j
                if w_background[j] >= 4150:
                    wlmax=w_background[j]
                    wlmax_j = j

            w_background = w_background[wlmax_j: wlmin_j]
            t_background = t_background[wlmax_j: wlmin_j]
            a_background = GFL.tr_to_ab(t_background)

        else:
            s_exp = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                 units={"transmittance": ""})
            w_exp, t_exp = s_exp.get("transmittance")

            w_exp = w_exp[wlmax_j: wlmin_j]
            t_exp = t_exp[wlmax_j: wlmin_j]
            a_exp = GFL.tr_to_ab(t_exp)

            bas = FIB.fit_baseline(w_exp, a_exp, a_background,True)
            a_exp_baseline_corrected = a_exp - bas
            t_exp_baseline_corrected = GFL.ab_to_tr(a_exp_baseline_corrected)

            for j in range(len(w_background)):
                if w_background[j] >= 2200:
                    wlmin2=w_background[j]
                    wlmin_j2 = j
                if w_background[j] >= 4150:
                    wlmax2=w_background[j]
                    wlmax_j2 = j

            w_exp = w_exp[wlmax_j2: wlmin_j2]
            t_exp_baseline_corrected = t_exp_baseline_corrected[wlmax_j2:wlmin_j2]

            dict_molecules3 = {}
            dict_molecules3["CO2"] = 0.001
            dict_molecules3["H2O"] = 0.001
            dict_molecules3["CH4"] = 0.001
            list_molecules3 = ["CH4", "H2O", "CO2"]

            for j in range(len(list_molecules3)):
                molecule = list_molecules3[j]
                spec = GFL.spectrum_in_air_creator(mol=molecule, mf=c_list[i]/1000000, pres=p_list[i], temp=t_list[i],
                                                       path_l=pathlength_0, wl_min=wlmin2, wl_max=wlmax2, step=0.002)
                dict_w3[molecule], dict_t3[molecule] = spec.get("transmittance_noslit")

            model3 = Model(spectra_linear_wavenumber_dependent_offset_3_molecules)
            params3 = model3.make_params(c_ch4=1, c_h2o=1, c_co2=1, offset_left=-0.075, offset_right=-0.14,
                                         slit_size=0.3)
            params3['c_ch4'].set(min=0.96, max=1.2)
            params3['c_h2o'].set(min=0.00001, max=100)
            params3["c_co2"].set(min=0.00001, max=100)
            params3['offset_left'].set(min=-0.2, max=0)
            params3['offset_right'].set(min=-0.2, max=0)
            params3['slit_size'].set(min=0.29, max=0.35)
            result3 = model3.fit(t_exp_baseline_corrected, params3, w=w_exp)
            print(result3.fit_report())

            plt.plot(w_exp, t_exp_baseline_corrected - result3.best_fit, alpha=0.7,  label=str(i))

            c_h2o_list.append(result3.best_values["c_h2o"])
            c_ch4_list.append(result3.best_values["c_ch4"])
            c_co2_list.append(result3.best_values["c_co2"])
            list_offset_left.append(result3.best_values["offset_left"])
            list_offset_right.append(result3.best_values["offset_right"])
            list_slitsize.append(result3.best_values["slit_size"])


            print(c_ch4_list)
            print(c_co2_list)
            print(c_h2o_list)
            print(list_offset_left)
            print(list_offset_right)
            print(list_slitsize)

plt.legend()
plt.show()