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
from radis import Spectrum, MergeSlabs, load_spec
import matplotlib.pyplot as plt
from lmfit import Model
from lmfit import minimize, Parameters, fit_report
from matplotlib.gridspec import GridSpec

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

    offset_a = (offset_right - offset_left)/(4000-1000)
    offset_b = offset_right-offset_a*4000

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
def spectra_molecules(pars: Parameters, w_meas, t_meas, dict_w, dict_t, test=False):
    """
    This fuction is used in order to fit a calculated spectra to the gained experimental data. This formula can only be
    used when one wants to fit a single molecule. Because of the way scipy.curve_fit() works, this function is repeated
    below, but with more c's added, so more molecules can be fitted at the same time. Also, due to the way
    scipy.curve_fit() works, we make use of a storage-formula.

    :param w: wavenumber range needed to match simulated with experimental data
    :param c1: mole fraction of the molecule one wants to fit
    :param slit_size: Size of the needed slit to match simulated with experimental data
    :param k0: The imported k0, which is the offset from 0
    :param k1: The imported k1, which is a linear change in the x-axis
    :return: return fitted transmittance
    """
    k0 = pars['k0'].value
    k1 = pars['k1'].value
    slit_size = pars['slit'].value
    c_list = []
    for i in range(len(pars)-3):
        c_list.append(pars['c' + str(i+1)])

    dict_t_new = {}
    dict_a_new = {}
    list_a_full = {}
    for i in range(len(dict_w)):
        # Convert all the old transmission spectra to new transmission spectra, using a factor change of c
        dict_t_new[list(dict_w.keys())[i]] = GFL.trold_to_trnew(dict_t[list(dict_w.keys())[i]], c_list[i])
        # Convert all the transmission spectra to absorption spectra, so that they can be added together
        dict_a_new[list(dict_w.keys())[i]] = GFL.tr_to_ab(dict_t_new[list(dict_w.keys())[i]])
        if i == 0:
            list_a_full = dict_a_new[list(dict_w.keys())[i]]
        else:
            list_a_full += dict_a_new[list(dict_w.keys())[i]]
    # Convert the full absorption spectra to a transmission spectra
    list_t_full = GFL.ab_to_tr(list_a_full)

    molecule = list(dict_w.keys())[0]
    
    w_new_temp = k0 + k1*dict_w[molecule]

    # Use the created transmission spectra to create the needed spectra object (for RADIS)
    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance_noslit": list_t_full}, wunit='cm-1',
                        units={"transmittance_noslit": ""}, conditions={"path_length": pathlength_0},)

    # Apply experimental slit (broadening coefficient term) to spectra object
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")
    # Get necessary list with wavenumbers and transmission
    #w_temp, t_temp = spec_new.get("transmittance")
    #w_new_temp = k_0 + k_1*w_temp

    # Re-create spectrum objectw
    #spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
    #                    units={"transmittance": ""}, conditions={"path_length": pathlength})

    # Resample spectrum object onto experimental spectrum, to be able to compare them
    spec_new.resample(w_meas, inplace=True)
    # Get necessary list with wavenumbers and transmission
    w_new, t_new = spec_new.get("transmittance")
    t_new[np.isnan(t_new)]=1
    t_new[t_new==0] = 1*10**-4
    # if any non-numbers exist within the spectra, change these into a no-molecule zone (transmission = 1)
    if test:
        return t_new
    else:
        return (t_new - t_meas)**2
def spectra_molecules2(pars: Parameters, w_meas, t_meas, s, test=False):
    """
    This fuction is used in order to fit a calculated spectra to the gained experimental data. This formula can only be
    used when one wants to fit a single molecule. Because of the way scipy.curve_fit() works, this function is repeated
    below, but with more c's added, so more molecules can be fitted at the same time. Also, due to the way
    scipy.curve_fit() works, we make use of a storage-formula.

    :param w: wavenumber range needed to match simulated with experimental data
    :param c1: mole fraction of the molecule one wants to fit
    :param slit_size: Size of the needed slit to match simulated with experimental data
    :param k0: The imported k0, which is the offset from 0
    :param k1: The imported k1, which is a linear change in the x-axis
    :return: return fitted transmittance
    """
    k0 = pars['k0'].value
    k1 = pars['k1'].value
    slit_size = pars['slit'].value
    dict_spec_new = {}
    dict_c = {}
    spectra = []
    for molecule in s.keys():
        dict_c[molecule] = pars['c_' + molecule].value 
        spectra.append(s[molecule].rescale_mole_fraction(dict_c[molecule]))

    print(dict_c)

    spec_new = MergeSlabs(*spectra, out="transparent")
    # Apply experimental slit (broadening coefficient term) to spectra object
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")
    # Get necessary list with wavenumbers and transmission
    w_temp, t_temp = spec_new.get("transmittance")
    w_new_temp = k0 + k1*w_temp

    # Re-create spectrum objectw
    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength_0})

    # Resample spectrum object onto experimental spectrum, to be able to compare them
    spec_new.resample(w_meas, inplace=True)
    # Get necessary list with wavenumbers and transmission
    w_new, t_new = spec_new.get("transmittance")
    t_new[np.isnan(t_new)]=1
    t_new[t_new==0] = 1*10**-4
    # if any non-numbers exist within the spectra, change these into a no-molecule zone (transmission = 1)
    if test:
        return t_new
    else:
        return (t_new - t_meas)**2

directory = "C:/Users/P70085588/Data/Invenio-R/2026/03/17/Methane_Check/"

files_in_directory_temp = [f for f in os.listdir(directory) if op.isfile(op.join(directory, f))]
print(files_in_directory_temp)
array_data = GFL.read_opus_data_from_folder_into_array_for_gui(directory)
c_ch4_list = []
pathlength_list = []
list_k0 = []
list_k1 = []
list_slitsize = []

for i in range(0, len(files_in_directory_temp)):
    data_text_file_selected = directory + files_in_directory_temp[i]
    print(files_in_directory_temp[i])
    data = array_data[i]

    if os.path.exists(data_text_file_selected) and os.stat(data_text_file_selected).st_size > 0:
        t_list = [0,298.15,298.45,298.15,297.55,297.55,297.55, 297.55]
        p_list = [0,0.052,0.0502,0.0506, 0.0506, 0.0506, 0.0506]
        c_list = [0,1,1,0.5,0.01,0.002,0.0005]
        c_list2 = [0,0,0,0,0,0,0,0]
        pathlength_0 = 20.062

        w_exp = []
        t_exp = []
        t_fit = []
        t_residual = []
        dict_w2 = {}
        dict_t2 = {}
        dict_w3 = {}
        dict_t3 = {}
        list_molecules = ["CH4", "CO"]
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
                if w_background[j] <= 1000:
                    wlmin=w_background[j]
                    wlmin_j = j
                    print(wlmin_j)
                if w_background[j] <= 3025:
                    wlmax=w_background[j]
                    wlmax_j = j


            w_background = w_background[wlmin_j: wlmax_j]
            t_background = t_background[wlmin_j: wlmax_j]
            a_background = GFL.tr_to_ab(t_background)


        else:
            s_exp = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                                 units={"transmittance": ""})
            w_exp, t_exp = s_exp.get("transmittance")

            w_exp = w_exp[wlmin_j: wlmax_j]
            t_exp = t_exp[wlmin_j: wlmax_j]

            for t_i in range(len(t_exp)):
                if np.isnan(t_exp[t_i]):
                    t_exp[t_i]=0
                    print(0, True)
            a_exp = GFL.tr_to_ab(t_exp)
            t_exp_min = min(t_exp)
            a_exp_max = max(a_exp)
            
            bas = GFL.fit_baseline(w_exp, a_exp, a_background,False)
            a_exp_baseline_corrected = a_exp - bas
            t_exp_baseline_corrected = GFL.ab_to_tr(a_exp_baseline_corrected)

            for j in range(len(w_background)):
                if w_background[j] <= 1000:
                    wlmin2=w_background[j]
                    wlmin_j2 = j
                if w_background[j] <= 3025:
                    wlmax2=w_background[j]
                    wlmax_j2 = j

            w_exp = w_exp[wlmin_j2:wlmax_j2]
            a_exp_baseline_corrected = a_exp_baseline_corrected[wlmin_j2:wlmax_j2]
            t_exp_baseline_corrected = t_exp_baseline_corrected[wlmin_j2:wlmax_j2]



            dict_molecules = {}
            for mol_i in range(len(list_molecules)):
                dict_molecules[list_molecules[mol_i]] = c_list[i]
                if mol_i == "CO":
                    dict_molecules[list_molecules[mol_i]] = 0
            print(dict_molecules)
            dict_spec = {}
            dict_w = {}
            dict_t = {}
            fitParameters = Parameters()
            fitParameters.add('k0', value=0, min=-2, max=2, vary=True)
            fitParameters.add('k1', value=1, min=0.99, max=1.01, vary=True)
            fitParameters.add('slit', value=0.25, min=0.1, max=0.5, vary=True)

            
            for molecule in dict_molecules.keys(): 
                dict_spec[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=dict_molecules[molecule],
                                                                                    pres=p_list[i], temp=t_list[i],
                                                                                    path_l=pathlength_0,
                                                                                    wl_min=wlmin2, wl_max=wlmax2, step=0.001)
                fitParameters.add('c_' + molecule, value=c_list[i], min=1*10**-6, max=1, vary=False)
            

            out = minimize(spectra_molecules2, fitParameters, args=(w_exp, t_exp_baseline_corrected, dict_spec), method='leastq')
            print(fit_report(out))
            fit = spectra_molecules2(out.params, w_exp, t_exp_baseline_corrected, dict_spec, test=True)

            c_list_old = {}
            c_list_minus = {}

            for molecule in dict_molecules.keys():
                c_list_old[molecule] = out.params['c_' + molecule].value

            refit_bool = True
            for fit_i in range(5):
                print(c_list_old)
                if refit_bool:
                    dict_spec_refit={}
                    dict_w_refit = {}
                    dict_t_refit = {}
                    
                    fitParameters_new = Parameters()
                    if out.params['k0'] >= 0:
                        fitParameters_new.add('k0', value=out.params['k0'].value, min=0.999*out.params['k0'].value, max=1.001*out.params['k0'].value, vary=True)
                    else:
                        fitParameters_new.add('k0', value=out.params['k0'].value, min=1.001*out.params['k0'].value, max=0.999*out.params['k0'].value, vary=True)
                    fitParameters_new.add('k1', value=out.params['k1'].value, min=0.999*out.params['k1'].value, max=1.001*out.params['k1'].value, vary=True)
                    fitParameters_new.add('slit', value=out.params['slit'].value, min=0.999*out.params['slit'].value, max=1.001*out.params['slit'].value, vary=True)
                    
                    
                    for molecule in dict_molecules.keys():
                        dict_spec_refit[molecule] = GFL.spectrum_in_air_creator(mol=molecule, mf=c_list_old[molecule],
                                                                                    pres=p_list[i], temp=t_list[i],
                                                                                    path_l=pathlength_0,
                                                                                    wl_min=wlmin2, wl_max=wlmax2, step=0.001)
                        fitParameters_new.add('c_' + molecule, value=c_list_old[molecule], min=0.9*c_list_old[molecule], max=1, vary=True)
                        #dict_w_refit[molecule] = \
                        #    spectrum_dictionary_refit[molecule].get("transmittance_noslit")[0]
                        #dict_t_refit[molecule] = \
                        #    spectrum_dictionary_refit[molecule].get("transmittance_noslit")[1]
                        #dict_t_refit[molecule][dict_t_refit[molecule]==0] = 1*10**-4

                    print(fitParameters_new)
                            
                            
                    out_refit = minimize(spectra_molecules2, fitParameters_new, args=(w_exp, t_exp_baseline_corrected, dict_spec_refit), method='leastq')
                    print(fit_report(out_refit))
                    refit = spectra_molecules2(out_refit.params, w_exp, t_exp_baseline_corrected, dict_spec_refit, test=True)
                    
                    c_list_new = {}
                    for molecule in dict_molecules.keys():
                        c_list_new[molecule] = out_refit.params['c_' + molecule].value

                    correctness_count = 0
                    for molecule in dict_molecules.keys():
                        print(np.abs(
                            c_list_new[molecule] - c_list_old[molecule]) /
                              c_list_old[molecule])
                        if np.abs(
                                c_list_new[molecule] - c_list_old[molecule]) /c_list_old[molecule] > 0.002:
                            correctness_count += 1
                            print("need refit")
                            c_list_old[molecule] = c_list_new[molecule]
                    if correctness_count == 0:
                        refit_bool = False
                        print("its ok")

            c_ch4_list.append(c_list_old["CH4"])
            list_k0.append(out_refit.params["k0"].value)
            list_k1.append(out_refit.params["k1"].value)
            list_slitsize.append(out_refit.params["slit"].value)
            residual = t_exp_baseline_corrected - refit
            print(c_ch4_list)
            print(list_k0)
            print(list_k1)
            print(list_slitsize)

            fig = plt.figure()
            gs = GridSpec(nrows=2, ncols=1, height_ratios=[4, 1])
            ax0 = fig.add_subplot(gs[0, 0])
            ax1 = fig.add_subplot(gs[1, 0], sharex=ax0) 

            ax0.plot(w_exp, t_exp_baseline_corrected, alpha=0.7, label="exp")
            ax0.plot(w_exp, refit, alpha=0.7, label="fit")
            ax1.plot(w_exp, residual, alpha=0.7)
            plt.legend()
            plt.show()
