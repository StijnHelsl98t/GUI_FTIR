from radis import Spectrum
import numpy as np
import Code_Gui.Gui_General_Code.General_Functions_Library as GFL
from lmfit import minimize, Parameters, fit_report

w_dict = {}
t_dict = {}
pathlength = 0
slit_size = 0
k_0 = 0
k_1 = 0


def storage_for_dict(temp, pres, pa_le,sl_si,k_1):
    """
    This function creates a storage place, so that the wavenumber and transmission gained for the spectra can be
    imported into this module. This was the only way we were be able to import all the needed parameters correctly,
    so that the fit could be taken.
    :param w: The imported dictionary containing the wavenumber-arrays
    :param t: The imported dictionary containing the transmission-arrays
    :param p: The imported pathlength, measured value for the specific gas-cell used
    :param s: The imported slit_size, calculated within calibration for the specific machine used
    :param k0: The imported k0, which is the offset from 0
    :param k1: The imported k1, which is a linear change in the x-axis
    :return: -
    """
    global temperature, pressure, pathlength, slit_size, k1
    
    temperature = temp
    pressure=pres
    pathlength = pa_le
    slit_size = sl_si
    k1 = k_1

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
    k0 = pars['k0']
    c_list = []
    for i in range(len(pars)-1):
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
                        units={"transmittance_noslit": ""}, conditions={"path_length": pathlength},)

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
    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    # if any non-numbers exist within the spectra, change these into a no-molecule zone (transmission = 1)
    if test:
        return t_new
    else:
        return (t_new - t_meas)**2

def spectra_molecules_s(pars: Parameters, w_meas, t_meas, s, test=False):
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
    k0 = pars['k0']
    dict_spec_new = {}
    dict_c = {}
    print(s.keys())
    for molecule in s.keys():
        dict_c[molecule] = pars['c_' + molecule].value 
        dict_spec_new[molecule] = s[molecule].rescale_mole_fraction(dict_c[molecule])
        try:
            spec_new = spec_new // dict_spec_new[molecule]
        except:
            spec_new = dict_spec_new[molecule]
    
    # Apply experimental slit (broadening coefficient term) to spectra object
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")
    # Get necessary list with wavenumbers and transmission
    w_temp, t_temp = spec_new.get("transmittance")
    w_new_temp = k0 + k1*w_temp

    # Re-create spectrum objectw
    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

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


    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4, c5, c6, c7, c8]
    k_0 = k0
    t_dict_new = {}
    a_dict_new = {}
    a_full = {}
    for i in range(len(w_dict)):
        t_dict_new[list(w_dict.keys())[i]] = GFL.trold_to_trnew(t_dict[list(w_dict.keys())[i]], c_list[i])
        a_dict_new[list(w_dict.keys())[i]] = GFL.tr_to_ab(t_dict_new[list(w_dict.keys())[i]])
        if i == 0:
            a_full = a_dict_new[list(w_dict.keys())[i]]
        else:
            a_full += a_dict_new[list(w_dict.keys())[i]]

    t_full = GFL.ab_to_tr(a_full)
    molecule = list(w_dict.keys())[0]
    spec_new = Spectrum({"wavenumber": w_dict[molecule], "transmittance_noslit": t_full}, wunit='cm-1',
                        units={"transmittance_noslit": ""}, conditions={"path_length": pathlength})
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True)

    w_temp, t_temp = spec_new.get("transmittance")
    w_new_temp = k_0 + k_1*w_temp
    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")

    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new