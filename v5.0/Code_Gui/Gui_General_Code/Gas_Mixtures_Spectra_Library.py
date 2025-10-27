from radis import Spectrum
import numpy as np
import Code_Gui.Gui_General_Code.General_Functions_Library as GFL

w_dict = {}
t_dict = {}
pathlength = 0
slit_size = 0
k_0 = 0
k_1 = 0


def storage_for_dict(w, t, p, s, k0, k1):
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
    global w_dict, t_dict, pathlength, slit_size, k_0, k_1
    w_dict = w
    t_dict = t
    pathlength = p
    slit_size = s
    k_0 = k0
    k_1 = k1


def spectra_fit_1_molecule(w,k0, c1):
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
    global w_dict, t_dict
    c_list=[c1]
    k_0 = k0
    t_dict_new = {}
    a_dict_new = {}
    a_full = {}
    for i in range(len(w_dict)):
        # Convert all the old transmission spectra to new transmission spectra, using a factor change of c
        t_dict_new[list(w_dict.keys())[i]] = GFL.trold_to_trnew(t_dict[list(w_dict.keys())[i]], c_list[i])
        # Convert all the transmission spectra to absorption spectra, so that they can be added together
        a_dict_new[list(w_dict.keys())[i]] = GFL.tr_to_ab(t_dict_new[list(w_dict.keys())[i]])
        if i == 0:
            a_full = a_dict_new[list(w_dict.keys())[i]]
        else:
            a_full += a_dict_new[list(w_dict.keys())[i]]
    # Convert the full absorption spectra to a transmission spectra
    t_full = GFL.ab_to_tr(a_full)
    molecule = list(w_dict.keys())[0]

    # Use the created transmission spectra to create the needed spectra object (for RADIS)
    spec_new = Spectrum({"wavenumber": w_dict[molecule], "transmittance_noslit": t_full}, wunit='cm-1',
                        units={"transmittance_noslit": ""}, conditions={"path_length": pathlength},)

    # Apply experimental slit (broadening coefficient term) to spectra object
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")
    # Get necessary list with wavenumbers and transmission
    w_temp, t_temp = spec_new.get("transmittance")
    w_new_temp = k_0 + k_1*w_temp

    # Re-create spectrum object
    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    # Resample spectrum object onto experimental spectrum, to be able to compare them
    spec_new.resample(w, inplace=True)
    # Get necessary list with wavenumbers and transmission
    w_new, t_new = spec_new.get("transmittance")
    # if any non-numbers exist within the spectra, change these into a no-molecule zone (transmission = 1)
    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new


def spectra_fit_2_molecules(w,k0, c1, c2):
    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2]
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
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")

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


def spectra_fit_3_molecules(w,k0, c1, c2, c3):
    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3]
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
    spec_new.apply_slit(slit_size, unit="cm-1", norm_by="area", inplace=True, shape="gaussian")

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


def spectra_fit_4_molecules(w,k0, c1, c2, c3, c4):
    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4]
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


def spectra_fit_5_molecules(w,k0, c1, c2, c3, c4, c5):
    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4, c5]
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


def spectra_fit_6_molecules(w,k0, c1, c2, c3, c4, c5, c6):
    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4, c5, c6]
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


def spectra_fit_7_molecules(w,k0, c1, c2, c3, c4, c5, c6, c7):
    """
    See explanation of "spectra_fit_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4, c5, c6, c7]
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


def spectra_fit_8_molecules(w,k0, c1, c2, c3, c4, c5, c6, c7, c8):
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