from radis import Spectrum
import numpy as np
import General_Functions_Library as GFL

w_dict = {}
t_dict = {}

def storage_for_dict(w,t, p):
    """
    This function creates a storage place, so that the wavenumber and transmission gained for the spectra can be
    imported into this module.
    :param w: The imported dictionary containing the wavenumber-arrays
    :param t: The imported dictionary containing the transmission-arrays
    :return: -
    """
    global w_dict, t_dict, pathlength
    w_dict = w
    t_dict = t
    pathlength = p

def spectra_linear_wavenumber_dependent_offset_1_molecule(w, c1, offset_min, offset_max, slit_size):
    """
    This fuction is used in order to fit a calculated spectra to the gained experimental data. This formula can only be
    used when one wants to fit a single molucele. Because of the way scipy.curve_fit() works, this function is repeated
    below, but with more c's, so more molecules can be fitted at the same time.

    :param w: wavenumber range needed to match simulated with experimental data
    :param c1: mole fraction of the molecule one wants to fit
    :param offset_min: Minimimal offset-error compared with simulated data
    :param offset_max: Maximal offset-error compared with simulated data
    :param slit_size: Size of the needed slit to match simulated with experimental data
    :return: return fitted transmittance
    """
    global w_dict, t_dict
    c_list=[c1]
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
    offset_max_temp = offset_max
    if offset_max < offset_min:
        offset_max = offset_min
        offset_min = offset_max_temp
    offset_lin_list = np.linspace(offset_min, offset_max, len(w_temp))

    w_new_temp = np.zeros(len(w_temp))
    for i in range(len(w_temp)):
        w_new_temp[i] = w_temp[i] + offset_lin_list[i]

    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")

    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new

def spectra_linear_wavenumber_dependent_offset_2_molecules(w, c1, c2, offset_min, offset_max, slit_size):
    """
    See explanation of "spectra_linear_wavenumber_dependent_offset_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2]
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
    offset_max_temp = offset_max
    if offset_max < offset_min:
        offset_max = offset_min
        offset_min = offset_max_temp
    offset_lin_list = np.linspace(offset_min, offset_max, len(w_temp))

    w_new_temp = np.zeros(len(w_temp))
    for i in range(len(w_temp)):
        w_new_temp[i] = w_temp[i] + offset_lin_list[i]

    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")

    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new

def spectra_linear_wavenumber_dependent_offset_3_molecules(w, c1, c2, c3, offset_min, offset_max, slit_size):
    """
    See explanation of "spectra_linear_wavenumber_dependent_offset_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3]
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

    offset_max_temp = offset_max
    if offset_max < offset_min:
        offset_max = offset_min
        offset_min = offset_max_temp
    offset_lin_list = np.linspace(offset_min, offset_max, len(w_temp))

    w_new_temp = np.zeros(len(w_temp))
    for i in range(len(w_temp)):
        w_new_temp[i] = w_temp[i] + offset_lin_list[i]

    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")

    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new

def spectra_linear_wavenumber_dependent_offset_4_molecules(w, c1, c2, c3,c4, offset_min, offset_max, slit_size):
    """
    See explanation of "spectra_linear_wavenumber_dependent_offset_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4]
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
    offset_max_temp = offset_max
    if offset_max < offset_min:
        offset_max = offset_min
        offset_min = offset_max_temp
    offset_lin_list = np.linspace(offset_min, offset_max, len(w_temp))

    w_new_temp = np.zeros(len(w_temp))
    for i in range(len(w_temp)):
        w_new_temp[i] = w_temp[i] + offset_lin_list[i]

    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")

    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new

def spectra_linear_wavenumber_dependent_offset_5_molecules(w, c1, c2, c3,c4, c5, offset_min, offset_max, slit_size):
    """
    See explanation of "spectra_linear_wavenumber_dependent_offset_1_molecule".
    """
    global w_dict, t_dict
    c_list=[c1, c2, c3, c4, c5]
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
    offset_max_temp = offset_max
    if offset_max < offset_min:
        offset_max = offset_min
        offset_min = offset_max_temp
    offset_lin_list = np.linspace(offset_min, offset_max, len(w_temp))

    w_new_temp = np.zeros(len(w_temp))
    for i in range(len(w_temp)):
        w_new_temp[i] = w_temp[i] + offset_lin_list[i]

    spec_new = Spectrum({"wavenumber": w_new_temp, "transmittance": t_temp}, wunit='cm-1',
                        units={"transmittance": ""}, conditions={"path_length": pathlength})

    spec_new.resample(w, inplace=True)
    w_new, t_new = spec_new.get("transmittance")

    for i in range(len(t_new)):
        if np.isnan(t_new[i]):
            t_new[i] = 1
    return t_new
