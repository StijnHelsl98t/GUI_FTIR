"""
File used to calculate the center, sigma and gamma of the ice-bands. This way, the created function can also be used in
other fitting procedures.
"""

from Code_Gui.Gui_General_Code import General_Functions_Library as GFL
import numpy as np
from lmfit.models import  ExponentialModel, SkewedGaussianModel, Linea
import os
import matplotlib.pyplot as plt
import os.path as op
from radis import Spectrum
import scipy

def fit_baseline(w, a, a_b, bool):
    def baseline_creator(x, A1, center, sigma, gamma, A2, decay):
        if A1 > 0:
            skewed = (A1/(sigma*np.sqrt(2*np.pi)))*np.exp(-1*(x-center)**2/2/sigma**2)*\
                     (1+scipy.special.erf(gamma*(x-center)/(sigma*np.sqrt(2))))
        else:
            skewed = np.zeros(len(x))

        exp_now = A2*np.exp(-1*x/decay)

        full_baseline = skewed + exp_now
        return full_baseline

    bas, b = GFL.baseline_correction_calculator(a, 0.0001, 10000000, 0.5)
    bas_b, b_b = GFL.baseline_correction_calculator(a_b, 0.0001, 1000000, 0.5)

    a_background_removed = a-bas_b
    a_zero = a_b - bas_b

    bas_zero, b_zero = GFL.baseline_correction_calculator(a_zero, 0.0001, 1000000, 0.5)
    bas_background_removed, b_background_removed = GFL.baseline_correction_calculator(
        a_background_removed, 0.0001, 1000000, 0.5)

    bas_line_temp = []
    bas_w = []

    for k in range(len(w)):
        if (np.abs(a_background_removed[k] - bas_background_removed[k])) <= 0.00005:
            bas_w.append(w[k])
            bas_line_temp.append(a_background_removed[k])

    bas_w = np.array(bas_w)
    bas_line_temp = np.array(bas_line_temp)

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

    baseline = baseline_creator(x=w, A1=g_amplitude, center=g_center, sigma=g_sigma, gamma=g_gamma, A2=e_amplitude,
                                 decay = e_decay)

    #plt.plot(w, a_background_removed)
    #plt.plot(w, baseline, alpha=0.7)
    #plt.plot(w, bas_zero, alpha=0.7)
    #plt.show()

    a_new = a-baseline
    bas_line_new, bas_new = GFL.baseline_correction_calculator(a_new, 0.0001, 10000000, 0.5)

    error = []
    for k in range(len(w)):
        if bas_line_new[k] <= bas_zero[k]:
            error.append(bas_line_new[k])

    a_new += np.abs(np.mean(error))

    #plt.plot(w, a_new, alpha=0.7)
    #plt.plot(w, a, alpha=0.7)
    #plt.show()

    return baseline

"""
directory = "C:\\Users\\P70085588\\Data\\Invenio-r\\2023\\07\\19\\MCT_time_effects\\"

files_in_directory_temp = [f for f in os.listdir(directory) if op.isfile(op.join(directory, f))]
array_data = GFL.read_opus_data_from_folder_into_array_for_gui(directory)

list_g_a = []
list_g_c = []
list_g_s = []
list_g_g = []
list_l_s = []
list_l_i = []
list_c = []
for i in range(0, len(files_in_directory_temp)):
    data_text_file_selected = directory + files_in_directory_temp[i]
    data = array_data[i]
    print(i, len(files_in_directory_temp)-1)
    if os.path.exists(data_text_file_selected) and os.stat(data_text_file_selected).st_size > 0:
        if i ==0:
            data_wavelength = np.array(data.get_range("AB")[0:-1])
            data_transmission = np.array(data["AB"][0:-1])

            s_background = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                             units={"transmittance": ""})
            w_background, t_background = s_background.get("transmittance")

            wlmin = min(data_wavelength)
            wlmax = max(data_wavelength)

            wlmin_j = len(w_background)
            wlmax_j=0
            for j in range(len(w_background)):
                if w_background[j] >= 1000:
                    wlmin = w_background[j]
                    wlmin_j = j
                if w_background[j] >= 4200:
                    wlmax = w_background[j]
                    wlmax_j = j

            a_background = GFL.tr_to_ab(t_background[wlmax_j:wlmin_j])

        else:
            data_wavelength = np.array(data.get_range("AB")[0:-1])
            data_transmission = np.array(data["AB"][0:-1])

            s_exp = Spectrum({"wavenumber": data_wavelength, "transmittance": data_transmission}, wunit='cm-1',
                             units={"transmittance": ""})
            w_exp, t_exp = s_exp.get("transmittance")

            a_exp = GFL.tr_to_ab(t_exp[wlmax_j:wlmin_j])

            bas, bas_b = fit_baseline(w_exp[wlmax_j:wlmin_j], a_exp, a_background, True)

            a_exp_baseline_corrected = a_exp - bas

            plt.plot(w_exp[wlmax_j:wlmin_j], a_exp, alpha=0.7, label="experiment")
            plt.plot(w_exp[wlmax_j:wlmin_j], a_exp_baseline_corrected, alpha=0.7, label="baseline_corrected")

            plt.legend()
            plt.show()
            #list_g_a.append(g_a)
            #list_g_c.append(g_c)
            #list_g_s.append(g_s)
            #list_g_g.append(g_g)
            #list_l_i.append(l_i)
            #list_l_s.append(l_s)

            #plt.show()

print("center:", list_g_c)
print("sigma:", list_g_s)

"""
