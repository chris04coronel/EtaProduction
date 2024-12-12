import numpy as np
import matplotlib.pyplot as plt
import uproot
import vector
import time
import math
from scipy.optimize import curve_fit

#def two_mu_two_e_check():

def two_lep_inv_mass_hist(vec1, vec2, vec3, vec4, parent, savefig):
    ''' Notes:'''
    mu_inv_mass_list = []
    e_inv_mass_list =[]
    for i in range(len(vec1)):
        if vec1[i].eta >=2 and vec1[i].eta <= 5 and vec2[i].eta >=2 and vec2[i].eta <= 5 and vec3[i].eta >=2 and vec3[i].eta <=5 and vec4[i].eta >=2 and vec4[i].eta <= 5:
            E_mu = vec1[i].e + vec2[i].e
            px_mu = vec1[i].px + vec2[i].px
            py_mu = vec1[i].py + vec2[i].py
            pz_mu = vec1[i].pz + vec2[i].pz
            P2_mu = px_mu**2 + py_mu**2 +pz_mu**2

            E_e = vec3[i].e +vec4[i].e
            px_e = vec3[i].px + vec4[i].px
            py_e = vec3[i].py + vec4[i].py
            pz_e = vec3[i].pz + vec4[i].pz
            P2_e = px_e**2 + py_e**2 +pz_e**2

            mu_inv_mass = math.sqrt(abs(E_mu**2 - P2_mu))
            mu_inv_mass_list.append(mu_inv_mass)

            e_inv_mass = math.sqrt(abs(E_e**2 - P2_e))
            e_inv_mass_list.append(e_inv_mass)
            

    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    axs.xaxis.set_tick_params(pad = 20) 
    axs.yaxis.set_tick_params(pad = 20)

    bin_heights1, bin_borders1, _ = plt.hist(mu_inv_mass_list, bins='auto', label='Muon Inv Mass', alpha=0.8, color = 'cornflowerblue', histtype = 'step')
    bin_centers1 = bin_borders1[:-1] + np.diff(bin_borders1) / 2

    plt.axvline(np.array(mu_inv_mass_list).mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(np.array(mu_inv_mass_list).mean()*1.0, max_ylim*0.95, 'Muon Mean: {:.2f}'.format(np.array(mu_inv_mass_list).mean()))


    bin_heights1, bin_borders1, _ = plt.hist(e_inv_mass_list, bins='auto', label='Electron Inv Mass', alpha=0.8, color = 'rebeccapurple', histtype = 'step')
    bin_centers1 = bin_borders1[:-1] + np.diff(bin_borders1) / 2

    plt.axvline(np.array(e_inv_mass_list).mean(), color='k', linestyle='dashed', linewidth=1)
    plt.text(np.array(e_inv_mass_list).mean()*1.0, max_ylim*0.80, '$Electron Mean: {:.2f}'.format(np.array(e_inv_mass_list).mean()))

    plt.title('Invariant Mass of lepton pairs')
    plt.xlabel('Mass (GeV)')
    plt.ylabel('Number of Candidates Per Mass')
    plt.legend()
    plt.savefig('pythia/Histograms/'+ parent + '/'+ savefig)