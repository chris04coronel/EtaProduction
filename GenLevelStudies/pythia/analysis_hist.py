# Imports
import sys
import os
import pdb

import numpy as np # type: ignore
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import EtaTools as et

import awkward as ak
import uproot
import vector
import ROOT as root

# Hard Code
# Example = [rootfile, particle, particle]
#particle = ['Eta4Mu.root', 'Eta', 'Eta']
#particle = ['EtaPrime4Mu.root', 'EtaPrime', 'EtaPrime']
#particle = ['JPsi4Mu.root', 'JPsi', 'JPsi]
#particle = ['Eta2Mu2E.root', 'Eta_2', 'Eta_2]
# phase space
particle = ['Eta4Mu_R3.root', 'Eta', 'Eta_R3', 'Run3']


# Open the file
#ifile = uproot.open("pythia/RootFiles/" + particle[0])
ifile = uproot.open("RootFiles/" + particle[0])
tree = ifile['Tree_Particle'].arrays()

# Eta and it's Muons: Four Vectors
particle_vec = vector.zip({
    'px': tree[particle[1]].px,
    'py': tree[particle[1]].py,
    'pz': tree[particle[1]].pz,
    'e': tree[particle[1]].energy
})
lep1 = vector.zip({
    'px': tree['Muon1'].px,
    'py': tree['Muon1'].py,
    'pz': tree['Muon1'].pz,
    'e': tree['Muon1'].energy
})
lep2 = vector.zip({
    'px': tree['Muon2'].px,
    'py': tree['Muon2'].py,
    'pz': tree['Muon2'].pz,
    'e': tree['Muon2'].energy
})
lep3 = vector.zip({
    'px': tree['Muon3'].px,
    'py': tree['Muon3'].py,
    'pz': tree['Muon3'].pz,
    'e': tree['Muon3'].energy
})
lep4 = vector.zip({
    'px': tree['Muon4'].px,
    'py': tree['Muon4'].py,
    'pz': tree['Muon4'].pz,
    'e': tree['Muon4'].energy
})

####################
####################
####################
####################
mb_arr = et.min_bias_arr(particle_vec, lep1, lep2, lep3, lep4)

# # Histogram 1/4 Invariant mass reconstruction of 4 muons
# particle_inv_mass = et.inv_mass_recon_list(lep1, lep2, lep3, lep4)
# arr1=np.array(particle_inv_mass)
# et.hist_mass_curve(arr1, particle[2], particle[1]+'_inv_mass', '$m_{4\mu}$ Invariant Mass ', 'teal', 200)

# # Historgram 2/4 Creates 4 histograms on one plot, for each lepton pT.
# plt.clf()
# muon_sorted_pt = et.eta_acc_sort_pt(particle_vec, lep1, lep2, lep3, lep4)
# et.hist_pT_4mu(particle[2], '4muons_mean_pT', particle[2]+' Muon Siblings in Acc pT Spread', muon_sorted_pt)

# # Histogram 3/4 pT
# muons_pt = et.LApt(lep1, lep2, lep3, lep4)
# et.hist(muons_pt, particle[2], '4Muons_pT', 'pT Among All Muons', 'Energy (GeV)', 'Number of Particles per pT', 'teal', 300)

# # Histogram 4/4 all 4 muons 
# muons_eta = et.LAeta(lep1, lep2, lep3, lep4)
# et.hist(muons_eta, particle[2], '4MuonsEta', 'eta Among Individual Muons', 'eta', 'Number pf particle per eta','teal',  200 )

# Histogram pT distributions
mu1pass = et.mb_pt_check1(mb_arr, .45)
mu12pass = et.mb_pt_check12(mb_arr, .45)
mu123pass = et.mb_pt_check123(mb_arr, .45)

et.hist(mu1pass[:,2], particle[2], particle[1]+'mu2_pt_distrib'+particle[3], particle[1]+'Second Fastest Muon pT distribution'+particle[3], 'pT (Gev)', 'Candidates', 'black', 200)
et.hist(mu1pass[:,1], particle[2],  particle[1]+'mu3_pt_distrib'+particle[3], particle[1]+'Third Fastest Muon pT distribution'+particle[3], 'pT (Gev)', 'Candidates', 'black', 200 )
et.hist(mu1pass[:,0], particle[2],  particle[1]+'mu4_pt_distrib'+particle[3], particle[1]+'Fourth Fastest Muon pT distribution'+particle[3], 'pT (Gev)', 'Candidates', 'black', 200 )


# Histogram XX (All Muons in Acceptance)
# all_muons_acc_pt = et.eta_acc_LApt(lep1, lep2, lep3, lep4)
# et.hist(all_muons_acc_pt, particle[1], '4muons_acc_pT', 'Muons (all in acceptance) pT', 'pT in GeV','Number of Particles per pT', 'teal', 65 )

#pdb.set_trace()

# # Histogram XX eta spread among sibling muons
# # 2d array each row is 4 sibling muons and their listed etas
# sibmuons = et.TwoDEta(lep1, lep2, lep3, lep4)
# counter = 0
# delta_eta_list = []
# for i in range(len(sibmuons[:,1])):
#     mineta = 0
#     maxeta = 0
#     sibs = sibmuons[i,:]
#     mineta = np.min(sibs)
#     maxeta = np.max(sibs)
#     delta_eta = maxeta-mineta
#     delta_eta_list.append(delta_eta)
#     if sibs[0] <= 5 and sibs[0] >= 2 and sibs[1] <= 5 and sibs[1] >= 2 and sibs[2] <= 5 and sibs[2] >= 2 and sibs[3] <= 5 and sibs[3] >= 2:
#         counter += 1
#         if counter == 1:
#             sibs_eta_acc = sibmuons[i,:]
#         else:
#             sibs_eta_acc = np.vstack((sibs_eta_acc,sibs))
# et.hist(delta_eta_list, '4MuonsEtaSpread', 'Eta Spread Among all eta siblings', 'Delta eta', 'Number of Particles per Delta eta','teal', 200)

# # Histogram XX 
# # Delta eta when all siblings are in Acceptance
# delta_eta_acc_list = []
# for i in range(len(sibs_eta_acc[:,1])):
#     accmineta = 0
#     accmaxeta = 0
#     accsibs = sibs_eta_acc[i,:]
#     accmineta = np.min(accsibs)
#     accmaxeta = np.max(accsibs)
#     acc_delta = accmaxeta - accmineta
#     delta_eta_acc_list.append(acc_delta)
# et.hist(delta_eta_acc_list, '4MuonsAccEtaSpread', 'Eta Spread Among all eta siblings in Acceptance', 'Delta eta', 'Number of Particles per Delta eta','teal', 200)