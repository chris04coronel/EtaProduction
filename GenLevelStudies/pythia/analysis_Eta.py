# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt
import EtaTools as et

import awkward as ak
import uproot
import vector
import ROOT as root

# Open the file
ifile = uproot.open("pythia/EtaProduction.root")
tree = ifile['Tree_eta'].arrays()

# four Vectors
eta_vec = vector.zip({
    'px': tree['Eta'].px,
    'py': tree['Eta'].py,
    'pz': tree['Eta'].pz,
    'e': tree['Eta'].energy
})
eta_mu_1_vec = vector.zip({
    'px': tree['Muon1'].px,
    'py': tree['Muon1'].py,
    'pz': tree['Muon1'].pz,
    'e': tree['Muon1'].energy
})
eta_mu_2_vec = vector.zip({
    'px': tree['Muon2'].px,
    'py': tree['Muon2'].py,
    'pz': tree['Muon2'].pz,
    'e': tree['Muon2'].energy
})
eta_mu_3_vec = vector.zip({
    'px': tree['Muon3'].px,
    'py': tree['Muon3'].py,
    'pz': tree['Muon3'].pz,
    'e': tree['Muon3'].energy
})
eta_mu_4_vec = vector.zip({
    'px': tree['Muon4'].px,
    'py': tree['Muon4'].py,
    'pz': tree['Muon4'].pz,
    'e': tree['Muon4'].energy
})

# Histogram 1 pT
muons_pt = et.LApt(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(muons_pt ,edgecolor='black', color='teal', bins=35, range=[0,1])
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('pT Among Muons')
plt.xlabel('Energy (GeV)')
plt.ylabel('Number of Particles per pT')
plt.savefig('4Muons_pT')

# Histogram 2 all 4 muons in acceptance
muons_eta = et.LAeta(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)

fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(muons_eta ,edgecolor='black', color='teal', bins=35)
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('eta Among Individual Muons')
plt.xlabel('eta')
plt.ylabel('Number of Particles per eta')
plt.savefig('4Muons_eta')

# Histogram 3 eta spread among sibling muons
#2d array each row is 4 sibling muons and their listed etas
sibmuons = et.TwoDEta(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
counter = 0
delta_eta_list = []
for i in range(len(sibmuons[:,1])):
    mineta = 0
    maxeta = 0
    sibs = sibmuons[i,:]
    mineta = np.min(sibs)
    maxeta = np.max(sibs)
    delta_eta = maxeta-mineta
    delta_eta_list.append(delta_eta)
    if sibs[0] <= 5 and sibs[0] >= 2 and sibs[1] <= 5 and sibs[1] >= 2 and sibs[2] <= 5 and sibs[2] >= 2 and sibs[3] <= 5 and sibs[3] >= 2:
        counter += 1
        if counter == 1:
            sibs_eta_acc = sibmuons[i,:]
        else:
            sibs_eta_acc = np.vstack((sibs_eta_acc,sibs))

fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(delta_eta_list ,edgecolor='black', color='teal', bins=35)
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('Greatest Delta Eta Among Sibling Muons')
plt.xlabel('Delta Eta')
plt.ylabel('Number of Particles per Delta eta')
plt.savefig('4Muons_delta_eta')

# Histogram 4 
# Delta eta when all siblings are in Acceptance
delta_eta_acc_list = []
for i in range(len(sibs_eta_acc[:,1])):
    mineta = 0
    maxeta = 0
    mineta = np.min(sibs_eta_acc[i,:])
    maxeta = np.max(sibs_eta_acc[i,:])
    delta_eta_acc = maxeta-mineta
    delta_eta_acc_list.append(delta_eta_acc)

fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(delta_eta_acc_list ,edgecolor='black', color='teal', bins=35)
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('Greatest Delta Eta Among Sibling Muons All in Acceptance')
plt.xlabel('Delta Eta')
plt.ylabel('Number of Particles per Delta eta')
plt.savefig('4Muons_delta_eta_acc')



print('The Highest pT among all muons is', muons_pt.max(),'GeV')
print('The Lowest pT among all muons is', muons_pt.min(),'GeV')
print('The avaerage pseudorapidty spread among siblings is', np.array(delta_eta_list).mean())
print('The greatest pseudorapidty spread among siblings is', np.array(delta_eta_list).max())
print('The smallest pseudorapidty spread among siblings is', np.array(delta_eta_list).min())

print('\nThe avaerage pseudorapidty spread among siblings all in acceptance is', np.array(delta_eta_acc_list).mean())
print('The greatest pseudorapidty spread among siblings all in acceptance is', np.array(delta_eta_acc_list).max())
print('The smallest pseudorapidty spread among siblings all in acceptance is', np.array(delta_eta_acc_list).min())
