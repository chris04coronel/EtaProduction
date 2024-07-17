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
treeprime = ifile['Tree_eta_prime'].arrays()

# Eta and it's Muons: Four Vectors
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

# Eta_Prime and it's Muons: Four Vectors
EtaPVec = vector.zip({
    'px': treeprime['EtaPrime'].px,
    'py': treeprime['EtaPrime'].py,
    'pz': treeprime['EtaPrime'].pz,
    'e': treeprime['EtaPrime'].energy
})
PMuVec1 = vector.zip({
    'px': treeprime['MuonPrime1'].px,
    'py': treeprime['MuonPrime1'].py,
    'pz': treeprime['MuonPrime1'].pz,
    'e': treeprime['MuonPrime1'].energy
})
PMuVec2 = vector.zip({
    'px': treeprime['MuonPrime2'].px,
    'py': treeprime['MuonPrime2'].py,
    'pz': treeprime['MuonPrime2'].pz,
    'e': treeprime['MuonPrime2'].energy
})
PMuVec3 = vector.zip({
    'px': treeprime['MuonPrime3'].px,
    'py': treeprime['MuonPrime3'].py,
    'pz': treeprime['MuonPrime3'].pz,
    'e': treeprime['MuonPrime3'].energy
})
PMuVec4 = vector.zip({
    'px': treeprime['MuonPrime4'].px,
    'py': treeprime['MuonPrime4'].py,
    'pz': treeprime['MuonPrime4'].pz,
    'e': treeprime['MuonPrime4'].energy
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
plt.savefig('pythia/Histograms/4Muons_pT')

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
plt.savefig('pythia/Histograms/4Muons_eta')

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
plt.savefig('pythia/Histograms/4Muons_delta_eta')

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
plt.savefig('pythia/Histograms/4Muons_delta_eta_acc')


print('The Highest pT among all muons is', muons_pt.max(),'GeV')
print('The Lowest pT among all muons is', muons_pt.min(),'GeV')
print('The avaerage pseudorapidty spread among siblings is', np.array(delta_eta_list).mean())
print('The greatest pseudorapidty spread among siblings is', np.array(delta_eta_list).max())
print('The smallest pseudorapidty spread among siblings is', np.array(delta_eta_list).min())

print('\nThe avaerage pseudorapidty spread among siblings all in acceptance is', np.array(delta_eta_acc_list).mean())
print('The greatest pseudorapidty spread among siblings all in acceptance is', np.array(delta_eta_acc_list).max())
print('The smallest pseudorapidty spread among siblings all in acceptance is', np.array(delta_eta_acc_list).min())

################
################
################
################ 

# Prime Histogram 1 pT
prime_muons_pt = et.LApt(PMuVec1, PMuVec2, PMuVec3, PMuVec4)
fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(prime_muons_pt ,edgecolor='black', color='rebeccapurple', bins=35, range=[0,1])
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('pT Among Muons from Eta Prime')
plt.xlabel('Energy (GeV)')
plt.ylabel('Number of Particles per pT')
plt.savefig('pythia/Histograms/prime_4Muons_pT')

# Prime Histogram 2 all 4 muons in acceptance
prime_muons_eta = et.LAeta(PMuVec1, PMuVec2, PMuVec3, PMuVec4)

fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(prime_muons_eta ,edgecolor='black', color='rebeccapurple', bins=35)
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('Eta Among Individual Muons from Eta Prime')
plt.xlabel('eta')
plt.ylabel('Number of Particles per eta')
plt.savefig('pythia/Histograms/prime_4Muons_eta')

# Prime Histogram 3 eta spread among sibling muons
# 2d array each row is 4 sibling muons and their listed etas
prime_sibmuons = et.TwoDEta(PMuVec1, PMuVec2, PMuVec3, PMuVec4)
counter = 0
prime_delta_eta_list = []
for i in range(len(prime_sibmuons[:,1])):
    mineta = 0
    maxeta = 0
    sibs = prime_sibmuons[i,:]
    mineta = np.min(sibs)
    maxeta = np.max(sibs)
    delta_eta = maxeta-mineta
    prime_delta_eta_list.append(delta_eta)
    if sibs[0] <= 5 and sibs[0] >= 2 and sibs[1] <= 5 and sibs[1] >= 2 and sibs[2] <= 5 and sibs[2] >= 2 and sibs[3] <= 5 and sibs[3] >= 2:
        counter += 1
        if counter == 1:
            prime_sibs_eta_acc = sibmuons[i,:]
        else:
            prime_sibs_eta_acc = np.vstack((prime_sibs_eta_acc,sibs))

fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(delta_eta_list ,edgecolor='black', color='rebeccapurple', bins=35)
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('Greatest Delta Eta Among Sibling Muons From Eta Prime')
plt.xlabel('Delta Eta')
plt.ylabel('Number of Particles per Delta eta')
plt.savefig('pythia/Histograms/prime_4Muons_delta_eta')

# Prime Histogram 4 
# Delta eta when all siblings are in Acceptance
prime_delta_eta_acc_list = []
for i in range(len(prime_sibs_eta_acc[:,1])):
    pmineta = 0
    pmaxeta = 0
    pmineta = np.min(prime_sibs_eta_acc[i,:])
    pmaxeta = np.max(prime_sibs_eta_acc[i,:])
    prime_delta_eta_acc = pmaxeta-pmineta
    prime_delta_eta_acc_list.append(prime_delta_eta_acc)

fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(b = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
counts, edges, bars = plt.hist(prime_delta_eta_acc_list ,edgecolor='black', color='rebeccapurple', bins=35)
plt.bar_label(bars)
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('Greatest Delta Eta Among Sibling Muons All in Acceptance from Eta Prime')
plt.xlabel('Delta Eta')
plt.ylabel('Number of Particles per Delta eta')
plt.savefig('pythia/Histograms/prime_4Muons_delta_eta_acc')

print('The Highest pT among all muons from eta_prime is', prime_muons_pt.max(),'GeV')
print('The Lowest pT among all muons from eta_pime is', prime_muons_pt.min(),'GeV')
print('The avaerage pseudorapidty spread among siblings from eta prime is', np.array(prime_delta_eta_list).mean())
print('The greatest pseudorapidty spread among siblings from eta prime is', np.array(prime_delta_eta_list).max())
print('The smallest pseudorapidty spread among siblings from eta prime is', np.array(prime_delta_eta_list).min())

print('\nThe avaerage pseudorapidty spread among siblings all in acceptance from eta prime is', np.array(prime_delta_eta_acc_list).mean())
print('The greatest pseudorapidty spread among siblings all in acceptance from eta prime is', np.array(prime_delta_eta_acc_list).max())
print('The smallest pseudorapidty spread among siblings all in acceptance from eta prime is', np.array(prime_delta_eta_acc_list).min())
#pdb.set_trace