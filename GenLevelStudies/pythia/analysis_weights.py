# Imports
import sys
import os
import pdb

import numpy as np # type: ignore
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tabulate import tabulate
import math
import EtaTools as et

import awkward as ak
import uproot
import vector
import ROOT as root
import pythia8
# 2>&1 | tee Logs/.log
# Hard Code
# Example = [rootfile, particle, card file]
#particle = ['JPsi4Mu.root', 'JPsi', 'jpsi4mu.card', 'Run 2']
#particle = ['Eta2Mu2E.root', 'Eta_2']
# Run 3
particle = ['Eta4Mu_R3.root', 'Eta', 'eta4mu_r3.card', 'Run 3', 'Eta_R3', 5*10**-9]
#particle = ['EtaPrime4Mu_R3.root', 'EtaPrime', 'etaprime4mu_r3.card', 'Run 3', 'EtaP_R3', 1.69*10**-8 ]

# Open the file
ifile = uproot.open("RootFiles/" + particle[0])
tree = ifile['Tree_Particle'].arrays()


# Personal Packages
sys.path.append(".")
import AnalysisTools_etameson as at
eta_path = at.find_ETA_path()
# Set up makefile configuration
# pythia tells us to include this 
cfg = open(eta_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/Applications/pythia8310/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)
# Read in Card File
# How to read in a card file. Look up in a tutorial.
pythia = pythia8.Pythia()
pythia.readFile(eta_path + "/GenLevelStudies/pythia/CardFiles/" + particle[2])
# Initialize Pythia
# Starting up pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents") 

# Eta and it's Muons: Four Vectors
parent_vec = vector.zip({
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
#Calculating the Weights
N = 9.09*(10**14)
br = particle[5]
weights1 = N*br/nEvent

# Creates Arrays to plotted in a Histogram
acc_arr = et.LHCb_acc(parent_vec, lep1, lep2, lep3, lep4)
sorted_pt = et.sort_pt(acc_arr)
#pdb.set_trace()
# Checks fraction that make it past the pT trigger and returns new Array
mu1234pass = et.mb_pt_check1234(sorted_pt, .45, .25)

#Print Values 
nDec = len(parent_vec)                                       
nPDAcc = len(acc_arr[:,0])    # Takes both Parent and Daughter in Eta Acceptance
print(particle[1],particle[3],'Statistics')
print('\nThe total number of events is', nEvent)
print('The total number of',particle[1],'decay to 4 muons is', nDec)
print('The total number of decays where both the daughters and the parent particle are in LHCb acceptance is', nPDAcc)
print('The average number of decays with parent and daughters in acceptance per event is', nPDAcc/nEvent)
print('\nThe Weighted value is equal to N', "{:.2e}".format(N),'* br', br,'/', nEvent,'=', weights1)

#Weighted pT spread of muons before pT trigger cuts
et.hist_weight(sorted_pt[:,3], particle[4], particle[1]+'_weighted_mu1_pt_distrib_'+particle[3]+'_no_cuts', particle[1]+' Largest Muon pT Distribution '+ particle[3], 'pT (Gev)', 'Candidates', 'darkblue', 200, 'LHCb Run 3 Before Cuts',weights1)
et.hist_weight(sorted_pt[:,2], particle[4], particle[1]+'_weighted_mu2_pt_distrib_'+particle[3]+'_no_cuts',particle[1]+' Second largest Muon pT Distribution ' + particle[3], 'pT (Gev)', 'Candidates', 'darkblue', 200, 'LHCb Run 3 Before Cuts',weights1) 
et.hist_weight(sorted_pt[:,1], particle[4], particle[1]+'_weighted_mu3_pt_distrib_'+particle[3]+'_no_cuts',particle[1]+' Third Largest Muon pT Distribution '+ particle[3], 'pT (Gev)', 'Candidates', 'darkblue', 200, 'LHCb Run 3 Before Cuts',weights1)
et.hist_weight(sorted_pt[:,0], particle[4], particle[1]+'_weighted_mu4_pt_distrib_'+particle[3]+'_no_cuts', particle[1]+' Fourth Largest Muon pT Distribution '+ particle[3], 'pT (Gev)', 'Candidates', 'darkblue', 200, 'LHCb Run 3 Before Cuts',weights1)

#Weighted pT spread of muons before pT trigger cuts
et.hist_weight(mu1234pass[:,3], particle[4], particle[1]+'_weighted_mu1_pt_distrib_'+particle[3]+'_with_cuts', particle[1]+' Largest Muon pT Distribution '+ particle[3], 'pT (Gev)', 'Candidates', 'black', 200, 'LHCb Run 3 After Cuts',weights1)
et.hist_weight(mu1234pass[:,2], particle[4], particle[1]+'_weighted_mu2_pt_distrib_'+particle[3]+'_with_cuts',particle[1]+' Second largest Muon pT Distribution ' + particle[3], 'pT (Gev)', 'Candidates', 'black', 200, 'LHCb Run 3 After Cuts',weights1) 
et.hist_weight(mu1234pass[:,1], particle[4], particle[1]+'_weighted_mu3_pt_distrib_'+particle[3]+'_with_cuts',particle[1]+' Third Largest Muon pT Distribution '+ particle[3], 'pT (Gev)', 'Candidates', 'black', 200, 'LHCb Run 3 After Cuts',weights1)
et.hist_weight(mu1234pass[:,0], particle[4], particle[1]+'_weighted_mu4_pt_distrib_'+particle[3]+'_with_cuts', particle[1]+' Fourth Largest Muon pT Distribution '+ particle[3], 'pT (Gev)', 'Candidates', 'black', 200, 'LHCb Run 3 After Cuts',weights1)
#pdb.set_trace()


# Events are those that passed pT trigger cuts. To reconstuct the total events we divide events by 
# fraction that made the pT trigger cuts in stand alone pythia. 
particle_inv_mass = et.inv_mass_recon_list(parent_vec, lep1, lep2, lep3, lep4)
et.hist_weight_mass_curve(particle_inv_mass, particle[1], particle[1]+'_weighted_inv_mass', 'Weighted Reconstructed Invariant Mass', 'teal', 100, None, weights1)

# deletes these lines below
counts, bins, _ = plt.hist(particle_inv_mass, bins=200, histtype = 'step')
bin_widths = np.diff(bins)
integral = np.sum(counts * bin_widths)

pdb.set_trace()