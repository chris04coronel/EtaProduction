# Imports
import sys
import os
import pdb

import numpy as np # type: ignore
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import EtaTools as et
import tltl_tools as tlt

import awkward as ak
import uproot
import vector
import ROOT as root

# Hard Code
# Line 26 to distingush the name of particle from root file line 14-15

# Example = [rootfile, particle]
particle = ['Eta4Mu.root', 'Eta']
#particle = ['EtaPrime4Mu.root', 'EtaPrime']
#particle = ['JPsi4Mu.root', 'JPsi']
#particle = ['Eta2Mu2E.root', 'Eta_2']

# Open the file
ifile = uproot.open("pythia/RootFiles/" + particle[0])
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

# Histogram 1 inv mass of individual muons and electrons pairs
tlt.two_lep_inv_mass_hist(lep1, lep2, lep3, lep4, particle[1], 'Inv_mass_lep_pairs')
