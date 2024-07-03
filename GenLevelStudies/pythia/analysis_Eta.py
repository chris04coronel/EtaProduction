# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot
import vector

# Personal Packages
sys.path.append(".")
import AnalysisTools as at

# Parse Inputs
parser = at.Parser(sys.argv[1:])        # Create the parser
args = parser.args                      # Parser parses the args lol
file_name =  args.input_files[0].name  
cross_section = args.cross_section # Cross section in fb

# Open the file
ifile = uproot.open(file_name)
tree = ifile['Tree'].arrays()

# Create Vectors, They open up branches
leps_vec=[lep1_vec, lep2_vec, lep3_vec, lep4_vec]
TargetLeptons=['TargetLepton1','TargetLepton2','TargetLepton3','TargetLepton4']
for i in range(len(leps_vec)):
    leps_vec[i] = vector.zip({
    'px': tree[TargetLeptons[i]].px,
    'py': tree[TargetLeptons[i]].py,
    'pz': tree[TargetLeptons[i]].pz,
    'e': tree[TargetLeptons[i]].energy,
    'pid': tree[TargetLeptons[i]].id
})

# lepton_vec = vector.zip({
#     'px': tree['TargetLepton'].px,
#     'py': tree['TargetLepton'].py,
#     'pz': tree['TargetLepton'].pz,
#     'e': tree['TargetLepton'].energy,
#     'pid': tree['TargetLepton'].id
# })

eta_vec = vector.zip({
    'px': tree['TargetEta'].px,
    'py': tree['TargetEta'].py,
    'pz': tree['TargetEta'].pz,
    'e': tree['TargetEta'].energy,
    'pid': tree['TargetEta'].id
})

# Masks
# Demanding that vecs meet some kind of criteria
# How many events have all four muons meeting LHCb acceptance 
#might need to make a for loop that loops through all leptons
both_product_loose_acc_mask = ((lepton_vec.eta>2)
                                & (lepton_vec.eta<5)
                                & (jet_vec.eta>2)
                                & (jet_vec.eta<5)) 
both_product_tight_acc_mask = ((lepton_vec.eta>2.2)
                                & (lepton_vec.eta<4.4)
                                & (jet_vec.eta>2.2)
                                & (jet_vec.eta<4.4)) 
high_pT_product_mask = ((lepton_vec.pt>20) & (jet_vec.pt>20))
low_pT_product_mask = ((lepton_vec.pt>5) & (jet_vec.pt>5))
one_product_gauss_mask = ((lepton_vec.eta>1.596) & (lepton_vec.pt>17))
product_mask = both_product_tight_acc_mask&high_pT_product_mask&low_pT_product_mask

# Print
print(f"Total number of events: {len(lepton_vec)}")
print(f"Gauss Cuts: {sum((one_product_gauss_mask))}")
print(f"DaVinci Cut: {sum((low_pT_product_mask))}")
print(f"Total Cuts: {sum((one_product_gauss_mask&low_pT_product_mask))}")
print(f"Loose Acceptance Cuts: {sum(both_product_loose_acc_mask)}")
print(f"Tight Acceptance Cuts: {sum(product_mask)}")
