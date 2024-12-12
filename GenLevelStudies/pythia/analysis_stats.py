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

# Hard Code
# Example = [rootfile, particle]
#particle = ['Eta4Mu.root', 'Eta']
particle = ['EtaPrime4Mu.root', 'EtaPrime']
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
muons_sort_pt = et.eta_acc_sort_pt(lep1, lep2, lep3, lep4)
#et.pT_check1(muons_sort_pt[:,0],.3)

# Table 1 muon Statistics
print('\nMuon pT information for',particle[1], 'parent particle')
table1=[['-'            ,'1st Muon', '2nd Muon', '3rd Muon', '4th Muon'], 
        ['Highest (GeV)', '{:.5f}'.format(muons_sort_pt[:,3].max()), '{:.5f}'.format(muons_sort_pt[:,2].max()), '{:.5f}'.format(muons_sort_pt[:,1].max()), '{:.5f}'.format(muons_sort_pt[:,0].max())],
        ['Lowest  (GeV)', '{:.5f}'.format(muons_sort_pt[:,3].min()), '{:.5f}'.format(muons_sort_pt[:,2].min()), '{:.5f}'.format(muons_sort_pt[:,1].min()), '{:.5f}'.format(muons_sort_pt[:,0].min())],
        ['Average (GeV)', '{:.5f}'.format(muons_sort_pt[:,3].mean()), '{:.5f}'.format(muons_sort_pt[:,2].mean()), '{:.5f}'.format(muons_sort_pt[:,1].mean()), '{:.5f}'.format(muons_sort_pt[:,0].mean())]]
print(tabulate(table1,tablefmt="fancy_grid", showindex="always"))



# Table 2 Muons meeting Prompt 1 Statistics
print('\nPrompt 1 \n Statistics for',particle[1], 'particle meeting prompt 1 J/Psi-->4mu conditions.')
table2 = [['-'                        , '1st Muon', '2nd Muon', '3rd Muon', '4th Muon', 'Percent'],
        ['All Muons pT >0.5GeV'       , 'x'       , 'x'       , 'x'       , 'x'       , et.pT_check1234(muons_sort_pt[:,3],muons_sort_pt[:,2],muons_sort_pt[:,1],muons_sort_pt[:,0], 0.5, 0.5, 0.5, 0.5)]]
print(tabulate(table2,tablefmt="fancy_grid", showindex="always"))



# Table 3 Muons meeting prompt 2 statistics
print('\nPrompt 2 \n Statistics for',particle[1], 'particle meeting prompt 2 J/Psi-->4mu conditions.')
table3 = [['-'                        , '1st Muon', '2nd Muon', '3rd Muon', '4th Muon', 'Percent'],
        ['pT > 0.50GeV and pT >0.3GeV', 'x'       , 'x'       , '-'      , '-'        , et.pT_check12(muons_sort_pt[:,3],muons_sort_pt[:,2],0.5,0.3)],
        ['pT > 0.50GeV and pT >0.3GeV', 'x'       , '-'       , 'x'      , '-'        , et.pT_check12(muons_sort_pt[:,3],muons_sort_pt[:,1],0.5,0.3)],
        ['pT > 0.50GeV and pT >0.3GeV', 'x'       , '-'       , '-'      , 'x'        , et.pT_check12(muons_sort_pt[:,3],muons_sort_pt[:,0],0.5,0.3)],
        ['pT > 0.50GeV and pT >0.3GeV', '-'       , 'x'       , 'x'      , '-'        , et.pT_check12(muons_sort_pt[:,2],muons_sort_pt[:,1],0.5,0.3)],
        ['pT > 0.50GeV and pT >0.3GeV', '-'       , 'x'       , '-'      , 'x'        , et.pT_check12(muons_sort_pt[:,2],muons_sort_pt[:,0],0.5,0.3)],
        ['pT > 0.50GeV and pT >0.3GeV', '-'       , '-'       , 'x'      , 'x'        , et.pT_check12(muons_sort_pt[:,1],muons_sort_pt[:,0],0.5,0.3)]]
print(tabulate(table3,tablefmt="fancy_grid", showindex="always"))
#pdb.set_trace()
#(et.pT_check12(muons_sort_pt[:,3],muons_sort_pt[:,2],0.5,0.3))