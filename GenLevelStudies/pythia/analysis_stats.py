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
# 2>&1 | tee /Logs/.log
# Hard Code
# Example = [rootfile, particle, card file]
#particle = ['Eta4Mu.root', 'Eta', 'eta4mu.card', 'Run 2']
#particle = ['EtaPrime4Mu.root', 'EtaPrime', 'etaprime4mu.card', 'Run 2']
particle = ['JPsi4Mu.root', 'JPsi', jpsi4mu.card, 'Run 2']
#particle = ['Eta2Mu2E.root', 'Eta_2']
#particle = ['Eta4Mu_phsp.root', 'Eta']

# Run 3
#particle = ['Eta4Mu_R3.root', 'Eta', 'eta4mu_r3.card', 'Run 3']
#particle = ['EtaPrime4Mu_R3.root', 'EtaPrime', 'etaprime4mu_r3.card', 'Run 3']

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

nDec = len(parent_vec)
nDecAcc = len(et.TwoDEtaAcc(lep1, lep2, lep3, lep4))
nPDAcc = et.PDEtaAcc(parent_vec, lep1, lep2, lep3, lep4)
mb_arr = et.min_bias_arr(parent_vec, lep1, lep2, lep3, lep4)
# J/Psi pt and p prompts
pt_prompt = 0.5
p_prompt = 10

print(particle[1],particle[3],'Statistics')
print('\nThe total number of events is', nEvent)
print('The total number of',particle[1],'decay to 4 muons is', nDec)
print('The total number of decays where all daughters are in LHCb acceptance',nDecAcc )
print('The total number of decays where both the daughters and the parent particle are in LHCb acceptance is', nPDAcc)
#print('The average number of decays with daughters in acceptance per event is', nDecAcc/nEvent)
print('The average number of decays with parent and daughters in acceptance per event is', nPDAcc/nEvent)

#Code Printing out percent of min bias muons above a pt threshold. 
muons_sort_pt = et.eta_acc_sort_pt(parent_vec, lep1, lep2, lep3, lep4)

# Table 1 muon Statistics
print('\nMuon pT in acc information. Children of',particle[1], 'parent particle')
table1=[['-'            ,'1st Muon', '2nd Muon', '3rd Muon', '4th Muon'], 
        ['Highest (GeV)', '{:.5f}'.format(muons_sort_pt[:,3].max()), '{:.5f}'.format(muons_sort_pt[:,2].max()), '{:.5f}'.format(muons_sort_pt[:,1].max()), '{:.5f}'.format(muons_sort_pt[:,0].max())],
        ['Lowest  (GeV)', '{:.5f}'.format(muons_sort_pt[:,3].min()), '{:.5f}'.format(muons_sort_pt[:,2].min()), '{:.5f}'.format(muons_sort_pt[:,1].min()), '{:.5f}'.format(muons_sort_pt[:,0].min())],
        ['Average (GeV)', '{:.5f}'.format(muons_sort_pt[:,3].mean()), '{:.5f}'.format(muons_sort_pt[:,2].mean()), '{:.5f}'.format(muons_sort_pt[:,1].mean()), '{:.5f}'.format(muons_sort_pt[:,0].mean())]]
print(tabulate(table1,tablefmt="fancy_grid", showindex="always"))

# Table 2 Muons meeting J/Psi Prompt 1 (pt > 0.5 GeV and p > 10 GeV) Statistics
print('\nPrompt 1 \n Statistics for',particle[1], 'daughter particles in acc meeting prompt 1 J/Psi-->4mu conditions (pt > 0.5 GeV and p > 10 GeV).')
table2 = [['-'                                              , '1st Muon', '2nd Muon', '3rd Muon', '4th Muon', 'Fraction'],
        ['All Muons pT > 0.5 GeV'                           , 'x'       , 'x'       , 'x'       , 'x'       , '{:.4f}'.format(et.mb_pt_min(mb_arr, pt_prompt))],
        ['All Muons p > 10 GeV'                            , 'x'      ,'x'        ,'x'        ,'x'        ,'{:.4f}'.format(et.mb_p_min(mb_arr, p_prompt))],
        ['All Muons pT > 0.5 GeV and p > 10 GeV'         , 'x'      ,'x'        ,'x'        ,'x'        ,'{:.4f}'.format(et.pt_p_min(mb_arr, pt_prompt, p_prompt))]]
print(tabulate(table2,tablefmt="fancy_grid", showindex="always"))

# Table 3 Muons meeting prompt 2 statistics
print('\nPrompt 2 \n Statistics for',particle[1], 'daughter particles in acc meeting prompt 2 J/Psi-->4mu conditions.')
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