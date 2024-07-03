# Imports
# STL Packages 
import sys
import pdb
import yaml
# Scikit Packages
import numpy as np
import matplotlib.pyplot as plt
# HEP Packages
import pythia8
import ROOT

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

# Inputs 
card_file_name = "etameson.card"
ofile_name = "EtaPro.root"

# Read in Card File
# How to read in a card file. Look up in a tutorial.
# If you have two string can "+"
pythia = pythia8.Pythia()
pythia.readFile(eta_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents") 

# Initialize array
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F' 
evt_array = np.array([0], dtype=np.float32)


targ_eta_array  = np.array([0]*12, dtype=np.float32)
acc_eta_array = np.array([0]*12, dtype=np.float32)
# targ_lep1_array = np.array([0]*12, dtype=np.float32)
# targ_lep2_array = np.array([0]*12, dtype=np.float32)
# targ_lep3_array = np.array([0]*12, dtype=np.float32)
# targ_lep4_array = np.array([0]*12, dtype=np.float32)
# all_lep_list    =[targ_lep1_array,targ_lep2_array,targ_lep3_array,targ_lep4_array]
eta = np.array([0]*2, dtype=np.float32)
muon = 0
# Set up ROOT
file = ROOT.TFile.Open(eta_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.

tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/I') # Only has 1 leaf
tree.Branch('Eta', targ_eta_array, var_str)
#tree.Branch('EtaInAcc', acc_eta_array, var_str)
tree.Branch('EtaPerEvent', eta, 'TotalEtaPerEvent/I:TotalEtaAccPerEvent/I')
tree.Branch('MuonsAccPerEvent', muon,'MuonsAccPerEvent/I')


# tree.Branch('TargetLepton1', targ_lep1_array, var_str) 
# tree.Branch('TargetLepton2', targ_lep2_array, var_str)
# tree.Branch('TargetLepton3', targ_lep3_array, var_str)
# tree.Branch('TargetLepton4', targ_lep4_array, var_str)

# Entering Event Loop
for iEvent in range(nEvent):
    eta[0]=0
    eta[1]=0
    muon = 0
    if not pythia.next():
        continue
    evt_array[:] = iEvent 
    # ^Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list 
    targ_eta_list = []
    targ_particle_list = []
    targ_lep_list = []
    targ_lep_particle = []

    # Particle Loop
    # Desired eta particles
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 221: 
            eta[0] += 1
            targ_eta_list.append(index)
            targ_particle_list.append(particle)
            targ_eta_array = at.fill_array(targ_eta_array, pythia.event, index)

            if particle.eta() <= 5 and particle.eta() >= 2:
                eta[1] += 1
    
    # Desired lep(muon) particles
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 13 and particle.eta() <= 5 and particle.eta() >= 2:
            targ_lep_list.append(index)
            targ_lep_particle.append(particle)
    for i in range(len(targ_lep_list)):
        if pythia.event[pythia.event[targ_lep_list[i]].mother1()].id() == 221:
            muon += 1

    tree.Fill()
#pdb.set_trace() 
# pythia.stat()

tree.Print()
file.Write()      
#tree.Draw("TotalEtaPerEvent:TotalEtaAccPerEvent")
tree.Scan()