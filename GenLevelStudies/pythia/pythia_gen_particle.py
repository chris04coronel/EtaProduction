# Imports
# STL Packages 
import sys
import pdb
import yaml
# Scikit Packages
import numpy as np
# HEP Packages
import pythia8
#import root
import ROOT

# Hard Code
# [carfile(Must exist), rootfile(create), particle_branch(Create), pID(Must exist)]
#p = ['eta4mu.card','Eta4Mu.root', 'Eta', 221]
#p = ['etaprime4mu.card','EtaPrime4Mu.root', 'EtaPrime', 331]
p = ['jpsi4mu.card', 'JPsi4Mu.root', 'JPsi', 443]
#p = ['eta2mu2e.card', 'Eta2Mu2E.root', 'Eta_2', 221]
# Phase space meMode
#p = ['eta4mu_phsp.card','Eta4Mu_phsp.root', 'Eta', 221]

# Run 3
#p = ['eta4mu_r3.card','Eta4Mu_R3.root', 'Eta', 221]
#p = ['etaprime4mu_r3.card','EtaPrime4Mu_R3.root', 'EtaPrime', 331]

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
card_file_name = p[0]
#parent Id selector
ofile_name = p[1]

# Read in Card File
# How to read in a card file. Look up in a tutorial.
pythia = pythia8.Pythia()
pythia.readFile(eta_path + "/GenLevelStudies/pythia/CardFiles/" + card_file_name)

# Initialize Pythia
# Starting up pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents") 

# Initialize Arrays
# Connect branches to python_array
# Hold information for root, but first have to create the arrays in python
# /F tells ROOT what type of variable, F for float i.e. boolean, int, etc. 
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F' # 12 inputs
evt_array = np.array([0], dtype=np.float32)

targ_particle_array = np.array([0]*12, dtype=np.float32)
mu1 = np.array([0]*12, dtype=np.float32)
mu2 = np.array([0]*12, dtype=np.float32)
mu3 = np.array([0]*12, dtype=np.float32)
mu4 = np.array([0]*12, dtype=np.float32)
AllMu = [mu1, mu2, mu3, mu4]

# Set up ROOT
file = ROOT.TFile.Open(eta_path + "/GenLevelStudies/pythia/RootFiles/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.


# Creating Tree for particle and its muons
tree = ROOT.TTree("Tree_Particle", "Tree_Particle")
tree.Branch('Event', evt_array, 'Event/F') # Only has 1 leaf
tree.Branch(p[2], targ_particle_array, var_str)
tree.Branch('Muon1', mu1, var_str) # Has 12 leavess. Shown on line 46
tree.Branch('Muon2', mu2, var_str)
tree.Branch('Muon3', mu3, var_str)
tree.Branch('Muon4', mu4, var_str)

# Entering Event Loop
for iEvent in range(nEvent):
    if not pythia.next(): continue
    evt_array[:] = iEvent 
    # ^Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list 
    targ_eta_list = []
    targ_particle_list = []

    # Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == p[3]: # Look for outgoing particle
            targ_eta_list.append(index)
            targ_particle_list.append(particle)        
    # ^Now that we have list of desired particles we need to find out order of particles.         
    
    # Find Daughters of particle Loop  
    for i in range(len(targ_eta_list)):
        daughter_index_list = [0]*(targ_particle_list[i].daughter2()-targ_particle_list[i].daughter1()+1)
        targ_particle_array = at.fill_array(targ_particle_array, pythia.event, targ_eta_list[i])  
        for j in range(len(daughter_index_list)):
            daughter_index_list[j] = j + targ_particle_list[i].daughter1()
            AllMu[j] = at.fill_array(AllMu[j], pythia.event, daughter_index_list[j])
        tree.Fill()  

# pdb.set_trace() 
pythia.stat()
tree.Print()
file.Write() 
#pdb.set_trace() 
