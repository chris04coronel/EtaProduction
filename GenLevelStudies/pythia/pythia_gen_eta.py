# Imports
# STL Packages 
import sys
import pdb
import yaml
# Scikit Packages
import numpy as np
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
ofile_name = "EtaProduction.root"

# Read in Card File
# How to read in a card file. Look up in a tutorial.
# If you have two string can "+"
pythia = pythia8.Pythia()
pythia.readFile(eta_path + "/GenLevelStudies/pythia/" + card_file_name)

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

targ_eta_array = np.array([0]*12, dtype=np.float32)
mu1 = np.array([0]*12, dtype=np.float32)
mu2 = np.array([0]*12, dtype=np.float32)
mu3 = np.array([0]*12, dtype=np.float32)
mu4 = np.array([0]*12, dtype=np.float32)
AllMu = [mu1, mu2, mu3, mu4]

targ_eta_prime_array = np.array([0]*12, dtype=np.float32)
p_mu1 = np.array([0]*12, dtype=np.float32)
p_mu2 = np.array([0]*12, dtype=np.float32)
p_mu3 = np.array([0]*12, dtype=np.float32)
p_mu4 = np.array([0]*12, dtype=np.float32)
AllPMu = [p_mu1, p_mu2, p_mu3, p_mu4]

# Set up ROOT
file = ROOT.TFile.Open(eta_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.


# Creating Tree for Eta and its muons
tree_eta = ROOT.TTree("Tree_eta", "Tree_eta")
tree_eta.Branch('Event', evt_array, 'Event/F') # Only has 1 leaf
tree_eta.Branch('Eta', targ_eta_array, var_str)
tree_eta.Branch('Muon1', p_mu1, var_str) # Has 12 leavess. Shown on line 46
tree_eta.Branch('Muon2', p_mu2, var_str)
tree_eta.Branch('Muon3', p_mu3, var_str)
tree_eta.Branch('Muon4', p_mu4, var_str)

# Creating Tree/Branches for Eta_prime and its muons
tree_eta_prime = ROOT.TTree("Tree_eta_prime", "Tree_eta_prime")
tree_eta_prime.Branch('Event', evt_array, 'Event/F') # Only has 1 leaf
tree_eta_prime.Branch('EtaPrime', targ_eta_array, var_str)
tree_eta_prime.Branch('MuonPrime1', mu1, var_str) # Has 12 leavess. Shown on line 46
tree_eta_prime.Branch('MuonPrime2', mu2, var_str)
tree_eta_prime.Branch('MuonPrime3', mu3, var_str)
tree_eta_prime.Branch('MuonPrime4', mu4, var_str)

# Entering Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    evt_array[:] = iEvent 
    # ^Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list 
    targ_eta_list = []
    targ_particle_list = []

    targ_eta_prime_list = []
    targ_prime_particle_list = []

    # Eta Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 221: # Look for outgoing particle eta id = 221 
            targ_eta_list.append(index)
            targ_particle_list.append(particle)        
    # ^Now that we have list of desired particles we need to find out order of particles.         
    
    # Find Daughters of Eta Loop  
    for i in range(len(targ_eta_list)):
        daughter_index_list = [0]*(targ_particle_list[i].daughter2()-targ_particle_list[i].daughter1()+1)
        targ_eta_array = at.fill_array(targ_eta_array, pythia.event, targ_eta_list[i])  
        for j in range(len(daughter_index_list)):
            daughter_index_list[j] = j + targ_particle_list[i].daughter1()
            AllMu[j] = at.fill_array(AllMu[j], pythia.event, daughter_index_list[j])
        tree_eta.Fill()  

    # Eta_Prime Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 331: # Look for outgoing particle eta id = 221 
            targ_eta_prime_list.append(index)
            targ_prime_particle_list.append(particle)        
    
    # Find Daughters of Eta_Prime Loop  
    for i in range(len(targ_eta_prime_list)):
        prime_daughter_index_list = [0]*(targ_prime_particle_list[i].daughter2()-targ_prime_particle_list[i].daughter1()+1)
        targ_eta_prime_array = at.fill_array(targ_eta_prime_array, pythia.event, targ_eta_prime_list[i])  
        for j in range(len(prime_daughter_index_list)):
            prime_daughter_index_list[j] = j + targ_prime_particle_list[i].daughter1()
            AllPMu[j] = at.fill_array(AllPMu[j], pythia.event, prime_daughter_index_list[j])
        tree_eta_prime.Fill()


# pdb.set_trace() 
# pythia.stat()
# tree.Print()
file.Write() 