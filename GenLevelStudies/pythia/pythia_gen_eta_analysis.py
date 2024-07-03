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
ofile_name = "EtaProduction1.root"

# Read in Card File
# How to read in a card file. Look up in a tutorial.
# If you have two string can "+"
pythia = pythia8.Pythia()
pythia.readFile(eta_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
# Starting up pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents") 
#This is calling number of events from line 27

# Initialize Arrays
# Connect branches to python_array
# Hold information for root, but first have to create the arrays in python
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F' # 12 inputs
evt_array = np.array([0], dtype=np.float32)

# I should be looking at children particles. Either lep1, lep2, ... and possibly eta(parent) 
target_eta_array  = np.array([0]*12, dtype=np.float32)
acc_eta_array = np.array([0]*12, dtype=np.float32)
target_lep1_array = np.array([0]*12, dtype=np.float32)
target_lep2_array = np.array([0]*12, dtype=np.float32)
target_lep3_array = np.array([0]*12, dtype=np.float32)
target_lep4_array = np.array([0]*12, dtype=np.float32)
all_lep_list      =[target_lep1_array,target_lep2_array,target_lep3_array,target_lep4_array]

# Set up ROOT
file = ROOT.TFile.Open(eta_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.

tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F') # Only has 1 leaf
tree.Branch('TargetEta', target_eta_array, var_str)
tree.Branch('EtaInAcc', acc_eta_array, var_str)

tree.Branch('TargetLepton1', target_lep1_array, var_str) 
tree.Branch('TargetLepton2', target_lep2_array, var_str)
tree.Branch('TargetLepton3', target_lep3_array, var_str)
tree.Branch('TargetLepton4', target_lep4_array, var_str)



#Testing
counter1 = 0
counter2 = 0

# Entering Event Loop
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    evt_array[:] = iEvent 
    # ^Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list 
    target_eta_list = []
    target_particle_list = []
    #acc eta
    acc_eta_list = []
    acc_particle_list = []
    

    # Particle Loop
    # Desired eta particles
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 221: # Look for outgoing particle eta id = 221 
            target_eta_list.append(index)
            target_particle_list.append(particle)
            if  particle.eta() <= 5 and 2 <= particle.eta():
                acc_eta_list.append(index)
                acc_particle_list.append(particle)
                acc_eta_array = at.fill_array(acc_eta_array, pythia.event, index)
                tree.Fill()
                #create an array with attributes of eta in acceptance so i can probe that requirments are met
                
    
    # Find Target Particles Loop  
    for i in range(len(target_eta_list)):
        daughter_index_list = [0]*(target_particle_list[i].daughter2()-target_particle_list[i].daughter1()+1)
        target_eta_array = at.fill_array(target_eta_array, pythia.event, target_eta_list[i])  
        for j in range(len(daughter_index_list)):
            daughter_index_list[j] = j + target_particle_list[i].daughter1()
            all_lep_list[j] = at.fill_array(all_lep_list[j], pythia.event, daughter_index_list[j])
        tree.Fill()   
#pdb.set_trace() 
# Now that we filled our arrays we also have to fill our connected leaves to the array. 
print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write() 