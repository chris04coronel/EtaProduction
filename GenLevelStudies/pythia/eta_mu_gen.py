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
ofile_name = "EtaMu.root"

# Read in Card File
pythia = pythia8.Pythia()
pythia.readFile(eta_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents") 

# Initialize array
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F' 
evt_array = np.array([0], dtype=np.float32)

targ_eta_array  = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(eta_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.

tree_eta = ROOT.TTree("Tree_eta", "Tree_eta")
tree_eta.Branch('Event', evt_array, 'Event/I')
tree_eta.Branch('Eta_Kin', targ_eta_array, var_str)
tree_eta.Branch('EtaPerEvent', eta, 'TEPE/I:TEAPE/I') #TotalEtaPerEvent(TEPE), TotalEtaAcceptancePerEvent(TEAPE)
tree_eta.Branch('EtaNA  ')

# Total Counters

# Entering Event Loop
for iEvent in range(nEvent):
     if not pythia.next():
        continue
    evt_array[:] = iEvent 
    # ^Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list 

    # Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 221: 
            targ_eta_list.append(index)
            targ_eta_particle_list.append(particle)

   

# pdb.set_trace() 
pythia.stat()
#tree.Print()
file.Write()  