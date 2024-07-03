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
import AnalysisTools as at

# Setting up Path
ww_path = at.find_WW_path()

#  Set up makefile configuration.
# pythia tells us to include this
cfg = open(ww_path + "/GenLevelStudies/pythia/Makefile.inc")
lib = "/Applications/pythia8310/lib"
for line in cfg:
    if line.startswith("PREFIX_LIB="): lib = line[11:-1]; break
sys.path.insert(0, lib)

# Inputs
card_file_name = "WjetsProduction.cmnd"
ofile_name = "WjetsProduction.root"

# Read in Card File
# How to read in a card file. Look up in a tutorial.
# If you have two string can "+" 
pythia = pythia8.Pythia()
pythia.readFile(ww_path + "/GenLevelStudies/pythia/" + card_file_name)

# Initialize Pythia
# Starting up pythia
pythia.init()
nEvent = pythia.mode("Main:numberOfEvents")

# Initialize Arrays
# Connect branches to python_array
# Hold information for root, but first have to create the arrays in python
var_str = 'px/F:py/F:pz/F:pT/F:p/F:eta/F:energy/F:phi/F:m0/F:id/F:charge/F:status/F' # 12 inputs
evt_array = np.array([0], dtype=np.float32)
# Look at the number if inputs in line 44. This is why line below has "12".
# I should be looking at daughter/children particles. Either lep1, lep2, ...,lepn and possibly eta 
target_lepton_array = np.array([0]*12, dtype=np.float32)
target_jet_array = np.array([0]*12, dtype=np.float32)

# Set up ROOT
file = ROOT.TFile.Open(ww_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.
#This is the first instance of creating a tree
#(Name, Title)
#create our branches. It can have a variable number of leaves. 
tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/F') # Only has 1 leaf
tree.Branch('TargetLepton', target_lepton_array, var_str) # Has 12 leafs. Shown on line 44
tree.Branch('TargetJet', target_jet_array, var_str)
#what is the name of the first leaf in the target lepton branch? It is px can you find it? Hint: line 44


# Testing
counter1 = 0
counter2 = 0

# Entering Event Loop
# Go through pythia events 1 by 1.
for iEvent in range(nEvent):
    if not pythia.next():
        continue
    evt_array[:] = iEvent # Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list
    target_indices_list = []
    # Particle Loop
    for index, particle in enumerate(pythia.event):
        if particle.statusAbs() == 23: # Look for outgoing particles of hard process
            target_indices_list.append(index)
   
    # Now that we have list of desired particles we need to find out order of particles.
    # Can gind muon particles, but then below can specify if we're looking for only muons whcih came from eta particles. 
    # Once have found eta, in the next loop can loop over decay products. 
    # Gabe This: there is no right way to find particles, explore. 
    # Find Target Particles Loop
    for itarget in target_indices_list:
        target_particle = pythia.event[itarget]
        itarget_mother = target_particle.mother1()
        if (pythia.event[itarget_mother].idAbs() == 24):
            if target_particle.idAbs() in [11, 13]: # The braket here is 11 or 13 not 11-13.
                target_lepton_array = at.fill_array(target_lepton_array, pythia.event, itarget)
        else:
            target_jet_array = at.fill_array(target_jet_array, pythia.event, itarget)
    tree.Fill()
# Nopw that we filled our arrays we also have to fill our connected leaves to the array. 
#?? Look up tree.fill in root
print(f'Counter1: {counter1}')
print(f'Counter2: {counter2}')
pythia.stat()
tree.Print()
file.Write() #writes the file. 
