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
targ_eta_prime_array = np.array([0]*12, dtype=np.float32)
eta = np.array([0]*2, dtype=np.int32)
eta_prime = np.array([0]*2, dtype=np.int32)
eta_muons = np.array([0]*1, dtype=np.int32)
eta_prime_muons = np.array([0]*1, dtype=np.int32)
# Set up ROOT
file = ROOT.TFile.Open(eta_path + "/GenLevelStudies/pythia/" + ofile_name,
                       "RECREATE")
                       # RECREATE if the file isnt there it will make it, if it is, it will overwrite.

tree = ROOT.TTree("Tree", "Tree")
tree.Branch('Event', evt_array, 'Event/I') # Only has 1 leaf
tree.Branch('Eta_Kin', targ_eta_array, var_str)
tree.Branch('Eta_Prime_Kin', targ_eta_prime_array, var_str)
tree.Branch('EtaPerEvent', eta, 'TEPE/I:TEAPE/I')
tree.Branch('EtaPrimePerEvent', eta_prime, 'TEPPE/I:TEPAPE/I')
tree.Branch('EtaMuons', eta_muons, 'EtaMuonsAcc/I')
tree.Branch('EtaPrimeMuons', eta_prime_muons, 'EtaPrimeMuonsAcc')
# Total counters
total_eta = 0
total_eta_acc = 0
total_eta_prime = 0
total_eta_prime_acc = 0
total_eta_muons = 0 
total_eta_prime_muons = 0 
# Entering Event Loop
for iEvent in range(nEvent):
    eta[0] = 0
    eta[1] = 0
    eta_prime[0] = 0
    eta_prime[1] = 0
    eta_muons[0] = 0
    eta_prime_muons[0] = 0
    if not pythia.next():
        continue
    evt_array[:] = iEvent 
    # ^Standard way to fill root arrays in python. Not very pythonic. 
    # Create empty list 
    targ_eta_list = []
    targ_eta_particle_list = []
    targ_eta_prime_list = []
    targ_eta_prime_particle_list = []

    targ_lep_list = []
    targ_lep_particle = []

    # Particle Loop
    # Desired eta particles
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 221: 
            eta[0] += 1
            total_eta += 1
            targ_eta_list.append(index)
            targ_eta_particle_list.append(particle)
            targ_eta_array = at.fill_array(targ_eta_array, pythia.event, index)

            if particle.eta() <= 5 and particle.eta() >= 2:
                eta[1] += 1
                total_eta_acc += 1

    # Desired eta_prime particles
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 331: 
            eta_prime[0] += 1
            total_eta_prime += 1
            targ_eta_prime_list.append(index)
            targ_eta_prime_particle_list.append(particle)
            targ_eta_prime_array = at.fill_array(targ_eta_prime_array, pythia.event, index)

            if particle.eta() <= 5 and particle.eta() >= 2:
                eta_prime[1] += 1
                total_eta_prime_acc += 1
    
    # Desired muons 
    for index, particle in enumerate(pythia.event):
        if particle.idAbs() == 13 and particle.eta() <= 5 and particle.eta() >= 2:
            targ_lep_list.append(index)
            targ_lep_particle.append(particle)
    # Muons: From eta parent
    for i in range(len(targ_lep_list)):
        if pythia.event[pythia.event[targ_lep_list[i]].mother1()].id() == 221 and pythia.event[pythia.event[targ_lep_list[i]].mother1()].eta() >= 2 and pythia.event[pythia.event[targ_lep_list[i]].mother1()].eta() <= 5:
            eta_muons[0] += 1
            total_eta_muons +=1
    # Muons: From eta_prime parent
    for i in range(len(targ_lep_list)):
        if pythia.event[pythia.event[targ_lep_list[i]].mother1()].id() == 331 and pythia.event[pythia.event[targ_lep_list[i]].mother1()].eta() >= 2 and pythia.event[pythia.event[targ_lep_list[i]].mother1()].eta() <= 5:
            eta_prime_muons[0] += 1
            total_eta_prime_muons += 1
    tree.Fill()

print("\n The total number of events for this script is ",nEvent,".\n")
print("There are a total of", total_eta, "eta produced.")
print("There are a total of", total_eta_acc, "eta in acceptance produced.")
print("There are a total of", total_eta_muons, "muons in acceptance produced from eta(also in acceptance).")
print(f"For every eta in acceptance there are {total_eta_muons/total_eta_acc:.2f} muons also in acceptance.")
print(f"For every {total_eta/total_eta_acc:.2f} eta produced one will be in the acceptance. \n")


print("There are a total of", total_eta_prime, "eta_primes produced.")
print("There are a total of", total_eta_prime_acc, "eta_primes in acceptance produced.")
print("There are a total of", total_eta_prime_muons, "muons in acceptance produced from eta_prime(also in acceptance).")
print(f"For every eta_prime in acceptance there are {total_eta_prime_muons/total_eta_prime_acc:.2f} muons also in acceptance.")
print(f"For every {total_eta_prime/total_eta_prime_acc:.2f} eta_prime produced one will be in the acceptance. \n")

print(f"The ratio of eta to eta_prime particle in acceptance is {total_eta_acc/total_eta_prime_acc:.2f}.")
print(f"The ratio of eta muons to eta_prime muons in acceptance is {total_eta_muons/total_eta_prime_muons:.2f}.")
print("From this ratio just above I conclude there are more muons in acceptance from the eta than eta_prime")
 
 
# pdb.set_trace() 
# pythia.stat()
#tree.Print()
file.Write()      
#tree.Draw("TotalEtaPerEvent:TotalEtaAccPerEvent")
#tree.Scan('EtaPerEvent.TEPE:EtaPerEvent.TEAPE:EtaMuons.EtaMuonsAcc')