# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt

import awkward as ak
import uproot3
import vector

import ROOT as root
# Personal Packages
sys.path.append(".")
import AnalysisTools_etameson as at

file = '/data/home/chris04coronel/EtaProduction/GenLevelStudies/EtaPro.root'
eta_data = uproot3.open(file)
tree = eta_data['Tree'].arrays()
tree.show()

# Create Vectors
