# Imports
import sys
import os
import pdb

import numpy as np
import matplotlib.pyplot as plt
import EtaTools as et

import awkward as ak
import uproot
import vector
import ROOT as root

# Open the file
ifile = uproot.open("pythia/EtaProduction.root")
tree = ifile['Tree_eta'].arrays()
treeprime = ifile['Tree_eta_prime'].arrays()

# Eta and it's Muons: Four Vectors
eta_vec = vector.zip({
    'px': tree['Eta'].px,
    'py': tree['Eta'].py,
    'pz': tree['Eta'].pz,
    'e': tree['Eta'].energy
})
eta_mu_1_vec = vector.zip({
    'px': tree['Muon1'].px,
    'py': tree['Muon1'].py,
    'pz': tree['Muon1'].pz,
    'e': tree['Muon1'].energy
})
eta_mu_2_vec = vector.zip({
    'px': tree['Muon2'].px,
    'py': tree['Muon2'].py,
    'pz': tree['Muon2'].pz,
    'e': tree['Muon2'].energy
})
eta_mu_3_vec = vector.zip({
    'px': tree['Muon3'].px,
    'py': tree['Muon3'].py,
    'pz': tree['Muon3'].pz,
    'e': tree['Muon3'].energy
})
eta_mu_4_vec = vector.zip({
    'px': tree['Muon4'].px,
    'py': tree['Muon4'].py,
    'pz': tree['Muon4'].pz,
    'e': tree['Muon4'].energy
})

# Eta_Prime and it's Muons: Four Vectors
EtaPVec = vector.zip({
    'px': treeprime['EtaPrime'].px,
    'py': treeprime['EtaPrime'].py,
    'pz': treeprime['EtaPrime'].pz,
    'e': treeprime['EtaPrime'].energy
})
PMuVec1 = vector.zip({
    'px': treeprime['MuonPrime1'].px,
    'py': treeprime['MuonPrime1'].py,
    'pz': treeprime['MuonPrime1'].pz,
    'e': treeprime['MuonPrime1'].energy
})
PMuVec2 = vector.zip({
    'px': treeprime['MuonPrime2'].px,
    'py': treeprime['MuonPrime2'].py,
    'pz': treeprime['MuonPrime2'].pz,
    'e': treeprime['MuonPrime2'].energy
})
PMuVec3 = vector.zip({
    'px': treeprime['MuonPrime3'].px,
    'py': treeprime['MuonPrime3'].py,
    'pz': treeprime['MuonPrime3'].pz,
    'e': treeprime['MuonPrime3'].energy
})
PMuVec4 = vector.zip({
    'px': treeprime['MuonPrime4'].px,
    'py': treeprime['MuonPrime4'].py,
    'pz': treeprime['MuonPrime4'].pz,
    'e': treeprime['MuonPrime4'].energy
})
####################
####################
####################
####################

# Histogram 1 pT
muons_pt = et.LApt(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
et.bar_hist(muons_pt, '4Muons_pT', 'pT Among All Muons', 'Energy (GeV)', 'Number of Particles per pT', 'teal', 70)

# Histogram 2 all 4 muons 
muons_eta = et.LAeta(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
et.bar_hist(muons_eta, '4MuonsEta', 'eta Among Individual Muons', 'eta', 'Number pf particle per eta','teal',  35 )

# Histogram 3 eta spread among sibling muons
# 2d array each row is 4 sibling muons and their listed etas
sibmuons = et.TwoDEta(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
counter = 0
delta_eta_list = []
for i in range(len(sibmuons[:,1])):
    mineta = 0
    maxeta = 0
    sibs = sibmuons[i,:]
    mineta = np.min(sibs)
    maxeta = np.max(sibs)
    delta_eta = maxeta-mineta
    delta_eta_list.append(delta_eta)
    if sibs[0] <= 5 and sibs[0] >= 2 and sibs[1] <= 5 and sibs[1] >= 2 and sibs[2] <= 5 and sibs[2] >= 2 and sibs[3] <= 5 and sibs[3] >= 2:
        counter += 1
        if counter == 1:
            sibs_eta_acc = sibmuons[i,:]
        else:
            sibs_eta_acc = np.vstack((sibs_eta_acc,sibs))
et.bar_hist(delta_eta_list, '4MuonsEtaSpread', 'Eta Spread Among all eta siblings', 'Delta eta', 'Number of Particles per Delta eta','teal', 35)

# Histogram 4 
# Delta eta when all siblings are in Acceptance
delta_eta_acc_list = []
for i in range(len(sibs_eta_acc[:,1])):
    accmineta = 0
    accmaxeta = 0
    accsibs = sibs_eta_acc[i,:]
    accmineta = np.min(accsibs)
    accmaxeta = np.max(accsibs)
    acc_delta = accmaxeta - accmineta
    delta_eta_acc_list.append(acc_delta)
et.bar_hist(delta_eta_acc_list, '4MuonsAccEtaSpread', 'Eta Spread Among all eta siblings in Acceptance', 'Delta eta', 'Number of Particles per Delta eta','teal', 35)

print('The Highest pT among all muons is', muons_pt.max(),'GeV')
print('The Lowest pT among all muons is', muons_pt.min(),'GeV')
print('The average pseudorapidity spread among siblings is', np.array(delta_eta_list).mean())
print('The greatest pseudorapidity spread among siblings is', np.array(delta_eta_list).max())
print('The smallest pseudorapidity spread among siblings is', np.array(delta_eta_list).min())

print('\nThe avaerage pseudorapidity spread among siblings all in acceptance is', np.array(delta_eta_acc_list).mean())
print('The greatest pseudorapidity spread among siblings all in acceptance is', np.array(delta_eta_acc_list).max())
print('The smallest pseudorapidity spread among siblings all in acceptance is', np.array(delta_eta_acc_list).min())

################
################
################
################ 

# Prime Histogram 1' muons pT
prime_muons_pt = et.LApt(PMuVec1, PMuVec2, PMuVec3, PMuVec4)
et.bar_hist(prime_muons_pt, 'prime_4Muons_pT', 'pT Among Muons from Eta Prime', 'Energy (GeV)', 'Number of Particles per pT', 'rebeccapurple', 45)

# Prime Histogram 2' all 4 muons in acceptance
prime_muons_eta = et.LAeta(PMuVec1, PMuVec2, PMuVec3, PMuVec4)
et.bar_hist(prime_muons_eta, 'prime_4Muons_eta','Eta Among Individual Muons from Eta Prime', 'eta' ,'Number of Particles per eta', 'rebeccapurple', 35 )

# Prime Histogram 3' eta spread among sibling muons
# 2d array each row is 4 sibling muons and their listed etas
prime_sibmuons = et.TwoDEta(PMuVec1, PMuVec2, PMuVec3, PMuVec4)
counter = 0
prime_delta_eta_list = []
for i in range(len(prime_sibmuons[:,1])):
    pmineta = 0
    pmaxeta = 0
    sibs = prime_sibmuons[i,:]
    pmineta = np.min(sibs)
    pmaxeta = np.max(sibs)
    prime_delta_eta = pmaxeta-pmineta
    prime_delta_eta_list.append(prime_delta_eta)
    if sibs[0] <= 5 and sibs[0] >= 2 and sibs[1] <= 5 and sibs[1] >= 2 and sibs[2] <= 5 and sibs[2] >= 2 and sibs[3] <= 5 and sibs[3] >= 2:
        counter += 1
        if counter == 1:
            prime_sibs_eta_acc = sibmuons[i,:]
        else:
            prime_sibs_eta_acc = np.vstack((prime_sibs_eta_acc,sibs))
et.bar_hist(prime_delta_eta_list, 'prime_4muons_delta_eta', 'Greatest Delta Eta Among Sibling Muons From Eta Prime', 'Delta Eta', 'Number of Particles per Delta eta', 'rebeccapurple', 35) 

# Prime Histogram 4' 
# Delta eta when all siblings are in Acceptance
prime_delta_eta_acc_list = []
for i in range(len(prime_sibs_eta_acc[:,1])):
    pmineta = 0
    pmaxeta = 0
    pmineta = np.min(prime_sibs_eta_acc[i,:])
    pmaxeta = np.max(prime_sibs_eta_acc[i,:])
    prime_delta_eta_acc = pmaxeta-pmineta
    prime_delta_eta_acc_list.append(prime_delta_eta_acc)

et.bar_hist(prime_delta_eta_acc_list, 'prime_4muons_acc_delta_eta', 'Greatest Delta Eta Among Sibling Muons (all in acceptance) From Eta Prime', 'Delta Eta', 'Number of Particles per Delta eta', 'rebeccapurple', 35)

############
############
############
#New Histogram (All Muons in Acceptance)
all_muons_acc_pt = et.eta_acc_LApt(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
et.bar_hist(all_muons_acc_pt, '4muons_acc_pT', 'Muons (all in acceptance) pT', 'pT in GeV','Number of Particles per pT', 'teal', 65 )


# # New Histogram
all_prime_muons_acc_pt = et.eta_acc_LApt(PMuVec1, PMuVec2, PMuVec3, PMuVec4)
et.bar_hist(all_prime_muons_acc_pt, 'prime_4muons_acc_pT', 'Eta_Prime Muons (all in acceptance) pT','pT in GeV', 'Number of Particles per pT', 'rebeccapurple', 65)

#NewHistogram
fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
#counts, edges, bars = plt.hist([all_muons_acc_pt,all_prime_muons_acc_pt], color=['teal','rebeccapurple'], bins=35, label=['Eta Muons', 'Eta Prime Muons'])
#plt.bar_label(bars)
plt.hist([all_muons_acc_pt,all_prime_muons_acc_pt], color=['teal','rebeccapurple'], bins=35, label=['Eta Muons', 'Eta Prime Muons'], alpha=0.5, range=[0,2])
axs.xaxis.set_tick_params(pad = 20) 
axs.yaxis.set_tick_params(pad = 20) 
plt.title('pT of Muons in Acceptance Decayed from Eta and Eta_Prime')
plt.legend(loc='upper right')
plt.xlabel('pT in GeV')
plt.ylabel('Number of Particles per pT')
plt.savefig('pythia/Histograms/E_EP_muons_acc_pT')

print('The Highest pT among all muons from eta_prime is', prime_muons_pt.max(),'GeV')
print('The Lowest pT among all muons from eta_pime is', prime_muons_pt.min(),'GeV')
print('The average pseudorapidity spread among siblings from eta prime is', np.array(prime_delta_eta_list).mean())
print('The greatest pseudorapidity spread among siblings from eta prime is', np.array(prime_delta_eta_list).max())
print('The smallest pseudorapidity spread among siblings from eta prime is', np.array(prime_delta_eta_list).min())

print('\nThe avaerage pseudorapidity spread among siblings all in acceptance from eta prime is', np.array(prime_delta_eta_acc_list).mean())
print('The greatest pseudorapidity spread among siblings all in acceptance from eta prime is', np.array(prime_delta_eta_acc_list).max())
print('The smallest pseudorapidity spread among siblings all in acceptance from eta prime is', np.array(prime_delta_eta_acc_list).min())

print('\nAll Muons in acceptance decayed from both eta and eta prime')
print('The highest pT among all muons from eta is', np.array(all_muons_acc_pt).max(),'GeV')
print('The highest pT among all muons from eta prime is', np.array(all_prime_muons_acc_pt).max(),'GeV')
print('The lowest pT among all muons from eta is', np.array(all_muons_acc_pt).min(),'GeV')
print('The lowest pT among all muons from eta prime is', np.array(all_prime_muons_acc_pt).min(),'GeV')
print('The average pT among all muons from eta is', np.array(all_muons_acc_pt).mean(),'GeV')
print('The average pT among all muons from eta prime is', np.array(all_prime_muons_acc_pt).mean(),'GeV')

muon_sorted_pt = et.eta_acc_sort_pt(eta_mu_1_vec, eta_mu_2_vec, eta_mu_3_vec, eta_mu_4_vec)
muon_prime_sorted_pt = et.eta_acc_sort_pt(PMuVec1, PMuVec2, PMuVec3, PMuVec4)

print('\n information about muons decayed from eta sorted pt')
print('The highest pT among the largets pt eta muons is', muon_sorted_pt[:,3].max(),'GeV')
print('The highest pT among the 2nd eta muons is', muon_sorted_pt[:,2].max(),'GeV')
print('The highest pT among the 3rd eta muons is', muon_sorted_pt[:,1].max(),'GeV')
print('The highest pT among eta 4th muons is', muon_sorted_pt[:,0].max(),'GeV')
print('The lowest pT among the largets pt eta muons is', muon_sorted_pt[:,3].min(),'GeV')
print('The lowest pT among the 2nd eta muons is', muon_sorted_pt[:,2].min(),'GeV')
print('The lowest pT among the 3rd eta muons is', muon_sorted_pt[:,1].min(),'GeV')
print('The lowest pT among the 4th eta muons is', muon_sorted_pt[:,0].min(),'GeV')
print('The average pT among the largets pt eta muons is', muon_sorted_pt[:,3].mean() ,'GeV')
print('The average pT among the 2nd eta muons is', muon_sorted_pt[:,2].mean(),'GeV')
print('The average pT among the 3rd eta muons is', muon_sorted_pt[:,1].mean(),'GeV')
print('The average pT among the 4th eta muons is', muon_sorted_pt[:,0].mean(),'GeV')

print('\nInformation about muons decayed from eta_prime sorted pt')
print('The highest pT among the largets pt eta_prime muons is', muon_prime_sorted_pt[:,3].max(),'GeV')
print('The highest pT among the 2nd eta_prime muons is', muon_prime_sorted_pt[:,2].max(),'GeV')
print('The highest pT among the 3rd eta_prime muons is', muon_prime_sorted_pt[:,1].max(),'GeV')
print('The highest pT among eta_prime 4th muons is', muon_prime_sorted_pt[:,0].max(),'GeV')
print('The lowest pT among the largets pt eta_prime muons is', muon_prime_sorted_pt[:,3].min(),'GeV')
print('The lowest pT among the 2nd eta_prime muons is', muon_prime_sorted_pt[:,2].min(),'GeV')
print('The lowest pT among the 3rd eta_prime muons is', muon_prime_sorted_pt[:,1].min(),'GeV')
print('The lowest pT among the 4th eta_prime muons is', muon_prime_sorted_pt[:,0].min(),'GeV')
print('The average pT among the largets pt eta_prime muons is', muon_prime_sorted_pt[:,3].mean() ,'GeV')
print('The average pT among the 2nd eta_prime muons is', muon_prime_sorted_pt[:,2].mean(),'GeV')
print('The average pT among the 3rd eta_prime muons is', muon_prime_sorted_pt[:,1].mean(),'GeV')
print('The average pT among the 4th eta_prime muons is', muon_prime_sorted_pt[:,0].mean(),'GeV')
#pdb.set_trace()