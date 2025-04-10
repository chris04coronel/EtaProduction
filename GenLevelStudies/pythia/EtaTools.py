import numpy as np
import matplotlib.pyplot as plt
import uproot
import vector
import time
import math
from scipy.optimize import curve_fit



def LApt(vec1, vec2, vec3, vec4):
    '''Creates a 1D array made from 4 independent arrays
    
    Take 4 four "vectors"  find their pT 
    Create a long 1D array that concantenates 4 four vector one after another 
    in the order they are put in. The goal is to extract the pT form the given arrays. 

    vec1 (vector): A numpy array to begin the Long array 
    vec2 (vector): A numpy array to be concatenated behind array1
    vec3 (vector): A numpy array to be concatenated behind array1 + array2
    vec4 (vector): A numpy array to be concatenated behind array1 + array2 + array 3
    '''
    import vector
    LongList =  [vec1.pt, vec2.pt, vec3.pt, vec4.pt]
    for i in range(len(LongList)):
        if i == 0:
            LongArray = np.array(LongList[0])
        else:
            LongArray = np.append(LongArray, np.array(LongList[i]))
    return(LongArray)

def LAeta(vec1, vec2, vec3, vec4):
    import vector
    LongList =  [vec1.eta, vec2.eta, vec3.eta, vec4.eta]
    for i in range(len(LongList)):
        if i == 0:
            LongArray = np.array(LongList[0])
        else:
            LongArray = np.append(LongArray, np.array(LongList[i]))
    return(LongArray)

def TwoDEta(vec1, vec2, vec3, vec4):
    TwoDimArray = np.transpose([np.array(vec1.eta),np.array(vec2.eta),np.array(vec3.eta),np.array(vec4.eta)])
    return(TwoDimArray)

def TwoDEtaAcc(vec1, vec2, vec3, vec4):
    counter = 0
    ndarr1 = np.transpose([np.array(vec1.eta),np.array(vec2.eta),np.array(vec3.eta),np.array(vec4.eta)])
    for i in range(len(ndarr1[:,1])):
        if ndarr1[i,0] <=5 and ndarr1[i,0] >=2 and ndarr1[i,1] <=5 and ndarr1[i,1] >=2 and ndarr1[i,2] <=5 and ndarr1[i,2] >=2 and ndarr1[i,3] <=5 and ndarr1[i,3] >=2:
            counter +=1
            if counter == 1:
                nd_eta_acc_arr1 = ndarr1[i,:]
            else:
                nd_eta_acc_arr1 = np.vstack((nd_eta_acc_arr1, ndarr1[i,:]))
    return(nd_eta_acc_arr1)

def PDEtaAcc(par1, vec1, vec2, vec3, vec4):
    ''' This function serves to determine the total amount of decays where both the parent and the daughter particles 
    are in the LHCb 2 < eta < 5 acceptance. Will return a number so set a variable equal to this equation.
    '''
    counter = 0
    for i in range(len(par1)):
        if par1[i].eta <=5.298 and par1[i].eta >=1.596 and vec1[i].eta <=5.298 and vec1[i].eta >=1.596 and vec2[i].eta <=5.298 and vec2[i].eta >=1.596 and vec3[i].eta <=5.298 and vec3[i].eta >=1.596 and vec4[i].eta <=5.298 and vec4[i].eta >=1.596:
            counter +=1
    return(counter)


def eta_acc_LApt(vec1, vec2, vec3, vec4):
    'takes four 4vectors cuts off those which dont meet eta(pseudorapidity) requierments and makes one list'
    import vector
    import time
    start_time=time.time()
    vec1_acc_pt = []
    vec2_acc_pt = [] 
    vec3_acc_pt = [] 
    vec4_acc_pt = []
    
    for i in range(len(vec1)):
        if vec1[i].eta >=2 and vec1[i].eta <= 5 and vec2[i].eta >=2 and vec2[i].eta <= 5 and vec3[i].eta >=2 and vec3[i].eta <=5 and vec4[i].eta >=2 and vec4[i].eta <= 5:
            vec1_acc_pt.append(vec1[i].pt)
            vec2_acc_pt.append(vec2[i].pt)
            vec3_acc_pt.append(vec3[i].pt)
            vec4_acc_pt.append(vec4[i].pt)
    vec_list_eta_acc_pt = vec1_acc_pt + vec2_acc_pt + vec3_acc_pt + vec4_acc_pt
    print("My eta acc function took", time.time() - start_time, "seconds to run")
    return(vec_list_eta_acc_pt)

def mysort(vec1, vec2, vec3, vec4):
    '''
    Takes in a vector class specific attribute like pt, eta, p, e, etc and sort them from smallest to largest
    '''
    import vector
    TwoDimArray = np.transpose([np.array(vec1),np.array(vec2),np.array(vec3),np.array(vec4)])
    for i in range(len(TwoDimArray[:,0])):
        if i == 0:
            sort_ndarray = np.sort(TwoDimArray[i,:])
        else:
            sort_ndarray = np.vstack((sort_ndarray, np.sort(TwoDimArray[i,:])))
    return(sort_ndarray)


def eta_acc_sort_pt(par1,vec1, vec2, vec3, vec4):
    'looks for all 4 muons and parent in LHCb acceptance and then sorts them by ascending muon pT returns nx4 array'
    start_time=time.time()
    counter = 0
    for i in range(len(vec1)):
        if vec1[i].eta >=1.596 and vec1[i].eta <= 5.298 and vec2[i].eta >=1.596 and vec2[i].eta <= 5.298 and vec3[i].eta >=1.596 and vec3[i].eta <=5.298 and vec4[i].eta >=1.596 and vec4[i].eta <= 5.298 and par1[i].eta <= 5.298 and par1[i].eta >=1.596:
            if counter == 0:
                counter += 1
                arr1 =np.array([vec1[i].pt, vec2[i].pt, vec3[i].pt, vec4[i].pt])
                sort_arr1 = np.sort(arr1)
            else:
                arr1 =np.array([vec1[i].pt, vec2[i].pt, vec3[i].pt, vec4[i].pt])
                sort_arr1 = np.vstack((sort_arr1, np.sort(arr1)))
    print("My eta acc pT sorting function took", '{:.5f}'.format(time.time() - start_time), "seconds to run")
    return(sort_arr1)

def eta_acc_sort_p(par1,vec1, vec2, vec3, vec4):
    'looks for all 4 muons and parent in LHCb acceptance and then sorts them by ascending muon pT returns nx4 array'
    start_time=time.time()
    counter = 0
    for i in range(len(vec1)):
        if vec1[i].eta >=2 and vec1[i].eta <= 5 and vec2[i].eta >=2 and vec2[i].eta <= 5 and vec3[i].eta >=2 and vec3[i].eta <=5 and vec4[i].eta >=2 and vec4[i].eta <= 5 and par1[i].eta <= 5 and par1[i].eta >=2:
            if counter == 0:
                counter += 1
                arr1 =np.array([vec1[i].p, vec2[i].p, vec3[i].p, vec4[i].p])
                sort_arr1 = np.sort(arr1)
            else:
                arr1 =np.array([vec1[i].p, vec2[i].p, vec3[i].p, vec4[i].p])
                sort_arr1 = np.vstack((sort_arr1, np.sort(arr1)))
    print("My eta acc pT sorting function took", '{:.5f}'.format(time.time() - start_time), "seconds to run")
    return(sort_arr1) 

def pT_check1(arr1, pt1):
    'takes in an nx1 dimensional array. Checks percentage'
    counter = 0
    for i in range(len(arr1)):
        if arr1[i] >= pt1:
            counter += 1
    print('\n{:.2f}'.format(counter/len(arr1)*100),'%','of muons have momentum >=', pt1,'Gev')

def pT_check12(arr1, arr2, pt1, pt2):
    'takes in two nx1 dimensional array. Checks if each meet pT lsiten'
    counter = 0
    for i in range(len(arr1)):
        if arr1[i] >= pt1 and arr2[i] >= pt2:
            counter += 1
    return('{:.2f}'.format(counter/len(arr1)*100))
    #print('{:.2f}'.format(counter/len(arr1)*100),'%')

def pT_check123(arr1, arr2, arr3, pt1, pt2, pt3):
    'takes in an nx4 dimensional array. Checks '
    counter = 0
    for i in range(len(arr1)):
        if arr1[i] >= pt1 and arr2[i] >= pt2 and arr3[i] >= pt3:
            counter += 1
    print(counter/len(arr1))

def pT_check1234(arr1, arr2, arr3, arr4, pt1, pt2, pt3, pt4):
    '''takes in 4 arrays and checks that if tranverse momentum 
    in each array satisfy transverse momenutm in the argument  '''
    counter = 0
    for i in range(len(arr1)):
        if arr1[i] >= pt1 and arr2[i] >= pt2 and arr3[i] >= pt3 and arr4[i] >= pt4:
            counter += 1       
    return(counter/len(arr1))

def pT_min(arr1, arr2, arr3, arr4, low_pt):
    '''takes in 4 arrays and checks that if tranverse momentum 
    in each array satisfy transverse momentum minimum requirement in GeV. returns percentage that does.  '''
    counter = 0
    for i in range(len(arr1)):
        if arr1[i] >= pt1 and arr2[i] >= pt2 and arr3[i] >= pt3 and arr4[i] >= pt4:
            counter += 1
    return('{:.2f}'.format((counter/len(arr1))*100))
#####
#-------------------------------------------------------------------------------------
#####
def LHCb_acc(par1, vec1, vec2, vec3, vec4): #use to be min_bias_arr
    '''Takes in 4 awkward array vectors that contain kinematics for each variable. Then
    it will cycle through and delete elements where not all are within the LHCb acceptance
    of 1.596 < eta <5.298 It will return a nX4 vector array.
    '''
    nvec1 = np.array([])
    nvec2 = np.array([])
    nvec3 = np.array([])
    nvec4 = np.array([])
    for i in range(len(par1)):
        if par1[i].eta > 1.596 and par1[i].eta < 5.298 and vec1[i].eta > 1.596 and vec1[i].eta < 5.298 and vec2[i].eta > 1.596 and vec2[i].eta < 5.298 and vec3[i].eta > 1.596 and vec3[i].eta < 5.298 and vec4[i].eta > 1.596 and vec4[i].eta < 5.298:
            nvec1 = np.append(nvec1,vec1[i])
            nvec2 = np.append(nvec2,vec2[i])
            nvec3 = np.append(nvec3,vec3[i])
            nvec4 = np.append(nvec4,vec4[i])
    new_vec_array = np.transpose([nvec1, nvec2, nvec3, nvec4])
    return(new_vec_array)

def min_bias_arr2(par1, vec1, vec2, vec3, vec4):
    '''Takes in 4 awkward array vectors that contain kinematics fir each variable. Then
    it will cycle through and delete elements where not all are within the LHCb acceptance
    of 2 < eta <5. It will return a nX4 vector array.
    '''
    import vector
    count = 0
    for i in range(len(par1)):
        if par1[i-count].eta > 1.596 and par1[i-count].eta < 5.298 and vec1[i-count].eta > 1.5962 and vec1[i-count].eta < 5.298 and vec2[i-count].eta > 1.596 and vec2[i-count].eta < 5.298 and vec3[i-count].eta > 1.596 and vec3[i-count].eta < 5.298 and vec4[i-count].eta > 1.596 and vec4[i-count].eta < 5.298:
            print('the index is', i)
            print('par eta is', par1[i-count].eta)
            print('lep1 eta is ', vec1[i-count].eta)
            print('lep2 eta is ', vec2[i-count].eta)
            print('lep3 eta is ', vec3[i-count].eta)
            print('lep4 eta is ', vec4[i-count].eta)

            vec1 = np.delete(vec1, i-count)
            vec2 = np.delete(vec2, i-count)
            vec3 = np.delete(vec3, i-count)
            vec4 = np.delete(vec4, i-count)
            count +=1
    new_vec_array = np.transpose([vec1, vec2, vec3, vec4])
    return(new_vec_array)

def pt_p_min(vec_darr,min_pt, min_p ):
    '''Takes in an nx4 vectory array. Prefferec it has already passed the min bias test to save on computation time. 
    This it will check that every muon candidate passes the min pt and min p. It will return fractional value of decays 
    that pass this filter. 
    '''
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,0].pt >= min_pt and vec_darr[i,0].p >= min_p and vec_darr[i,1].pt >= min_pt and vec_darr[i,1].p >= min_p and vec_darr[i,2].pt >= min_pt and vec_darr[i,2].p >= min_p and vec_darr[i,3].pt >= min_pt and vec_darr[i,3].p >= min_p:
            counter += 1
    return(counter/len(vec_darr[:,0]))
    
def mb_pt_min(vec_darr, pt_min):
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,0].pt >= pt_min and vec_darr[i,1].pt >= pt_min and vec_darr[i,2].pt >= pt_min and vec_darr[i,3].pt >= pt_min:
            counter += 1
    return(counter/len(vec_darr[:,0]))

def mb_p_min(vec_darr, p_min):
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,0].p >= p_min and vec_darr[i,1].p >= p_min and vec_darr[i,2].p >= p_min and vec_darr[i,3].p >= p_min:
            counter += 1
    return(counter/len(vec_darr[:,0]))

def sort_pt(narr1):
    'Take 4 four vectors and sorts them in a nx4. Where n is the length of the vectors (Assuming)'
    'they are all the same length. the the 1st entry is the smallest pT in ascending order'
    for i in range(len(narr1[:,0])):
        if i == 0:
            arr1 = np.array([narr1[i,0].pt, narr1[i,1].pt, narr1[i,2].pt, narr1[i,3].pt])
            sort_arr1 = np.sort(arr1)
        else:
            arr1 = np.array([narr1[i,0].pt, narr1[i,1].pt, narr1[i,2].pt, narr1[i,3].pt])
            sort_arr1 = np.vstack((sort_arr1, np.sort(arr1)))
    return(sort_arr1)

def mb_pt_check1(vec_darr, pt_cut1):
    '''
    Will take in an array that has already passed the min bias it will check the lepton
    with the highest pT passing a certain threshhold and keep them. Then it will check to
    see if the muon with the second highest pT also passes that thresh hold and return pT ONLY narray 
    where both the largest and second largeest pT muosn meet the threshold. sorted
    '''
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,3] > pt_cut1:
            counter +=1
            if counter == 1:
                new_pt_arr = np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])
            else:
                new_pt_arr = np.vstack((new_pt_arr,np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])))
    print("Fraction is", counter/len(vec_darr[:,0]))
    return(new_pt_arr)

def mb_pt_check12(vec_darr, pt_cut1, pt_cut2):
    '''
    Will take in an array that has already passed the min bias it will check the lepton
    with the highest pT passing a certain threshhold and keep them. Then it will check to
    see if the muon with the second highest pT also passes that thresh hold and return pT ONLY narray 
    where both the largest and second largeest pT muosn meet the threshold. sorted
    '''
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,3] > pt_cut1 and vec_darr[i,2] > pt_cut2:
            counter +=1
            if counter == 1:
                new_pt_arr = np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])
            else:
                new_pt_arr = np.vstack((new_pt_arr,np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])))
    print("Fraction is", counter/len(vec_darr[:,0]))
    return(new_pt_arr)

def mb_pt_check123(vec_darr, pt_cut1, pt_cut2):
    '''
    Will take in an array that has already passed the min bias it will check the lepton
    with the highest pT passing a certain threshhold and keep them. Then it will check to
    see if the muon with the second highest pT also passes that thresh hold and return pT ONLY narray 
    where both the largest and second largeest pT muosn meet the threshold. sorted
    '''
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,3] > pt_cut1 and vec_darr[i,2] > pt_cut2 and vec_darr[i,1] > pt_cut2:
            counter +=1
            if counter == 1:
                new_pt_arr = np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])
            else:
                new_pt_arr = np.vstack((new_pt_arr,np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])))
    print("Fraction is", counter/len(vec_darr[:,0]))
    return(new_pt_arr)

def mb_pt_check1234(vec_darr, pt_cut1, pt_cut2):
    '''
    Neeeds to already be sorted
    Will take in an array that has already passed the min bias it will check the lepton
    with the highest pT passing a certain threshhold and keep them. Then it will check to
    see if the muon with the second highest pT also passes that thresh hold and return pT ONLY narray 
    where both the largest and second largeest pT muosn meet the threshold. sorted
    '''
    counter = 0
    for i in range(len(vec_darr[:,0])):
        if vec_darr[i,3] > pt_cut1 and vec_darr[i,2] > pt_cut2 and vec_darr[i,1] > pt_cut2 and vec_darr[i,0] > pt_cut2:
            counter +=1
            if counter == 1:
                new_pt_arr = np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])
            else:
                new_pt_arr = np.vstack((new_pt_arr,np.array([vec_darr[i,0], vec_darr[i,1], vec_darr[i,2], vec_darr[i,3]])))
    print("Fraction is", counter/len(vec_darr[:,0]))
    return(new_pt_arr)
    

#####
#-------------------------------------------------------------------------------------
#####
def pT_nx4_min(arr1, low_pt):
    '''takes in nx4 array that has already passed the min bias LHCb acceptance and checks that the
    pT of each muon passes this low_pt threshhold.
    '''
    counter = 0
    for i in range(len(arr1[:,0])):
        if arr1[i,0] >= low_pt and arr1[i,1] >= low_pt and arr1[i,2] >= low_pt and arr1[i,3] >= low_pt:
            counter += 1
    return((counter/len(arr1[:,0])))

def inv_mass_recon_list(par_vec, vec1, vec2, vec3, vec4):
    '''Put in 4 vector classes and recombines to find the invariant mass of those with daughts in hte LHCb acceptance. 
    '''
    inv_mass_list = []
    for i in range(len(vec1)):
        if par_vec[i].eta >=1.596 and par_vec[i].eta <=5.298 and vec1[i].eta >= 1.596 and vec1[i].eta <= 5.298 and vec2[i].eta >=1.596 and vec2[i].eta <= 5.298 and vec3[i].eta >=1.596 and vec3[i].eta <=5.298 and vec4[i].eta >=1.596 and vec4[i].eta <= 5.298:
            E = vec1[i].e + vec2[i].e + vec3[i].e + vec4[i].e
            px = vec1[i].px + vec2[i].px + vec3[i].px + vec4[i].px
            py = vec1[i].py + vec2[i].py + vec3[i].py + vec4[i].py 
            pz = vec1[i].pz + vec2[i].pz + vec3[i].pz + vec4[i].pz
            P2 = px**2 + py**2 + pz**2
            E2 = E**2
            inv_mass = math.sqrt(abs(E2 - P2))
            inv_mass_list.append(inv_mass)
    return(inv_mass_list)


def frac_mass_recon_check(arr1, arr2, arr3, arr4, low_mass, high_mass):
    ''' Take in 4 arrays and will reconstruct their masses and returns fraction of
    how many values are with isn the mass rnage from the original aray length.
    Mass (GeV).
    '''
    counter =0
    list1 = inv_mass_recon_list(arr1, arr2, arr3, arr4)
    for i in range(len(list1)):
        if list1[i] >= low_mass and list1[i] <= high_mass:
            counter +=1
    return(counter/len(list1))

def hist(list1, parent, savefig, title, xlabel, ylabel, color, bins, label):
    arr1=np.array(list1)
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    plt.hist(arr1, color=color, bins=bins, histtype = 'step')
    axs.xaxis.set_tick_params(pad = 20) 
    axs.yaxis.set_tick_params(pad = 20)
    plt.axvline(arr1.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(arr1.mean()*1.1, max_ylim*0.9, 'Mean: {:.3f}'.format(arr1.mean()))
    plt.yscale('log')
    plt.legend([label],loc='upper center')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('Histograms/' + parent +'/' + savefig)

def hist_weight(list1, parent, savefig, title, xlabel, ylabel, color, bins, label, weights1):
    weights2 = np.full(len(list1), weights1)
    arr1=np.array(list1)
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    counts,bins, _ = plt.hist(arr1, color=color, bins=bins, histtype = 'step', weights=weights2)   
    bin_widths = np.diff(bins)
    integral = np.sum(counts * bin_widths)

    axs.xaxis.set_tick_params(pad = 20) 
    axs.yaxis.set_tick_params(pad = 20)
    # Plot the mean
    plt.axvline(arr1.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(arr1.mean()*1.1, max_ylim*0.9, 'Mean: {:.3f}'.format(arr1.mean()))
    #plt.yscale('log')

    plt.legend([label],loc='upper center')
    plt.title(title + '\nHistogram with Integral = {:.2f}'.format(integral))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('Histograms/' + parent +'/' + savefig)

def hist_weight_mass_curve(list1, parent, savefig, title1, color1, bins1, range1, weights1):
    weights2 = np.full(len(list1), weights1)
    arr1 = np.array(list1)
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7))#, 
                        #tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    counts, bins, _ = plt.hist(arr1, color=color1, bins=bins1, alpha=0.65, range =range1, histtype = 'step', weights = weights2)
    bin_widths = np.diff(bins)
    integral = np.sum(counts * bin_widths)

    axs.xaxis.set_tick_params(pad = 5) 
    axs.yaxis.set_tick_params(pad = 5)

    # plt.axvline(arr1.mean(), color='k', linestyle='dashed', linewidth=1)
    # min_ylim, max_ylim = plt.ylim()
    # plt.text(arr1.mean()*1.0001, max_ylim*0.9, 'Mass Mean: {:.3f}'.format(arr1.mean()))
    # plt.yscale('log')
    axs.legend(['Run 3'], loc= 'center right') #title="Hello"
    plt.title(title1 + '\nHistogram with Integral = {:.2f}'.format(integral))
    plt.xlabel('m($\mu^+$$\mu^-$$\mu^+$$\mu^-$) [MeV/C^2]')
    plt.ylabel('Number of Candidates')
    plt.savefig('Histograms/'+ parent +'/' + savefig)

def hist_mass_curve(list1, parent, savefig, title1, color1, bins1):
    arr1 = np.array(list1)
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7))#, 
                        #tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    plt.hist(arr1, color=color1, bins=bins1, alpha=0.65)
    axs.xaxis.set_tick_params(pad = 5) 
    axs.yaxis.set_tick_params(pad = 5)
    plt.axvline(arr1.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(arr1.mean()*1.0001, max_ylim*0.9, 'Mass Mean: {:.3f}'.format(arr1.mean()))
    plt.yscale('log')
    plt.title(title1)
    plt.xlabel('m($\mu^+$$\mu^-$$\mu^+$$\mu^-$) [GeV/C^2]')
    plt.ylabel('Number of Candidates')
    plt.savefig('Histograms/'+ parent +'/' + savefig)



def hist_pT_4mu(parent, savefig, title1, arr1):
    def gaussian(x, mean, amplitude, standard_deviation):
        return amplitude * np.exp( - (x - mean)**2 / (2*standard_deviation ** 2))

    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    axs.xaxis.set_tick_params(pad = 20) 
    axs.yaxis.set_tick_params(pad = 20)

    bin_heights1, bin_borders1, _ = plt.hist(arr1[:,3], bins='auto', label='Muon1', alpha=0.8, color = 'cornflowerblue', histtype = 'step')
    bin_centers1 = bin_borders1[:-1] + np.diff(bin_borders1) / 2
    # popt1, _ = curve_fit(gaussian, bin_centers1, bin_heights1, p0=[1., 0., 1.])
    plt.axvline(arr1[:,3].mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(arr1[:,3].mean()*1.1, max_ylim*0.95, '$\mu_1$ Mean: {:.2f}'.format(arr1[:,3].mean()))

    bin_heights2, bin_borders2, _ = plt.hist(arr1[:,2], bins='auto', label='Muon2', alpha=0.8, color = 'forestgreen', histtype = 'step')
    bin_centers2 = bin_borders2[:-1] + np.diff(bin_borders2) / 2
    # popt2, _ = curve_fit(gaussian, bin_centers2, bin_heights2, p0=[1., 0., 1.])
    plt.axvline(arr1[:,2].mean(), color='k', linestyle='dashed', linewidth=1)
    plt.text(arr1[:,2].mean()*1.1, max_ylim*0.80, '$\mu_2$ Mean: {:.2f}'.format(arr1[:,2].mean()))


    bin_heights3, bin_borders3, _ = plt.hist(arr1[:,1], bins='auto', label='Muon3', alpha=0.8, color = 'goldenrod', histtype = 'step')
    bin_centers3 = bin_borders3[:-1] + np.diff(bin_borders3) / 2
    # popt3, _ = curve_fit(gaussian, bin_centers3, bin_heights3, p0=[1., 0., 1.])
    plt.axvline(arr1[:,1].mean(), color='k', linestyle='dashed', linewidth=1)
    plt.text(arr1[:,1].mean()*1.1, max_ylim*0.65, '$\mu_3$ Mean: {:.2f}'.format(arr1[:,1].mean()))


    bin_heights4, bin_borders4, _ = plt.hist(arr1[:,0], bins='auto', label='Muon4', alpha=0.8, color = 'firebrick', histtype = 'step')
    bin_centers4 = bin_borders4[:-1] + np.diff(bin_borders4) / 2
    # popt4, _ = curve_fit(gaussian, bin_centers4, bin_heights4, p0=[1., 0., 1.])
    plt.axvline(arr1[:,0].mean(), color='k', linestyle='dashed', linewidth=1)
    plt.text(arr1[:,0].mean()*1.1, max_ylim*0.5, '$\mu_4$ Mean: {:.2f}'.format(arr1[:,0].mean()))

    # muon_fit1 = np.linspace(bin_borders1[0], bin_borders1[-1], 1000)
    # muon_fit2 = np.linspace(bin_borders2[0], bin_borders1[-1], 1000)
    # muon_fit3 = np.linspace(bin_borders3[0], bin_borders1[-1], 1000)
    # muon_fit4 = np.linspace(bin_borders4[0], bin_borders1[-1], 1000)
    #plt.yscale('log')
    plt.title(title1)
    plt.xlabel('pT (GeV)')
    plt.ylabel('Number of Candidates Per pT')
    # plt.plot(muon_fit1, gaussian(muon_fit1, *popt1), label='Muon1 fit', color='royalblue')
    # plt.plot(muon_fit2, gaussian(muon_fit2, *popt2), label='Muon2 fit', color='g')
    # plt.plot(muon_fit3, gaussian(muon_fit3, *popt3), label='Muon3 fit', color='darkorange')
    # plt.plot(muon_fit4, gaussian(muon_fit4, *popt4), label='Muon4 fit', color='brown')
    plt.legend()
    plt.savefig('Histograms/'+ parent + '/'+ savefig)
