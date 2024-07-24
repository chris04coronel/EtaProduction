import numpy as np
import matplotlib.pyplot as plt
import uproot
import vector
import time


def LApt(vec1, vec2, vec3, vec4):
    '''Creates a 1D array made from 4 independent arrays
    
    Take 4 four vectors find their pT 
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
    print("My eta acc function took", time.time() - start_time, "to run")
    return(vec_list_eta_acc_pt)

def sort_pt(vec1, vec2, vec3, vec4):
    'Take 4 four vectors and sorts them in a nx4. Where n is the length of the vectors (Assuming)'
    'they are all the same length. the the 1st entry is the smallest pT in ascending order'
    TwoDimArray = np.transpose([np.array(vec1.pt),np.array(vec2.pt),np.array(vec3.pt),np.array(vec4.pt)])
    for i in range(len(TwoDimArray[:,1])):
        if i == 0:
            sort_ndarray = np.sort(TwoDimArray[i,:])
        else:
            sort_ndarray = np.vstack((sort_ndarray, np.sort(TwoDimArray[i,:])))
    return(sort_ndarray)

def eta_acc_sort_pt(vec1, vec2, vec3, vec4):
    'looks for all 4 muons in acceptance and then sorts them by ascending pT returns nx4 array'
    start_time=time.time()
    counter = 0
    for i in range(len(vec1)):
        if vec1[i].eta >=2 and vec1[i].eta <= 5 and vec2[i].eta >=2 and vec2[i].eta <= 5 and vec3[i].eta >=2 and vec3[i].eta <=5 and vec4[i].eta >=2 and vec4[i].eta <= 5:
            if counter == 0:
                counter += 1
                arr1 =np.array([vec1[i].pt, vec2[i].pt, vec3[i].pt, vec3[i].pt])
                sort_arr1 = np.sort(arr1)
            else:
                arr1 =np.array([vec1[i].pt, vec2[i].pt, vec3[i].pt, vec3[i].pt])
                sort_arr1 = np.vstack((sort_arr1, np.sort(arr1)))
    print("My eta acc pT sorting function took", time.time() - start_time, "to run")
    return(sort_arr1)    

def bar_hist(list1, savefig, title, xlabel, ylabel, color, bins):
    fig, axs = plt.subplots(1, 1,
                        figsize =(10, 7), 
                        tight_layout = True)
    axs.grid(visible = True, color ='grey', 
        linestyle ='-.', linewidth = 0.5, 
        alpha = 0.6)
    counts, edges, bars = plt.hist(list1 ,edgecolor='black', color=color, bins=bins)
    plt.bar_label(bars)
    axs.xaxis.set_tick_params(pad = 20) 
    axs.yaxis.set_tick_params(pad = 20)
    array1=np.array(list1)
    plt.axvline(array1.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(array1.mean()*1.1, max_ylim*0.9, 'Mean: {:.2f}'.format(array1.mean()))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig('pythia/Histograms/'+savefig)