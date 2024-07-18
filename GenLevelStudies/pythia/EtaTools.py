import numpy as np
import matplotlib.pyplot as plt
import uproot
import vector


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

def eta_acc_LApt(vec1, vec2, vec3, vec4):
    'takes four 4vectors cuts off those which dont meet eta(pseudorapidity) requierments'
    import vector
    import time
    start_time=time.time()
    veclist = [vec1, vec2, vec3, vec4]
    vec1_acc_pt = []
    vec2_acc_pt = [] 
    vec3_acc_pt = [] 
    vec4_acc_pt = []
    #veclist_eta_acc = [vec1_acc, vec2_acc, vec3_acc, vec4_acc]

    for i in range(len(vec1)):
        if vec1[i].eta >=2 and vec1[i].eta <= 5:
            vec1_acc_pt.append(vec1[i].pt)
    for i in range(len(vec2)):
        if vec1[i].eta >=2 and vec1[i].eta <= 5:
            vec2_acc_pt.append(vec2[i].pt)
    for i in range(len(vec3)):
        if vec3[i].eta >=2 and vec3[i].eta <= 5:
            vec3_acc_pt.append(vec3[i].pt)
    for i in range(len(vec4)):
        if vec4[i].eta >=2 and vec4[i].eta <= 5:
            vec4_acc_pt.append(vec4[i].pt)
    vec_list_eta_acc_pt = vec1_acc_pt + vec2_acc_pt + vec3_acc_pt + vec4_acc_pt
    print("My program took", time.time() - start_time, "to run")
    return(vec_list_eta_acc_pt)
    