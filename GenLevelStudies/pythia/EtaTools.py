import numpy as np
import matplotlib.pyplot as plt
import uproot
import vector


def LApt(a, b, c, d):
    import vector
    LongList =  [a.pt, b.pt, c.pt, d.pt]
    for i in range(len(LongList)):
        if i == 0:
            LongArray = np.array(LongList[0])
        else:
            LongArray = np.append(LongArray, np.array(LongList[i]))
    return(LongArray)

def LAeta(a, b, c, d):
    import vector
    LongList =  [a.eta, b.eta, c.eta, d.eta]
    for i in range(len(LongList)):
        if i == 0:
            LongArray = np.array(LongList[0])
        else:
            LongArray = np.append(LongArray, np.array(LongList[i]))
    return(LongArray)

def TwoDEta(a, b, c, d):
    TwoDimArray = np.transpose([np.array(a.eta),np.array(b.eta),np.array(c.eta),np.array(d.eta)])
    return(TwoDimArray)