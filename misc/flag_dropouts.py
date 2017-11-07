import numpy as np
import matplotlib.pyplot as plt
#import statistics as stats
import spacepy.time as spt
import spacepy.datamodel as dm
import datetime
import matplotlib.gridspec as gridspec

def offsetDeriv(counts, DERIV_OFFSET):
    """
    This function takes a derivative of counts array with a DERIV_OFFSET
    derivative offset, ie
    dc/dt[i] = (counts[i] - counts[i+DERIV_OFFSET])/DERIV_OFFSET
    """
    dcdt = np.zeros_like(counts)
    kernel = np.zeros(DERIV_OFFSET + 1)
    kernel[0] = 1
    kernel[-1] = -1
    if DERIV_OFFSET == 1:
        lbound = 0
        ubound = -1
    else:
        lbound = int(DERIV_OFFSET/2)
        ubound = -int(DERIV_OFFSET/2)
    
    dcdt[lbound:ubound] = np.convolve(kernel, counts, mode= 'valid')/DERIV_OFFSET
    
    return dcdt
    
def dropOutFlag(counts, DERIV_OFFSET = 1, DERIV_THRESH = 200, FLAG_WIDTH = 5):
    """
    This function identifies dropouts and flags it with a 1. 0 means it is good
    data.
    """
    dcdt = offsetDeriv(counts, DERIV_OFFSET)
    dropoutFlag = np.zeros_like(counts)
    
    # Set the dropout flag.             
    for i in range(len(counts)-DERIV_OFFSET):
        if dcdt[i-DERIV_OFFSET] <= -DERIV_THRESH or \
                dcdt[i+DERIV_OFFSET] >= DERIV_THRESH:
            # Now check if at least two channels drop to the same count value.
            #if len(set(hr['hr0'][i, :])) < 6:
            dropoutFlag[i-FLAG_WIDTH:i+FLAG_WIDTH] = 1
    
    return dropoutFlag
    
def dropOutInterp(counts, dropoutFlag):
    """
    This function does a linear interpolation over the data where the dropout
    flag is set to 1.
    """
    dataInterpol = np.array(counts, copy = True)
    starti = 0
    endi = starti
    while starti < len(dropoutFlag):
        if dropoutFlag[starti] == 1:
            endi = starti
            while dropoutFlag[endi] == 1:
                endi += 1
            else:
                # Now Interpolate over the dropout
                dataInterpol[starti:endi] = \
                        np.interp(ts[starti:endi], \
                        xp = [ts[starti], ts[endi]], \
                        fp = [CH1[starti], CH1[endi]])
                starti = endi
        else:
            starti += 1
    
if __name__ == '__main__':
    #Test harness
    if 'hr' not in globals():
        folder = "/home/mike/FIREBIRD/Datafiles/FB_FD/FU_3/Level1/"
        fname = 'FU_3_Hires_2015-02-02_L1_v02.txt'
        hr = dm.readJSONheadedASCII(folder + fname)
    
        CH1 = hr["hr0"][:, 0]
        td = hr["Epoch"]
        tp = spt.Ticktock(td)
        time = tp.UTC
        
        ts = np.zeros(len(time))
        for i in range(len(time)):
               ts[i] =  (time[i] - datetime.datetime(2015, 1, 1, 00, 00, 00, 0000)).total_seconds()    
    
    dropoutFlag = dropOutFlag(CH1)

    gs = gridspec.GridSpec(3, 3)
    data = plt.subplot(gs[0:2, :])
    derivative = plt.subplot(gs[2, :], sharex = data)
    
    derivative.plot(time, 200*dropoutFlag, label = 'Dropout flag', linewidth = 3)
    #derivative.plot(time, dcdt, label = 'dc/dt', linewidth = 3)
    derivative.set_ylim([-200,200])
    derivative.legend()
    for i in hr['Channel']:
        data.plot(time, hr['hr0'][:, int(i)], label = 'CH' + str(int(i)))
    #data.plot(time, dataInterpol, label = 'Interpolated CH 1')
    data.legend()    
    
    plt.show()
    plt.tight_layout()
