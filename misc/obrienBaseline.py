import numpy as np
import time
import multiprocessing

def obrienBaseline(data, timeWidth = 1.0, cadence = 18.75E-3, PERCENTILE_THRESH = 10):
    """
    NAME:    obrienBaseline(data, timeWidth = 1.0, cadence = 18.75E-3, PERCENTILE_THRESH = 10)
    USE:     Implements the baseline algorithm described by O'Brien et al. 2004, GRL. 
             The timeWidth parameter is the half of the range from which to calcualte the 
             baseline from, units are seconds. Along the same lines, cadence is used
             to determine the number of data points which will make the timeWidth
    RETURNS: A bseline array that represents the bottom 10th precentile of the data. 
    MOD:     2016-08-02
    """
    start_time = time.time()
    dataWidth = int(np.floor(timeWidth/cadence))

    assert 2*dataWidth < len(data), 'Data too short for timeWidth specified.'
    #if 2*dataWidth > len(data):
    #   raise IndexError('Data too short for timeWidth specified.')
    baseline= np.zeros_like(data, dtype = float)
    
    for i in range(int(len(baseline)-dataWidth)):
        i = int(i)
        dataSlice = data[i:2*dataWidth+i]
        baseline[dataWidth+i] = np.percentile(dataSlice, PERCENTILE_THRESH)
    print("O'Brien Baseline time: {}".format(time.time() - start_time))
    return baseline
    
def obrienBaselineMultiprocess(data, timeWidth, cadence, PERCENTILE_THRESH, channel, out_q):
    """
    NAME:    obrienBaseline_multiprocess(data, timeWidth = 1.0)
    USE:     Implements the baseline algorithm described by O'Brien et al. 2004, GRL. 
             The timeWidth parameter is the half of the range from which to calcualte the 
             baseline from, units are seconds. Along the same lines, cadence is used
             to determine the number of data points which will make the timeWidth.
             Note: This is for the multiprocessing module!
    RETURNS: A bseline array that represents the bottom 10th precentile of the data. 
    MOD:     2016-08-02
    """
    dataWidth = int(np.floor(timeWidth/cadence))
    if 2*dataWidth > len(data):
        raise IndexError('Data too short for timeWidth specified.')
    baseline= np.zeros_like(data, dtype = float)
    
    for i in range(int(len(baseline)-dataWidth)):
        i = int(i)
        dataSlice = data[i:2*dataWidth+i]
        baseline[dataWidth+i] = np.percentile(dataSlice, PERCENTILE_THRESH)
    out_q.put({channel:baseline})
    
def baselineMultiprocessWrapper(counts, baslineTimeWidth, PERCENTILE_THRESH, cadence = 18.75E-3):
    """
    NAME: firebirdBaselineMultiprocessWrapper(counts, baslineTimeWidth, PERCENTILE_THRESH, cadence = 18.75E-3)
    USE:   Calculates Paul's baseline using Python's multiprocessing module. 
    RETURNS: Baseline dictionay with channel number keys, with arrays of len nT. (This could be fixed up a little)
    MOD:     2017-10-25 
    """
    start_time = time.time()
    numProcs = counts.shape[1]
    out_q = multiprocessing.Queue()
    num_cores = multiprocessing.cpu_count()
    procs = []
    baseline = {}
    
    batch = 0
    
    # The logic behind dividing up the HiRes data between CPU cores
    if num_cores >= numProcs:
        print("Calculating baseline on a machine with more cores than processes.") 
        # Do in one swoop
        while batch < numProcs:
            p = multiprocessing.Process(target=obrienBaselineMultiprocess, \
            args=(counts[:, batch], baslineTimeWidth, cadence, \
            PERCENTILE_THRESH, batch, out_q))
                
            procs.append(p)
            p.start()
            batch += 1
             
        for i in range(numProcs):
            baseline.update(out_q.get())    
        
        # Wait for all worker processes to finish
        for p in procs:
            p.join()
    else:
        # Do in chunks
        print("Calculating baseline on a machine with less cores than processes.")
        while batch < numProcs:
            # Use all but one core.
            if batch + num_cores - 1 <= numProcs:
                lbound = batch
                ubound = batch + num_cores - 1
            else:
                lbound = batch
                ubound = numProcs

            for i in range(lbound, ubound):
                p = multiprocessing.Process(target=obrienBaselineMultiprocess, \
                args=(counts[:, batch], baslineTimeWidth, cadence, \
                PERCENTILE_THRESH, batch, out_q))
                    
                procs.append(p)
                p.start()
                batch += 1
                
            for i in range(lbound, ubound):
                baseline.update(out_q.get())                     
            batch = ubound

    print("O'Brien Multiprocess Baseline time: {}".format(time.time() - start_time))
    return baseline
   
## TEST SCRIPT FOR MULTIPROCESS STARTS HERE ##
#if __name__ == '__main__': 
#    import spacepy.datamodel as dm
#    folder = "/home/mike/FIREBIRD/Datafiles/FB_FD/FU_3/Level1/"
#    fname = "FU_3_Hires_2015-02-02_L1_v02.txt"
#    hr = dm.readJSONheadedASCII(folder + fname)
#     
#    baslineTimeWidth = 1.0 # s
#    PERCENTILE_THRESH = 20
#    cadence = 18.75E-3 #s
#    
#    baseline = firebirdBaselineMultiprocessWrapper(hr, baslineTimeWidth, PERCENTILE_THRESH, cadence)
