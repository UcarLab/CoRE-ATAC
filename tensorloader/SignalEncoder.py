import numpy as np
import pandas as pd
from tensorloader import EncoderUtil as u
import scipy
import os

#Author: Asa Thibodeau

def normalizeSignalTensor(m):
    """Normalizes a signal tensor after loading.
        This is useful for partitioning data into
        training/test sets without leaking test
        data into the training set through 
        normalization.
        
        Parameters
        ----------
        m : (#peaks, 5, #len) numpy array 
        Non-normalized tensor of signal features.
        
        Returns
        ----------
        rv : (#peaks, 5, #len) numpy array 
        Normalized tensor of signal features.
        """
    rv = np.zeros(np.shape(m))
    rv[:,0,:] = allBasesStandardize(m[:,0,:])
    rv[:,1,:] = allBasesStandardize(m[:,1,:])
    rv[:,2,:] = allBasesStandardize(m[:,2,:])
    rv[:,3,:] = minMaxMedianInserts(m[:,3,:])
    rv[:,4,:] = minMaxMedianInserts(m[:,4,:])
    return rv


def getSignalData(infile, spidx, l):
    """Obtains a single vector of
        signal features.
    
        Parameters
        ----------
        infile : str
        Path of the signal data file.

        spidx : (#numpeaks) numpy array
        The sorted peak id index.
        
        l : int
        The total length of each region.
        
        Returns
        -------
        rv : tensor
        A (#peaks, 1, #len) tensor of signal
        features.
        
        idx : (#peaks) numpy array
        A vector containing the ids of each region.
        """    
    numpeaks = len(spidx)
    rv = np.zeros((numpeaks, 1, l),dtype=float)
    idx = np.zeros((numpeaks), dtype=object)
    curidx = 0
    with open(infile, 'r') as f:
        while True:
            try:
                peakinfo = next(f)
                fv = u.getFrequencyVector(next(f).strip())
                peakid = int(peakinfo.split("\t")[3].strip())
                if peakid in spidx:
                    idx[curidx] = peakid
                    rv[curidx, 0] = fv
                    curidx += 1
            except StopIteration as se:
                break
    sortedidx = np.argsort(idx)
    idx = idx[sortedidx]
    rv = rv[sortedidx]
    assert np.all(spidx == idx)
    
    return rv

def minMaxMedianInserts(m):
    """Minmax normalization method for
        median insert sizes.
    
        Parameters
        ----------
        m : (#peaks, #len, 1) numpy array
        Vector of median inserts.

        Returns
        -------
        m' : (#peaks, 1, #len) numpy array 
        Minmax normalized vector of median
        inserts.
        """
    minmedian = np.min(m)
    maxmedian = np.max(m)
    return (m-minmedian)/maxmedian

def normNMedianInserts(m,n):
    """N normalization method for
        median insert sizes.

        Values are divided by a n
        and trimmed to a max value of
        1 if values are larger than n.

        Parameters
        ----------
        m : (#peaks, #len, 1) numpy array
        Vector of median inserts.

        n : int
        Constant used to divide the values.

        Returns
        -------
        m' : (#peaks, 1, #len) numpy array 
        N normalized vector of median
        inserts.
        """
    rv = m/n
    rv[rv[:,:,:] > 1.0] = 1.0
    return rv;

def allBasesStandardize(m):
    """Standardizes signal related features
        by counting all bases.
        
        Parameters
        ----------
        m : (#peaks, #len, 1) numpy array
        Vector of pileup information such
        as inserts or cuts.
        
        Returns
        -------
        m' : (#peaks, 1, #len) 
        Standardized vector of values.
        """
    mu = np.mean(m)
    dev = np.std(m)
    return (m-mu)/dev


def standardizeAndSmooth(m):
    """Smooths the cut sites and standardizes by total number of cuts observed in the region.
        
        Parameters
        ----------
        m : (#peaks, #len, 1) numpy array
        Vector of pileup information such
        as inserts or cuts.
        
        smooth : smoothing parameter for scipy.interpolate.Rbf
        
        Returns
        -------
        m' : (#peaks, 1, #len) 
        Standardized vector of values.
        """
    mu = np.mean(m)
    dev = np.std(m)
    m2 = (m-mu)/dev
    
    mshape = np.shape(m)
    rv = np.zeros(mshape)
    
    for i in range(mshape[0]):
        rv[i,0,:] = scipy.signal.savgol_filter(m2[i,0,:], 15, 2)
    
    return rv

def allPeaksStandardize(m, summitl, peaks, mlen):
    """Standardizes signal related features
        by counting only the bases within
        the original peak calls.
        
        This is implemented just in case
        the padded zeros cause issues
        when adjusting for read depth.
        
        Parameters
        ----------
        m : (#peaks, #len, 1) numpy array
        Vector of pileup information such
        as inserts or cuts.
        
        summitl : (#numpeaks) numpy array
        The extended locations.
        
        peaks : (#numpeaks) numpy array
        The original peak locations.
        
        mlen : int
        The total length of each region.
        
        Returns
        -------
        m' : (#peaks, 1, #len) 
        Standardized vector of values based
        on original peaks.
        """
    selectedvalues = []
    for i in range(len(midx)):
        lstart = summitl[i][1]
        peakstart = peaks[i][1]
        peakend = peaks[i][2]
        sidx = max(peakstart-lstart,0)
        eidx = min(peakend-lstart, maxlen)
        selectedvalues.extend(list(m[i,0,sidx:eidx]))
    
    selectedarray = np.array(selectedvalues)
    mu = np.mean(selectedvalues)
    dev = np.std(selectedarray)
    return (m-mu)/dev

def getCompiledSignalFeatures(inserts, fcuts, rcuts, fmedins, rmedins, summitpeaks, peaks, mlen):
    """Compiles signal features into a
        single tensor.
        
        Parameters
        ----------
        insert : (#peaks, 1, #len) numpy array
        Vector of insert pileup counts.
        
        fcuts : (#peaks, 1, #len) numpy array
        Vector of forward cut pileup counts.
        
        fcuts : (#peaks, 1, #len) numpy array
        Vector of reverse cut pileup counts.
        
        fmedins : (#peaks, 1, #len) numpy array
        Vector of forward median inserts.
        
        rmedins : (#peaks, 1, #len) numpy array
        Vector of reverse median inserts.
        
        mlen : int
        The total length of each region.
        
        Returns
        -------
        rv : (#peaks, 5, #len) numpy array 
        Tensor combining all vectors.
        """
    numpeaks = len(inserts)
    rv = np.zeros((numpeaks, 6, mlen), dtype=float)
    for i in range(numpeaks):
        rv[i,0] = inserts[i,0]
        rv[i,1] = fcuts[i,0]
        rv[i,2] = rcuts[i,0]
        rv[i,3] = fmedins[i,0]
        rv[i,4] = rmedins[i,0]
    
    rv[:, 5, :] = u.getPeakPositions(summitpeaks, peaks, mlen)
    
    return rv

def getSignalShape(insertpileups, mlen, numrows):
    rv = []
    for curpileup in insertpileups:
        curpileup = curpileup[0]
        curmax = np.max(curpileup)
        curmax = max(curmax,1)
        percentage = (curpileup*numrows/curmax)
        currv = np.zeros((numrows, mlen))
        for curcol in range(mlen):
            height = int(percentage[curcol]+1)
            currv[-height:,curcol] = 1
        rv.append(currv)
    return np.array(rv)

def getSignalFeatures(name, directory, spidx, summitpeaks, peaks, mlen, normalized=True):
    """Reads signal features into a tensor.
        
        Parameters
        ----------
        name : str
        Prefix provided to the feature files,
        labeling them.
        
        directory : str
        Directory containing the files.
        
        spidx : (#numpeaks) numpy array
        The sorted peak id index.
        
        l : int
        The total length of each region.
        
        normalized : bool (default=True)
        Sets whether or not to normalize the
        signal related features. Useful for
        normalization after partitioning.
        
        Returns
        -------
        rv : (#peaks, 5, #len) numpy array 
        Tensor containing signal features of
        all peaks.
        
        idx : (#peaks) numpy array
        A vector containing the ids of each region.
        """
    insertfile = os.path.join(directory, name+"_insertpileups.txt")
    fcutfile = os.path.join(directory, name+"_forwardcutpileups.txt")
    rcutfile = os.path.join(directory, name+"_reversecutpileups.txt")
    fmedinsfile = os.path.join(directory, name+"_forwardmedianinsert.txt")
    rmedinsfile = os.path.join(directory, name+"_reversemedianinsert.txt")
    
    inserts = getSignalData(insertfile, spidx, mlen)
    fcuts = getSignalData(fcutfile, spidx, mlen)
    rcuts = getSignalData(rcutfile, spidx, mlen)
    fmedins = getSignalData(fmedinsfile, spidx, mlen)
    rmedins = getSignalData(rmedinsfile, spidx, mlen)
    
    if normalized:
        inserts = allBasesStandardize(inserts)
        fcuts = standardizeAndSmooth(fcuts)#allBasesStandardize(fcuts)
        rcuts = standardizeAndSmooth(rcuts)#allBasesStandardize(rcuts)
        fmedins = standardizeAndSmooth(fmedins) #minMaxMedianInserts(fmedins)
        rmedins = standardizeAndSmooth(rmedins) #minMaxMedianInserts(rmedins)

    return getCompiledSignalFeatures(inserts, fcuts, rcuts, fmedins, rmedins, summitpeaks, peaks, mlen)
