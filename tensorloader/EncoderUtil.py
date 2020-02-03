import numpy as np
import pandas as pd

#Author: Asa Thibodeau

def getFrequencyVector(line):
    """Returns a numpy array of values from a comma
        delimited string.
        
        Parameters
        ----------
        line : str
        A comma delimited string encoding numeric
        values
        
        Returns
        -------
        rv : (#length) numpy array
        A number array of numeric values equal to the
        number of values delimited by commas.
        
        """
    line = line.strip()
    return np.array(line.split(","), dtype=float)

def getPeakPositions(summitpeaks, peaks, l):
    """Gets a binary vector corresponding to the
        peak position within the specified window
        of length l.
        
        Parameters
        ----------
        summitpeaks : (#numpeaks, 3) numpy array
        The extended peak locations
        chr start end
        
        peaks : (#numpeaks, 3) numpy array
        The original peak locations.
        chr start end
        
        l : int
        The total length of each region.
        
        Returns
        -------
        rv : (#length) numpy array
        A binary array corresponding to whether or
        not the original peak overlaps the location.
        """
    #Get the original peak positions relative to the padded peak size 
    positions = np.array(peaks[:,1:3], dtype=int)
    positions[:,0] = positions[:,0]-summitpeaks[:,1]
    positions[:,1] = positions[:,1]-summitpeaks[:,1]+1
    positions[positions[:,:] >= l] = l
    positions[positions[:,:] < 0] = 0
    #Initialize all positions to 0 (i.e., No peak)
    rv = np.zeros((len(summitpeaks), l))
    #assign original peak positions 1, 
    for i in range(len(peaks)):
        rv[i, positions[i,0]:positions[i,1]] = 1
    return rv