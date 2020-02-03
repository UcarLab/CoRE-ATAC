import numpy as np
import pandas as pd
from tensorloader import EncoderUtil as u

#Author: Asa Thibodeau
class BaseToInt:
    """Initializes a dictionary which maps characters
        that correspond to one or more bases.
        """
    def __init__(self):
        self.basetoint = dict()
        self.basetoint["A"] = (0,)
        self.basetoint["C"] = (1,)
        self.basetoint["M"] = (0,1,)
        self.basetoint["G"] = (2,)
        self.basetoint["R"] = (0,2,)
        self.basetoint["S"] = (1,2,)
        self.basetoint["V"] = (0,1,2,)
        self.basetoint["T"] = (3,)
        self.basetoint["W"] = (0,3,)
        self.basetoint["Y"] = (1,3,)
        self.basetoint["H"] = (0,1,3,)
        self.basetoint["K"] = (2,3,)
        self.basetoint["D"] = (0,2,3,)
        self.basetoint["B"] = (1,2,3,)
        self.basetoint["N"] = (0,1,2,3,)
    
    def getPositions(self, base):
        """List the bases in integers that correspond
            to a specific character.

            0 - A
            1 - C
            2 - G
            3 - T

            Parameters
            ----------
            base : str
            The character to convert to base pair integers.

            Returns
            -------
            rv : tuple
            A tuple of all bases (in integers) that correspond
            to the character provided.
            """
        return self.basetoint[base]

def updateReferenceAndNormalize(m, ref, thresh):
    """Updates a sequence matrix, filling in values
        from the reference when the number of reads for
        a position is less than the threshold provided.

        Parameters
        ----------
        m : (5, #len) numpy array
        The matrix containing the raw counts that will
        be normalized to [0,1] by columnwise minmax
        normalization and updated with the reference
        sequence provided.
        """
    ref = list(ref)
    thresh = max(thresh,0)
    totals = np.sum(m[:4,], axis=0)
    idx = 0;
    b2i = BaseToInt()
    for i in totals:
        if i < thresh:
            bases = np.array(b2i.getPositions(ref[idx].capitalize()))
            m[:4, idx] = 0
            m[bases, idx] = 1.0/len(bases)
        else:
            m[:4,idx] = m[:4,idx]/i

        #DEBUG CODE#
        if (m[:4,idx] > 1).any():
            print(i)
            print (m[:4,idx])
            print(totals)
        #END DEBUG CODE#
        
        idx += 1
        
        
def getSequenceFeatures(seqfile, summitpeaks, peaks, spidx, l, rthresh):
    """Reads dataset information corresponding to
        dataset file location and dataset labels.
        
        Parameters
        ----------
        seqfile : str
        Path of the file containing sequence features.
        
        summitpeaks : (#numpeaks, 3) numpy array
        The extended peak locations
        chr start end
        
        peaks : (#numpeaks, 3) numpy array
        The original peak locations.
        chr start end
        
        spidx : (#numpeaks) numpy array
        An sorted index of peak ids for properly matching
        data between files.
        
        l : int
        The total length of each region.
        
        rthresh : int
        The number of reads at a position required to
        override reference genome base.
        
        Returns
        -------
        rv : (#peaks, 4, #len) numpy array
        A (#peaks, 4, #len) tensor of sequence features.
        
        0:Peak position
        1:A
        2:C
        3:G
        4:T
        
        idx : (#peaks) numpy array
        A vector containing the ids of each region.
        """
    numpeaks = len(peaks)
    rv = np.zeros((numpeaks, 5, l),dtype=float)
    idx = np.zeros((numpeaks), dtype=object)
    curidx = 0
    with open(seqfile, 'r') as f:
        while True:
            try:
                metadata = next(f)
                peakid = int(metadata.split("\t")[3].strip())
                reference = next(f).strip()
                fv1 = u.getFrequencyVector(next(f))
                fv2 = u.getFrequencyVector(next(f))
                fv3 = u.getFrequencyVector(next(f))
                fv4 = u.getFrequencyVector(next(f))
                if peakid in spidx:
                    idx[curidx] = peakid
                    rv[curidx, 0, :] = fv1 #A
                    rv[curidx, 1, :] = fv2 #C
                    rv[curidx, 2, :] = fv3 #G
                    rv[curidx, 3, :] = fv4 #T
                    updateReferenceAndNormalize(rv[curidx], reference, rthresh)
                    curidx += 1
            except StopIteration as se:
                break

    sortedidx = np.argsort(idx)
    idx = idx[sortedidx]
    assert np.all(spidx == idx)
    rv = rv[sortedidx]
    
    rv[:, 4, :] = u.getPeakPositions(summitpeaks, peaks, l)
    
    return rv
