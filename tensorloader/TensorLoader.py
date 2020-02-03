import numpy as np
import pandas as pd
from tensorloader import PEASUtil
from sklearn import preprocessing
from sklearn import impute
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import PolynomialFeatures
import scipy
from scipy import signal
from tensorloader import SequenceEncoder as seqenc
from tensorloader import SignalEncoder as sigenc
import os



def getMetaData(peakfile):
    """Gets the extended peak, original peak,
        peak ids, and annotations.
        
        Parameters
        ----------
        peakfile : str
        Path of the meta data file.

        Returns
        -------
        locations : (#numpeaks) numpy array
        The extended locations sorted by id.
        
        peaks : (#numpeaks) numpy array
        The original peak locations sorted
        by id.
        
        spidx : (#numpeaks) numpy array
        The sorted peak id index.
        
        annotations : (#numpeaks) numpy array
        The corresponding annotations (zeros
        if not provided)
        
        """
    paddedpeaks = pd.read_csv(peakfile, sep="\t", header=None).values #TODO , header=None
    filter = (paddedpeaks[:,0] != "chrX") & (paddedpeaks[:,0] != "chrY") #TODO Allow for chrX and Y
    paddedpeaks = paddedpeaks[filter,:] #TODO Allow for chrX and Y
    locations = paddedpeaks[:,(0,1,2)]
    peaks = paddedpeaks[:,(0,4,5)]
    spidx = paddedpeaks[:, 3].astype(int)
    sortedidx = np.argsort(spidx)
    annotations = np.zeros((len(spidx)))
    if(np.shape(paddedpeaks)[1] > 6):
        annotations = paddedpeaks[:, 6]
    return locations[sortedidx], peaks[sortedidx], spidx[sortedidx], annotations[sortedidx]

def getSeqSigTensor(seq, sig):
    return np.concatenate((seq[:,:4,:], sig), axis=1)

def readTensors(name, directory, numbases, sequence=True, signal=True, rthresh=10, normalized=True):
    """Reads ATAC-seq tensors into python for deep learning.
    
        Parameters
        ----------
        name : str
        Prefix given to the data files.
        
        directory : str
        Directory where files are stored.
        
        numbases : int
        The total length of each region.
        Should be the same for each.
        
        sequence : bool (default=True)
        Whether to process sequence features.
        
        signal : bool (default=True)
        Whether to process signal features.
        
        cutmatrix : bonormalizeSignalTensorol (default=True)
        Whether to process cut matrices.
        
        rthresh : int (default=10)
        The minimum number of reads needed
        to override reference genome.
        
        normalized : bool (default=True)
        Sets whether or not to normalize the
        signal related features. Useful for
        normalization after partitioning.
        
        Returns
        -------
        rv1 : (#peaks, 4, #len) numpy array
        A (#peaks, 4, #len) tensor of sequence
        features.
        
        rv2 : (#peaks, 5, #len) numpy array 
        Tensor containing signal features of
        all peaks
        
        annotations : (#numpeaks) numpy array
        The corresponding annotations (zeros
        if not provided)
        
        paddedpeaks : (#numpeaks) numpy array
        The extended locations sorted by id.
        
        peaks : (#numpeaks) numpy array
        The original peak locations sorted
        by id.
        """
    rv1 = None
    rv2 = None
    rv3 = None
    rv4 = None
    annotations = None
    paddedpeaks = None
    peaks = None
    try:
        paddedpeaks, peaks, spidx, annotations = getMetaData(os.path.join(directory, name+"_peaks.txt")) #_summit_peaks_roadmap4.txt" #TODO
                
        if sequence:
            rv1 = seqenc.getSequenceFeatures(os.path.join(directory, name+"_sequencefreq.txt"), paddedpeaks, peaks, spidx, numbases, rthresh)
            print("Sequence data loaded.")
            
        if signal:
            rv2 = sigenc.getSignalFeatures(name, directory, spidx, paddedpeaks, peaks, numbases, normalized)
            print("Signal data loaded.")
            
    except AssertionError as ae:
        print("Peak Ids failed to match.")

    return rv1, rv2, annotations, paddedpeaks, peaks

def getPEASFeatures(filepath, featurefile, labelencoderfile, peaks):
    data = pd.read_csv(filepath, sep="\t")
    featurecolumns = PEASUtil.getFeatureColumnData(featurefile)
    labelencoder = PEASUtil.getLabelEncoder(labelencoderfile)
    chrpositions = data.values[:,:3]
    selecteddata = data.values[:, featurecolumns]

    onehotindices = []
    onehotvalues = []
    for i in range(len(labelencoder)):
        curonehotindex = np.argwhere(featurecolumns == labelencoder[i][0])[0][0]
        onehotindices.append(curonehotindex)
        curonehotvalues = np.array(labelencoder[i][1:])
        curonehotencoder = OneHotEncoder(categories=[curonehotvalues], sparse=False, handle_unknown='ignore')
        curvector = selecteddata[:,curonehotindex].reshape(-1,1)
        for i in range(len(curvector)):
            if str(curvector[i,0]) == 'nan':
                curvector[i,0] = ''
        onehotvalues.append(curonehotencoder.fit_transform(curvector))

    testX = selecteddata[:,~np.in1d(np.arange(np.shape(selecteddata)[1]), onehotindices)]

    imputer = impute.SimpleImputer(missing_values=np.nan, strategy='mean')
    sscaleddata = preprocessing.StandardScaler().fit_transform(imputer.fit_transform(testX))
    
    combinedvalues = [sscaleddata] + onehotvalues

    sdata = np.concatenate(combinedvalues, axis=1).astype(float)
    
    peakmap = dict()
    for i in range(len(sdata)):
        curchr = chrpositions[i,0]
        curstart = chrpositions[i,1]
        curend = chrpositions[i,2]
        peakmap[curchr+":"+str(curstart)+"-"+str(curend)] = sdata[i,:]
        
    #step 2: use dictionary to populate a tensor with peaks
    rv = []
    for i in range(len(peaks)):
        curkey = peaks[i][0]+":"+str(peaks[i][1])+"-"+str(peaks[i][2])
        if curkey in peakmap:
            rv.append(peakmap[curkey])
        else:
            print("Error: Couldn't find "+curkey+".")
    
    
    #return np.array(PolynomialFeatures(include_bias=False).fit_transform(rv))
    return np.array(rv)

