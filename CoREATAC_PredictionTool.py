from __future__ import print_function
import tensorflow as tf
import keras
from tensorflow.keras.models import load_model
from keras import backend as K
from keras.layers import Input
import numpy as np
import subprocess
from tensorloader import TensorLoader as tl
import partitioning_util as part
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_recall_curve,average_precision_score, confusion_matrix
import pandas as pd
from sklearn import impute
import argparse
import os
import time


#Step 0: Process arguments
parser = argparse.ArgumentParser(description='CoRE-ATAC Prediction Tool')
parser.add_argument("datadirectory")
parser.add_argument("model")
parser.add_argument("outputfile")
args = parser.parse_args()

datadirectory = args.datadirectory
model = args.model
outputfile = args.outputfile

def predict(datadirectory, model, outputfile):
    basename = os.path.basename(os.path.normpath(datadirectory))
    model = load_model(model)
    
    featurefile = "./PEAS/features.txt"
    labelencoder = "./PEAS/labelencoder.txt"

    #Step 1: Load the data
    start_time = time.time()
    seqdata,sigdata,annot,summitpeaks,peaks = tl.readTensors(basename, datadirectory, 600, sequence=True, signal=True)
    peasfeatures = tl.getPEASFeatures(datadirectory+"/peak_features/"+basename+"_features.txt", featurefile, labelencoder, peaks)
    num_classes = 4
    peasfeatures = np.expand_dims(peasfeatures, axis=2)
    sigseqdata = tl.getSeqSigTensor(seqdata, sigdata)

    print("--- Data loaded in %s seconds ---" % (time.time() - start_time))

    x_test_sigseq = sigseqdata
    x_test_sigseq = np.moveaxis(x_test_sigseq, 1, -1) #Originally had channels first, but CPU tensorflow requires channels last
    x_test_peas = peasfeatures

    #Step 2: Make predictions
    start_time = time.time()
    sig_predictions, peas_predictions, predictions = model.predict([x_test_sigseq, x_test_peas])

    print("--- Data predicted in %s seconds ---" % (time.time() - start_time))

    #Write the output file:
    columns = ["Chr", "Start", "End", "Promoter Probability", "Enhancer Probability", "Insulator Probability", "Other Probability"]
    pd.DataFrame(np.concatenate((peaks, predictions), axis=1), columns=columns).to_csv(outputfile, header=None, index=None, sep="\t")


predict(datadirectory, model, outputfile)
