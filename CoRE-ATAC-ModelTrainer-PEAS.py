from __future__ import print_function
import keras
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten, Input
from keras.layers import Conv1D, Conv2D, MaxPooling1D, MaxPooling2D, BatchNormalization, Activation, SeparableConv1D
from keras.callbacks import ModelCheckpoint
from keras import backend as K
from tensorloader import TensorLoader as tl
import partitioning_util as part
import sys
import numpy as np 
import pandas as pd
import pickle
import h5py
from sklearn.utils import class_weight
import argparse
import os


batch_size = 32
epochs = 100
args = sys.argv
seq_shape = (4,600)
sig_shape = (5,600)
sigshape_shape = (50,600,1)
cut_shape = (600,600,1)
peas_shape = (30,1)
num_classes = 4

featurefile = "/CoRE-ATAC/PEAS/features.txt"
labelencoder = "/CoRE-ATAC/PEAS/labelencoder.txt"


#Step 0: Process arguments
parser = argparse.ArgumentParser(description='CoRE-ATAC Model Trainer - PEAS')
parser.add_argument("basenames")
parser.add_argument("datadirectories")
parser.add_argument("peasfiles")
parser.add_argument("trainchrs")
parser.add_argument("outputmodel")


args = parser.parse_args()

basef = args.basenames
dirf = args.datadirectories
peasff =  args.peasfiles
trainchrf = args.trainchrs
outmodelf =  args.outputmodel

def readList(file):
    return list(pd.read_csv(file, header=None).values.flatten())

def getTrainChr(file):
    data = pd.read_csv(file, header=None, sep="\t").values
    
    trainchr = []
    valchr = []
    testchr = []
    
    for i in range(len(data)):
        if(data[i,0] == "train"):
            trainchr.extend(data[i,1].split(","))
        if(data[i,0] == "val"):
            valchr.extend(data[i,1].split(","))
        if(data[i,0] == "test"):
            testchr.extend(data[i,1].split(","))
            
    return trainchr, valchr, testchr


datanames = readList(basef)
datadirectories = readList(dirf)
peasdata = readList(peasff)
trainchr, valchr, testchr = getTrainChr(trainchrf)

bestmodelsaver = ModelCheckpoint(outmodelf, monitor='val_output_loss', verbose=1, save_best_only=True, mode='min')
callbacks = [bestmodelsaver]


###################
#Step 1: Load Data#
###################

_,_,annot,summitpeaks,peaks = tl.readTensors(datanames[0], datadirectories[0], 600, sequence=False, signal=False)
peasfeatures = tl.getPEASFeatures(peasdata[0], featurefile, labelencoder, peaks)

testindices = part.chrSelect(peaks, valchr)
ftestindices = part.chrSelect(peaks, testchr)
trainindices = part.chrSelect(peaks, trainchr)
peasfeatures = np.expand_dims(peasfeatures, axis=2)

x_train_peas = peasfeatures[trainindices]
y_train_0 = annot[trainindices]

x_test_peas = peasfeatures[testindices]
y_test_0 = annot[testindices]

x_ftest_peas = peasfeatures[ftestindices]
y_ftest_0 = annot[ftestindices]

y_train = y_train_0
y_test = y_test_0
y_ftest = y_ftest_0

for curindex in range(1, len(datadirectories)):
    _,_,annot,summitpeaks,peaks = tl.readTensors(datanames[curindex], datadirectories[curindex], 600, sequence=False, signal=False)
    peasfeatures = tl.getPEASFeatures(peasdata[curindex], featurefile, labelencoder, peaks)
                                                                                     
    testindices = part.chrSelect(peaks, valchr)
    ftestindices = part.chrSelect(peaks, testchr)
    trainindices = part.chrSelect(peaks, trainchr)

    peasfeatures = np.expand_dims(peasfeatures, axis=2)

    x_train_peas = np.concatenate((x_train_peas,peasfeatures[trainindices]), axis=0)
    y_train_0 = annot[trainindices]
    
    x_test_peas = np.concatenate((x_test_peas, peasfeatures[testindices]), axis=0)
    y_test_0 = annot[testindices]

    x_ftest_peas = np.concatenate((x_ftest_peas, peasfeatures[ftestindices]), axis=0)
    y_ftest_0 = annot[ftestindices]

    y_train = np.concatenate((y_train, y_train_0), axis=0)
    y_test = np.concatenate((y_test, y_test_0), axis=0)
    y_ftest = np.concatenate((y_ftest, y_ftest_0), axis=0)



peas_shape = (np.shape(x_train_peas)[1],1)

#####################
#Step 2: Build Model#
#####################

#Inputs
peas_input = Input(shape=peas_shape, dtype='float', name='peas_input')

#PEAS Features
peas_flatten = Flatten(name='peas_flatten')(peas_input)

peas_dense = Dense(2048, name='peas_dense', activation='relu')(peas_flatten)


output = Dense(num_classes, activation='softmax', name='output')(peas_dense) 


model = Model(inputs=[peas_input], outputs=[output])
model.compile(loss=keras.losses.sparse_categorical_crossentropy,
              optimizer=keras.optimizers.Adam(),
              metrics=['accuracy'])


######################################################
#Step 3: Create generator for reading sparse matrices#
######################################################

def listToDense(positions, nrows, ncols):
    num_samples = len(positions)
    rv = np.zeros((num_samples, nrows, ncols,1), dtype=bool)
    for i in range(num_samples):
        for coord in positions[i]:
            if len(coord) == 2:
                rv[i,coord[0],coord[1],0] = 1
    return rv

def sparse_to_dense_generator(X, y, batch_size, shuffle, rstate=None):
    num_samples = X[0].shape[0]
    mat_len = X[0].shape[2]
    num_batches = int(num_samples/batch_size) + 1 - int(num_samples % batch_size == 0)
    idx = np.arange(num_samples)
    
    rs = np.random.RandomState(seed=rstate)
    if shuffle:
        rs.shuffle(idx)
    
    curbatch = 0 
    while True:
        subset = idx[batch_size*curbatch:batch_size*(curbatch+1)]
        rv_X = [X[0][subset]]
        rv_y = y[subset]
        yield rv_X, rv_y
        curbatch += 1
        if curbatch == num_batches:
            if shuffle:
                rs.shuffle(idx)
            curbatch = 0



#####################
#Step 4: Train Model#
#####################

steps_per_epoch = int(np.ceil(x_train_peas.shape[0]/batch_size))
generator = sparse_to_dense_generator([x_train_peas], y_train, batch_size=batch_size, shuffle=True, rstate=929)

print(np.shape(x_train_peas))

class_weights = class_weight.compute_class_weight('balanced', np.unique(y_train), y_train)
class_weights_dict = dict()
for i in range(num_classes):
    class_weights_dict[i] = class_weights[i]

model.fit_generator(generator,
          steps_per_epoch=steps_per_epoch,
          epochs=epochs, callbacks=callbacks,
          verbose=1, class_weight=class_weights_dict,
          validation_data=([x_test_peas], y_test))


