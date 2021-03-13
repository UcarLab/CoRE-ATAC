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
epochs = 20
args = sys.argv

sig_shape = (600,10)
num_classes = 4

featurefile = "/CoRE-ATAC/PEAS/features.txt"
labelencoder = "/CoRE-ATAC/PEAS/labelencoder.txt"


#Step 0: Process arguments
parser = argparse.ArgumentParser(description='CoRE-ATAC Model Trainer - SigSeq')
parser.add_argument("basenames")
parser.add_argument("datadirectories")
parser.add_argument("trainchrs")
parser.add_argument("outputmodel")


args = parser.parse_args()

basef = args.basenames
dirf = args.datadirectories
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
trainchr, valchr, testchr = getTrainChr(trainchrf)


bestmodelsaver = ModelCheckpoint(outmodelf, monitor='val_output_loss', verbose=1, save_best_only=True, mode='min')
callbacks = [bestmodelsaver]

###################
#Step 1: Load Data#
###################

seqdata,sigdata,annot,summitpeaks,peaks = tl.readTensors(datanames[0], datadirectories[0], 600, sequence=True, signal=True)

testindices = part.chrSelect(peaks, valchr)
ftestindices = part.chrSelect(peaks, testchr)
trainindices = part.chrSelect(peaks, trainchr)

sigdata = tl.getSeqSigTensor(seqdata, sigdata)

x_train_sig = sigdata[trainindices]
y_train_0 = annot[trainindices]

x_test_sig = sigdata[testindices]
y_test_0 = annot[testindices]

x_ftest_sig = sigdata[ftestindices]
y_ftest_0 = annot[ftestindices]

y_train = y_train_0
y_test = y_test_0
y_ftest = y_ftest_0


print(datanames[0]+" loaded.")

for curindex in range(1, len(datadirectories)):
    print("Begin loading "+datanames[curindex]+".")

    seqdata,sigdata,annot,summitpeaks,peaks = tl.readTensors(datanames[curindex], datadirectories[curindex], 600, sequence=True, signal=True)
                                                                                     
    testindices = part.chrSelect(peaks, valchr)
    ftestindices = part.chrSelect(peaks, testchr)
    trainindices = part.chrSelect(peaks, trainchr)

    sigdata = tl.getSeqSigTensor(seqdata, sigdata)


    x_train_sig = np.concatenate((x_train_sig, sigdata[trainindices]), axis=0)
    y_train_0 = annot[trainindices]
    
    x_test_sig = np.concatenate((x_test_sig, sigdata[testindices]), axis=0)
    y_test_0 = annot[testindices]

    x_ftest_sig = np.concatenate((x_ftest_sig, sigdata[ftestindices]), axis=0)
    y_ftest_0 = annot[ftestindices]

    y_train = np.concatenate((y_train, y_train_0), axis=0)
    y_test = np.concatenate((y_test, y_test_0), axis=0)
    y_ftest = np.concatenate((y_ftest, y_ftest_0), axis=0)

    print(datanames[curindex]+" loaded.")



#####################
#Step 2: Build Model#
#####################

#Inputs
sig_input = Input(shape=sig_shape, dtype='float', name='sig_input')


#Signal
sig_conv1 = SeparableConv1D(name='sig_conv_1',
                 filters=256, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_input)
sig_batch1 = BatchNormalization(name='sig_batch_normalization1', axis=-1)(sig_conv1)

sig_conv2 = SeparableConv1D(name="sig_conv_2",
                 filters=256, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch1)
sig_maxpool1 = MaxPooling1D(pool_size=2, name='sig_maxpool1',
                 data_format='channels_last')(sig_conv2)
sig_batch2 = BatchNormalization(name='sig_batch_normalization2', axis=-1)(sig_maxpool1)


sig_conv3 = SeparableConv1D(name='sig_conv_3',
                 filters=256, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch2)
sig_batch3 = BatchNormalization(name='sig_batch_normalization3', axis=-1)(sig_conv3)

sig_conv4 = SeparableConv1D(name="sig_conv_4",
                 filters=256, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch3)
sig_maxpool2 = MaxPooling1D(pool_size=2, name='sig_maxpool2',
                 data_format='channels_last')(sig_conv4)
sig_batch4 = BatchNormalization(name='sig_batch_normalization4', axis=-1)(sig_maxpool2)


sig_conv5 = SeparableConv1D(name='sig_conv_5',
                 filters=512, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch4)
sig_batch5 = BatchNormalization(name='sig_batch_normalization5', axis=-1)(sig_conv5)

sig_conv6 = SeparableConv1D(name="sig_conv_6",
                 filters=512, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch5)
sig_maxpool3 = MaxPooling1D(pool_size=2, name='sig_maxpool3',
                 data_format='channels_last')(sig_conv6)
sig_batch6 = BatchNormalization(name='sig_batch_normalization6', axis=-1)(sig_maxpool3)


#
sig_conv7 = SeparableConv1D(name='sig_conv_7',
                 filters=512, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch6)
sig_batch7 = BatchNormalization(name='sig_batch_normalization7', axis=-1)(sig_conv7)

sig_conv8 = SeparableConv1D(name="sig_conv_8",
                 filters=512, kernel_size=19, activation='relu',
                 data_format='channels_last')(sig_batch7)
sig_maxpool4 = MaxPooling1D(pool_size=2, name='sig_maxpool4',
                 data_format='channels_last')(sig_conv8)
sig_batch8 = BatchNormalization(name='sig_batch_normalization8', axis=-1)(sig_maxpool4)


sig_flatten = Flatten(name='sig_flatten')(sig_batch8)
sig_dense = Dense(2048, name='sig_dense', activation='relu')(sig_flatten)

output = Dense(num_classes, activation='softmax', name='output')(sig_dense) 


model = Model(inputs=[sig_input], outputs=[output])
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
        #cuts_to_dense = listToDense(X[3][subset], mat_len, mat_len)
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

x_train_sig = np.moveaxis(x_train_sig, 1, -1)
x_test_sig = np.moveaxis(x_test_sig, 1, -1)

steps_per_epoch = int(np.ceil(x_train_sig.shape[0]/batch_size))
generator = sparse_to_dense_generator([x_train_sig], y_train, batch_size=batch_size, shuffle=True, rstate=929)

class_weights = class_weight.compute_class_weight('balanced', np.unique(y_train), y_train)
class_weights_dict = dict()
for i in range(num_classes):
    class_weights_dict[i] = class_weights[i]

model.fit_generator(generator,
          steps_per_epoch=steps_per_epoch,
          epochs=epochs, callbacks=callbacks,
          verbose=1, class_weight=class_weights_dict,
          validation_data=([x_test_sig], y_test))


