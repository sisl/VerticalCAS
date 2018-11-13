import numpy as np
import theano
import theano.tensor as T
import math
import keras
from keras.models import Sequential, load_model
from keras.layers import Dense, Dropout, Activation
import h5py
from keras.optimizers import Adamax, Nadam
import sys
from writeNNet import saveNNet

######## OPTIONS #########
ver = 4            # Neural network version
hu = 45            # Number of hidden units in each hidden layer in network
saveEvery = 10     # Epoch frequency of saving
totalEpochs = 200  # Total number of training epochs
trainingDataFiles = "../TrainingData/VertCAS_TrainingData_v2_%02d.h5" # File format for training data
nnetFiles = "../networks/VertCAS_pra%02d_v%d_45HU_%03d.nnet" # File format for .nnet files
##########################

# The previous RA should be given as a command line input
if len(sys.argv) > 1:
    pra = int(sys.argv[1])
    print "Loading Data for VertCAS, pra %02d, Network Version %d" % (pra, ver)
    f       = h5py.File(trainingDataFiles % pra,'r')
    X_train = np.array(f['X'])
    Q       = np.array(f['y'])
    means = np.array(f['means'])
    ranges=np.array(f['ranges'])
    min_inputs = np.array(f['min_inputs'])
    max_inputs = np.array(f['max_inputs'])
                       
    N,numOut = Q.shape
    print "Setting up Model"
    
    # Asymmetric loss function
    lossFactor = 40.0
    def asymMSE(y_true, y_pred):
        d = y_true-y_pred
        maxes = T.argmax(y_true,axis=1)
        maxes_onehot = T.extra_ops.to_one_hot(maxes,numOut)
        others_onehot = maxes_onehot-1
        d_opt = d*maxes_onehot 
        d_sub = d*others_onehot
        a = lossFactor*(numOut-1)*(d_opt**2+keras.backend.abs(d_opt))
        b = d_opt**2
        c = lossFactor*(d_sub**2+keras.backend.abs(d_sub))
        d = d_sub**2
        loss = T.switch(d_sub>0,c,d) + T.switch(d_opt>0,a,b)
        return loss

    # Define model architecture
    model = Sequential()
    model.add(Dense(hu, init='uniform', activation='relu', input_dim=4))
    model.add(Dense(hu, init='uniform', activation='relu'))
    model.add(Dense(hu, init='uniform', activation='relu'))
    model.add(Dense(hu, init='uniform', activation='relu'))
    model.add(Dense(hu, init='uniform', activation='relu'))
    model.add(Dense(hu, init='uniform', activation='relu'))
    model.add(Dense(numOut,   init='uniform'))
    opt = Nadam(lr = 0.0003)
    model.compile(loss=asymMSE, optimizer=opt, metrics=['accuracy'])

    # Train and write nnet files
    epoch= saveEvery
    while epoch <= totalEpochs:
        model.fit(X_train, Q, nb_epoch=saveEvery, batch_size=2**8,shuffle=True)
        saveFile = nnetFiles % (pra, ver,epoch)
        saveNNet(model,saveFile,means,ranges,min_inputs,max_inputs)
        epoch += saveEvery
