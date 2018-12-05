import h5py
import numpy as np

###### OPTIONS #######
tableFile = '../GenerateTable/VertCAS_qvals_parallel_v7_tauMax40.h5'
trainingDataFiles = '../TrainingData/VertCAS_TrainingData_v2_%02d.h5'
######################

# Define state space. Make sure this matches up with the constants used to generate the MDP table!
acts = [0,1,2,3,4,5,6,7,8]
vels = np.concatenate((np.linspace(-100,-60,5),np.linspace(-50,-35,4),np.linspace(-30,30,21),np.linspace(35,50,4),np.linspace(60,100,5)))
hs   = np.concatenate((np.linspace(-8000,-4000,5),np.linspace(-3000,-1250,8),np.linspace(-1000,-800,3),np.linspace(-700,-150,12),np.linspace(-100,100,9),np.linspace(150,700,12),np.linspace(800,1000,3),np.linspace(1250,3000,8),np.linspace(4000,8000,5)))
vowns = vels
vints = vels 
taus  = np.linspace(0,40,41)

# Get table cutpoints
X = np.array([[h,vo,vi,t] for t in taus for vi in vints for vo in vowns for h in hs])

# Compile table values
f = h5py.File(tableFile,'r')
Q = np.array(f['q'])
f.close()
Q = Q.T
ns2 = len(hs)*len(vowns)*len(vints)*len(acts)*2
Qtaus = [Q[i*ns2:(i*ns2+ns2/2)] for i in range(len(taus))]
ns = len(hs)*len(vowns)*len(vints)
Qacts = [np.concatenate([Qtaus[i][ns*j:(j+1)*ns] for i in range(len(taus))]) for j in range(len(acts))]

# Compute means, ranges, mins and maxes
means = np.mean(X,axis=0)
ranges = np.max(X,axis=0)-np.min(X,axis=0)
X =(X-means)/ranges

meanQ = np.mean([np.mean(Qacts[i]) for i in range(len(Qacts))])
rangeQ = np.max([np.max(Qacts[i]) for i in range(len(Qacts))])-np.min([np.min(Qacts[i]) for i in range(len(Qacts))]) 
Qacts = [(Qacts[i]-meanQ)/rangeQ for i in range(len(Qacts))]

means = np.concatenate((means,[meanQ]))
ranges = np.concatenate((ranges,[rangeQ]))

min_inputs = np.array([hs[0],vowns[0],vints[0],taus[0]])
max_inputs = np.array([hs[-1],vowns[-1],vints[-1],taus[-1]])

#Save the Training Data
for ra in range(1,10):
  with h5py.File(trainingDataFiles%ra,'w') as H:
      H.create_dataset('X',data=X)
      H.create_dataset('y',data=Qacts[ra-1])
      H.create_dataset('means',data=means)
      H.create_dataset('ranges',data=ranges)
      H.create_dataset('min_inputs',data=min_inputs)
      H.create_dataset('max_inputs',data=max_inputs)