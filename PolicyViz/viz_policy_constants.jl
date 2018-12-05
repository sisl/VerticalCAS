export relhs,dh0s,dh1s,taus,pas,nstates,actions, action_names

# State dimension cutpoints
const vels = vcat(linspace(-100,-60,5),linspace(-50,-35,4),linspace(-30,30,21),linspace(35,50,4),linspace(60,100,5))
const relhs   = vcat(linspace(-8000,-4000,5),linspace(-3000,-1250,8),linspace(-1000,-800,3),linspace(-700,-150,12),linspace(-100,100,9), linspace(150,700,12),linspace(800,1000,3),linspace(1250,3000,8),linspace(4000,8000,5))
const dh0s = vels
const dh1s = vels
const taus  = linspace(0,40,41)
const pas = [1,2,3,4,5,6,7,8,9]

# Number of states
const nstates = length(relhs)*length(dh0s)*length(taus)*length(dh1s)

# Actions
const actions = [1, 2, 3, 4, 5, 6, 7, 8, 9]
const action_names = ["COC","DNC","DND","DES1500","CL1500","SDES1500","SCL1500","SDES2500","SCL2500"]