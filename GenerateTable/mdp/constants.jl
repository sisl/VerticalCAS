export COC,DNC,DND,DES1500,CL1500,SDES1500,SCL1500,SDES2500,SCL2500, stateType, actType, acts, discount_f, hMin, hMax,hs, vMin, vMax, vowns, vints, tauMin, tauMax, numTau, taus, resps, interp, accels, velRanges, allowedTrans

# ADVISORY INDICES
COC=0
DNC=1
DND=2
DES1500 = 3
CL1500 = 4
SDES1500=5
SCL1500=6
SDES2500=7
SCL2500=8

# State Type:
stateType = Tuple{Float64,Float64,Float64,Int,Bool}
actType = Int
acts = [0,1,2,3,4,5,6,7,8]

# Default parameters
discount_f = 1.0

### STATE CUTPOINTS ###
# H: Altitude
hMin = -8000.0
hMax = 8000.0
hs   = vcat(LinRange(-8000,-4000,5),LinRange(-3000,-1250,8),LinRange(-1000,-800,3),LinRange(-700,-150,12),LinRange(-100,100,9),LinRange(150,700,12),LinRange(800,1000,3),LinRange(1250,3000,8),LinRange(4000,8000,5))

# Velocities
vMin = -100.
vMax = 100.
vels = vcat(LinRange(-100,-60,5),LinRange(-50,-35,4),LinRange(-30,30,21),LinRange(35,50,4),LinRange(60,100,5))
vowns = vels #LinRange(vMin,vMax,6)
vints = vels #LinRange(vMin,vMax,6)

# Taus
tauMin = 0.0
tauMax = 40.0
numTau = 41
taus  = LinRange(tauMin,tauMax,numTau)

# Responding
resps = [false,true]

interp = LocalGIFunctionApproximator(RectangleGrid(hs,vowns,vints,acts,resps))  # Create the local function approximator using the grid


### Dictionaires to define transitions ###
# Tuple of probabilities and corresponding acceleration in ft/s^2
accels = Dict(COC=>([0.5,0.25,0.25],[0.0,3.0,-3.0]),
              DNC=>([0.5,0.25,0.25],[-8.33,-9.33,-7.33]),
              DND=>([0.5,0.25,0.25],[8.33,9.33,7.33]),
              DES1500=>([0.5,0.25,0.25],[-8.33,-9.33,-7.33]),
              CL1500=>([0.5,0.25,0.25],[8.33,9.33,7.33]),
              SDES1500=>([0.5,0.25,0.25],[-10.7,-11.7,-9.7]),
              SCL1500=>([0.5,0.25,0.25],[10.7,11.7,9.7]),
              SDES2500=>([0.5,0.25,0.25],[-10.7,-11.7,-9.7]),
              SCL2500=>([0.5,0.25,0.25],[10.7,11.7,9.7]))

# Velocity range where aircraft is NON-compliant with advisory
velRanges = Dict(COC=>(-100.0,100.0),
                DNC=>(0.0,100.0),
                DND=>(-100.0,0.0),
                DES1500=>(-25.0,100.0),
                CL1500=>(-100.0,25.0),
                SDES1500=>(-25.0,100.0),
                SCL1500=>(-100.0,25.0),
                SDES2500=>(-41.67,100.0),
                SCL2500=>(-100.0,41.67))

# Allowed transitions between advisories
allowedTrans = Dict(COC=>[1,1,1,1,1,0,0,0,0],
                   DNC=>[1,1,1,1,1,0,0,0,0],
                   DND=>[1,1,1,1,1,0,0,0,0],
                   DES1500=>[1,1,1,1,1,1,1,0,0],
                   CL1500=>[1,1,1,1,1,1,1,0,0],
                   SDES1500=>[1,1,1,1,1,1,1,1,1],
                   SCL1500=>[1,1,1,1,1,1,1,1,1],
                   SDES2500=>[1,1,1,1,1,1,1,1,1],
                   SCL2500=>[1,1,1,1,1,1,1,1,1])