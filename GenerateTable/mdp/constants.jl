export COC,DNC,DND,DES1500,CL1500,SDES1500,SCL1500,SDES2500,SCL2500, stateType, actType, acts, discount_f, hMin, hMax,hs, vMin, vMax, vowns, vints, interp, accels, velRanges, allowedTrans

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
stateType = Tuple{Float64,Float64,Float64,Int}
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
vowns = vels
vints = vels 

# Create the local function approximator using the grid
interp = LocalGIFunctionApproximator(RectangleGrid(hs,vowns,vints,acts))  


### Dictionaries to define transitions ###
# Tuple of probabilities and corresponding acceleration in ft/s^2
probs = [0.5,0.3,0.2] #[0.5,0.25,0.25]
g = 32.2
accels = Dict(COC=>([0.34,0.33,0.33],[0.0,-g/3,g/3]),
              DNC=>(probs,[-g/3,-g/2,g/3]),
              DND=>(probs,[g/3,g/2,-g/3]),
              DES1500=>(probs,[-g/3,-g/2,g/3]),
              CL1500=>(probs,[g/3,g/2,-g/3]),
              SDES1500=>(probs,[-g/2.5,-g/2,g/3]),
              SCL1500=>(probs,[g/2.5,g/2,-g/3]),
              SDES2500=>(probs,[-g/2.5,-g/2,g/3]),
              SCL2500=>(probs,[g/2.5,g/2,-g/3]),
              -1=>([0.34,0.33,0.33],[0,-g/8,g/8]))

# Velocity range where aircraft is NON-compliant with advisory (ft/s)
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
allowedTrans = Dict(COC     => [1,1,1,1,1,0,0,0,0],
                   DNC      => [1,1,1,0,0,1,1,0,0],
                   DND      => [1,1,1,0,0,1,1,0,0],
                   DES1500  => [1,1,1,1,0,0,1,1,0],
                   CL1500   => [1,1,1,0,1,1,0,0,1],
                   SDES1500 => [1,1,1,0,0,1,1,1,0],
                   SCL1500  => [1,1,1,0,0,1,1,0,1],
                   SDES2500 => [1,1,1,0,0,0,1,1,0],
                   SCL2500  => [1,1,1,0,0,1,0,0,1])