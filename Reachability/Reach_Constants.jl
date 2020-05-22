# Advisory indices
COC = 0
DNC = 1
DND = 2
DES1500 = 3
CL1500 = 4
SDES1500 = 5
SCL1500 = 6
SDES2500 = 7
SCL2500 = 8

ACTIONS = [COC,DNC,DND,DES1500,CL1500,SDES1500,SCL1500,SDES2500,SCL2500]
G = 32.2

# Grid discretization
HS = vcat(range(-3000,stop=-1600,length=15),
                       range(-1500,stop=-950,length=12),
                       range(-900,stop=-725,length=8),
                       range(-700,stop=-510,length=20),
                       range(-500,stop=-205,length=60),
                       range(-200,stop=-53,length=50),
                       range(-50,stop=50,length=51),
                       range(53,stop=200,length=50),
                       range(205,stop=500,length=60),
                       range(510,stop=700,length=20),
                       range(725,stop=900,length=8),
                       range(950,stop=1500,length=12),
                       range(1600,stop=3000,length=15))

HS = vcat(range(-3000,stop=-1600,length=15),
                       range(-1500,stop=-950,length=12),
                       range(-900,stop=-725,length=8),
                       range(-700,stop=-510,length=20),
                       range(-500,stop=-203,length=100),
                       range(-200,stop=-102,length=50),
                       range(-100,stop=100,length=201),
                       range(102,stop=200,length=50),
                       range(203,stop=500,length=100),
                       range(510,stop=700,length=20),
                       range(725,stop=900,length=8),
                       range(950,stop=1500,length=12),
                       range(1600,stop=3000,length=15))

VOWNS = vcat([-41.6667,-41.666],
       range(-41,stop=-26,length=16) , #5
       [-25.001,-24.999],
       range(-24,stop=-1,length=24), #12
       [-0.001,0.001],
       range(1,stop=24,length=24), #12
       [24.999,25.001],
       range(26,stop=41,length=16), #5
       [41.666,41.667])

#vintVec = LinRange(-25,25,17)
VINTS = LinRange(-25,25,51)

TAUS = LinRange(0,40,41)

NUMH = Int32(length(HS)-1)
NUMVOWN = Int32(length(VOWNS)-1)
NUMVINT = Int32(length(VINTS)-1)
NUMACTION= Int32(length(ACTIONS))
NUMTAU = Int32(length(TAUS))
NUMREGIONS = NUMH*NUMVOWN*NUMVINT