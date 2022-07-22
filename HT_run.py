from pickle import load

#################################################

def dimOfTup(tup):
    tp = 1
    for l in tup: tp *= 2*l+1
    return tp

uenergyOfTuple = lambda tup : sum([l+1/2 for l in tup])

#################################################

fn = "tuples_20.pickle"
with open(fn,'rb') as f:
    cutoff, evenIrrepList, oddIrrepList = load(f)
maxSpin = int(cutoff)+1

countRep = lambda L : [x[1].get(L,0) for x in evenIrrepList]
evenIrrepTable= [(L,sum(countRep(L))) for L in range(0,maxSpin)]
mTableEven = []
for M in range(0,maxSpin):
    ct = 0
    for L,N in evenIrrepTable:
        if L >= abs(M):
            ct += N
    mTableEven.append((M,ct))
evenZipped = zip(evenIrrepTable, mTableEven)
print("A list of the number of all even-parity irreps below cutoff {}:\n".format(cutoff))
print("{:>3} {:>10} | {:<10}".format("J", "# irreps", "# of states with Jz = J"))
for x in evenZipped:
    (L, nIr), (_, nM) = x
    print("{:>3} {:>10} | {:>7}".format(L, nIr, nM))
totStatesEven = sum([(2*x[0]+1)*x[1] for x in evenIrrepTable])
print("The total number of parity-even states is {}.\n".format(totStates))
tot2 = mTableEven[0][1] + sum([2*x[1] for x in mTableEven[1:]])
if tot2 != totStatesEven: print("Mismatch in the number of states.")

##############
### If you care about odd-parity states, uncomment this:
##############
#
# countRep = lambda L : [x[1].get(L,0) for x in oddIrrepList]
# oddIrrepTable= [(L,sum(countRep(L))) for L in range(0,maxSpin)]
# mTableOdd = []
# for M in range(0,maxSpin):
#     ct = 0
#     for L,N in oddIrrepTable:
#         if L >= abs(M):
#             ct += N
#     mTableOdd.append((M,ct))
# oddZipped = zip(oddIrrepTable, mTableOdd)
# print("A list of the number of all odd-parity irreps below cutoff {}:\n".format(cutoff))
# print("{:>3} {:>10} | {:<10}".format("J", "# irreps", "# of states with Jz = J"))
# for x in oddZipped:
#     (L, nIr), (_, nM) = x
#     print("{:>3} {:>10} | {:>7}".format(L, nIr, nM))
# totStatesOdd = sum([(2*x[0]+1)*x[1] for x in oddIrrepTable])
# print("The total number of parity-odd states is {}.\n".format(totStatesOdd))
# tot2 = mTableOdd[0][1] + sum([2*x[1] for x in mTableOdd[1:]])
# if tot2 != totStates: print("Mismatch in the number of states.")

##################################

# set an effective cutoff and a spin L, and take all tuples below the effective cutoff
# that have an irrep of spin L
effCutoff = 18
L = 0
retainedTuples = []
for x in evenIrrepList:
    t, n = x[0], x[1].get(L,0)
    if n > 0 and energyOfState(t) <= effCutoff: retainedTuples.append((t,n))
    
stripZeroes = lambda tup : tuple(filter((0).__ne__,tup))
stripOnes = lambda tup : tuple(filter((1).__ne__,tup))
# a list of all tuples with 0s removed
noZeroes = {stripZeroes(x[0]) for x in retainedTuples}
print(noZeroes)

#print("All tensor products that don't contain factors of 0:", noZeroes)
maxFactorsOf1 = max({x[0].count(1) for x in retainedTuples})
print("There are at most {} factors of 1.".format(maxFactorsOf1))

# set of tuples without 0s and 1s
noOnes = set(map(stripOnes,noZeroes))
noOnes = sorted([(tup,dimOfTup(tup)) for tup in noOnes],key = lambda x : x[1])
print(noOnes)
