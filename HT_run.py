from pickle import load

def dimOfTup(tup):
    tp = 1
    for l in tup: tp *= 2*l+1
    return tp

energyOfState = lambda tup : sum([l+1/2 for l in tup])


fn = "tuples_20.pickle"
with open(fn,'rb') as f:
    cutoff, evenIrrepList, oddIrrepList = load(f)

print("A list of the number of all even-parity irreps below cutoff {}:".format(cutoff))
countRep = lambda L : [x[1].get(L,0) for x in evenIrrepList]
evenIrrepTable= [(L,sum(countRep(L))) for L in range(0,19)]
print("L   | # of irreps")
for x in evenIrrepTable:
    print("{:<3} | {:<10}".format(x[0],x[1]))
totStates = sum([(2*x[0]+1)*x[1] for x in evenIrrepTable])
print("The total number of states is {}.".format(totStates))

mTable = []
for M in range(0,19):
    ct = 0
    for L,N in evenIrrepTable:
        if L >= abs(M):
            ct += N
    mTable.append((M,ct))
print("A list of the number of all even-parity states with Jz = M:")
print("M   | # of states")
for x in mTable:
    print("{:<3} | {:<10}".format(x[0],x[1]))
print("\n")
tot2 = mTable[0][1] + sum([2*x[1] for x in mTable[1:]])
if tot2 != totStates: print("Mismatch in the number of states.")


print("A list of the number of all odd-parity irreps below cutoff {}:".format(cutoff))
countRep = lambda L : [x[1].get(L,0) for x in oddIrrepList]
oddIrrepTable= [(L,sum(countRep(L))) for L in range(0,20)]
print("L   | # of irreps")
for x in oddIrrepTable:
    print("{:<3} | {:<10}".format(x[0],x[1]))
totStates = sum([(2*x[0]+1)*x[1] for x in oddIrrepTable])
print("The total number of parity-odd states is {}.".format(totStates))

mTable = []
for M in range(0,20):
    ct = 0
    for L,N in oddIrrepTable:
        if L >= abs(M):
            ct += N
    mTable.append((M,ct))
print("A list of the number of all parity-odd states with Jz = M:")
print("M   | # of states")
for x in mTable:
    print("{:<3} | {:<10}".format(x[0],x[1]))
print("\n")
tot2 = mTable[0][1] + sum([2*x[1] for x in mTable[1:]])
if tot2 != totStates: print("Mismatch in the number of states.")
        
L = 0
tupList = []
for x in evenIrrepList:
    t, n = x[0], x[1].get(L,0)
    if n > 0: tupList.append((t,n))
#print(tupList)

noZeroes = [x for x in tupList if x[0].count(0) == 0]
#print("All tensor products that don't contain factors of 0:", noZeroes)
tups = [x[0] for x in noZeroes]
stripOnes = lambda x : tuple(filter((1).__ne__,x[0]))
noOnes = set(map(stripOnes,noZeroes))
noOnes = sorted([(tup,dimOfTup(tup)) for tup in noOnes],key = lambda x : x[1])
#print(noOnes)
