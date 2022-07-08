import chains
from bases import genAllData
from sympy import Matrix
import pickle
#from group_theory import countDimensionProduct

clebschDict = {}
## if you want to load the precomputed CG coefficients:
with open("clebsch_memo.pkl",'rb') as f:
    clebschDict = pickle.load(f)
    
spins = (0,0,0,0)
lMax = sum(spins)

rep = lambda l : range(-l,l+1)

numStates = sum([rawDec.get(L,0)*(2*L+1) for L in range(0,lMax+1)])
numStates2 = 1
for l in spins: numStates2 *= 2*l+1
#print(numStates, numStates2)
print("We are expecting {} tensors in total.".format(numStates))
if numStates != numStates2: print("Mismatch in the number of states.")

for L in range(0,lMax+1):
    tensorLabels = chains.generateChains(spins,L)
    print("Computing spin {}; there should be {} states.".format(L,rawDec.get(L,0)))
    for M in rep(L):
        #print(L,M)
        [[rawStates, rawLookup, rawDec], [symmStates, symmLookup, symmDec], rawToSymmMatrix] = genAllData(spins,M)
        if len(tensorLabels) != rawDec.get(L,0):
            print("Mismatch.")
        #print("Predicting {} spin-{} multiplets; measuring {}.".format(len(tensorLabels),L,rawDec.get(L,0)))
        coefTable = [[compCoef(chain,state) for state in rawStates] for chain in tensorLabels]
        #mat = Matrix(coefTable)
        #print(mat)


        
# numStates = sum([rawDec.get(L,0)*(2*L+1) for L in range(0,lMax+1)])
# numStates2 = 1
# for l in spins: numStates2 *= 2*l+1
# #print(numStates, numStates2)
# print("We are expecting {} tensors in total.".format(numStates))
# if numStates != numStates2: print("Mismatch in the number of states.")

# for L in range(0,lMax+1):
#     tensorLabels = chains.generateChains(spins,L)
#     print("Computing spin {}; there should be {} states.".format(L,rawDec.get(L,0)))
#     for M in rep(L):
#         #print(L,M)
#         [[rawStates, rawLookup, rawDec], [symmStates, symmLookup, symmDec], rawToSymmMatrix] = genAllData(spins,M)
#         if len(tensorLabels) != rawDec.get(L,0):
#             print("Mismatch.")
#         #print("Predicting {} spin-{} multiplets; measuring {}.".format(len(tensorLabels),L,rawDec.get(L,0)))
#         coefTable = [[compCoef(chain,state) for state in rawStates] for chain in tensorLabels]
#         #mat = Matrix(coefTable)
#         #print(mat)


# l_max = sum(spins)
# M = -1 ### needs to obey abs(M) <= l_max
# basis, lookUp = spinBasis(spins,M,False)

# l_range = range(abs(M),l_max+1)
# ct = 0
# mat = []
# for L in l_range:
    
#     ct += len(tensorLabels)
#     print("There are {0} states in this sector, and {1} different tensors with L = {2}.".format(len(basis), len(tensorLabels),L))
#     print("The different tensors are labeled as follows: \n",tensorLabels)
#     # print(basis)


#     compCoef = lambda chain, state : chains.chainToCoef(chain,state,M,clebschDict)
#     print("Computing matrix: ")
#     coefTable = [[compCoef(chain,state) for state in basis] for chain in tensorLabels]
#     mat += coefTable
#     print("Done.")
#     # print(len(coefTable), coefTable[0])
#     #coefMat = Matrix(coefTable)
#     #coefMatTranspose = coefMat.transpose()
#     #print("Check that the tensors are orthonormal: ",coefMat*coefMatTranspose == Matrix.eye(len(tensorLabels)))
#     # print("Computed CG coefficients: \n",clebschDict)
# print("Is the number of tensors OK? ",ct == len(basis))

#print("Is this a unitary change of basis? ",mat*mat.transpose() == Matrix.eye(len(basis)))
