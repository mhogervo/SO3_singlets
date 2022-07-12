import chains
from bases import genAllData, casimirMatrix, symmProjector
from sympy import Matrix, ImmutableMatrix, ImmutableSparseMatrix, diag, eye, sqrt
from random import randrange
import pickle
#from group_theory import countDimensionProduct

clebschDict = {}
## if you want to load the precomputed CG coefficients:
#with open("clebsch_memo.pkl",'rb') as f:
#    clebschDict = pickle.load(f)
    
spins = (1,1,1,1,1,1,1,2,2,4)
spins = (1,1,2,2,4)
lMax = sum(spins)
L = 3
M = -3
print("Computing bases etc.")
[rawStates, [symmStates, symmLookup, symmDec, halfOfGram], symmetrizer] = genAllData(spins,M)
print(len(symmStates), symmetrizer.shape)
numSymm = symmDec.get(L,0)
print(symmDec)

#print(symmNorms)
#gramMatrix = ImmutableSparseMatrix(diag(*symmNorms))
#invHalfOfGram = halfOfGram&&
#print(halfOfGram*halfOfGram == gramMatrix)
#proj = symmProjector(rawLookup, ymmLookup) 
print("Number of states: {} vs {}.".format(len(rawStates), len(symmStates)))
print("Number of tensors (among symmetrized states): {}.".format(numSymm))
print("Done.")

print("Computing chains...")
tensorLabels = chains.generateChains(spins,L)
print("Done. There are {} chains in total.".format(len(tensorLabels)))

compCoef = lambda chain, state : chains.chainToCoef(chain,state,M,clebschDict)
symmMat = ImmutableSparseMatrix(0,len(symmStates),{})
tls = tensorLabels.copy()
while symmMat.rows < numSymm:
    print("Current size: ",symmMat.shape," and remaining tensors: ",len(tls))
    i = randrange(len(tls))
    ch = tls.pop(i)
    #ch = tls.pop(0)
    vec = Matrix([[compCoef(ch,state) for state in rawStates]])
    vec = vec * symmetrizer
    symmMat = (symmMat.col_join(vec)).rref(pivots=False)
    rk = symmMat.rank()
    symmMat = symmMat[:rk,:]
print("Orthonormalizing.")
# work in orthonormal basis
symmMat = symmMat*halfOfGram
print(symmMat.shape)
# make sure that the final states are orthonormal
chol = (symmMat * symmMat.transpose()).cholesky()
symmMat = chol**(-1) * symmMat
symmGram = symmMat * symmMat.transpose()
print(symmGram == eye(numSymm))

#print(symmGram.is_positive_definite, symmGram)
print("Checks:")
symmCas = casimirMatrix(spins,M,True)
altCas = halfOfGram*symmCas*(halfOfGram**(-1))
print("Is the new Casimir symmetric? ", altCas == altCas.transpose())
shiftedCas = altCas - L*(L+1)*eye(len(symmStates))
out = shiftedCas * symmMat.transpose()
print(out.norm() == 0)


# print("Computing coefficients...")
# rawTensors = [[compCoef(chain,state) for state in rawStates] for chain in tensorLabels]
# rawTensors = ImmutableMatrix(rawTensors)
# print("Done.")

# print("Symmetrizing...")
# symmTensors = (rawTensors * symmetrizer.transpose())
# symmTensors = symmTensors.rref(pivots=False)
# stRank = symmTensors.rank()
# symmTensors = symmTensors[:stRank,:]
# print("Done.")
# print("Is the number of tensors correct?", rawTensors.rows == rawDec.get(L,0), symmTensors.rows == symmDec.get(L,0))
# altST = symmTensors*halfOfGram


# symmGram = symmTensors * gramMatrix * symmTensors.transpose()
# print(symmGram)
# print(altST * altST.transpose() == symmGram)

# print("Final checks:")
# print("Computing Casimir matrices...")
# rawCas = casimirMatrix(spins,M,False)
# 

# shiftedCas = rawCas - L*(L+1)*ImmutableSparseMatrix(eye(len(rawStates)))
# diff = shiftedCas * rawTensors.transpose()
# print(diff)
# shiftedCas = symmCas - L*(L+1)*ImmutableSparseMatrix(eye(len(symmStates)))
# diff = shiftedCas * symmTensors.transpose()
# print(diff)
# print("Done.")

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
        
#         if len(tensorLabels) != rawDec.get(L,0):
#             print("Mismatch.")
#         #print("Predicting {} spin-{} multiplets; measuring {}.".format(len(tensorLabels),L,rawDec.get(L,0)))
#         coefTable = [[compCoef(chain,state) for state in rawStates] for chain in tensorLabels]
#         #mat = Matrix(coefTable)
#         #print(mat)


        
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


#     
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
