import chains
from bases import spinBasis
from sympy import Matrix

clebschDict = {}

spins = (1,2)
l_max = sum(spins)
M = -1 ### needs to obey abs(M) <= l_max
basis, lookUp = spinBasis(spins,M,False)

l_range = range(abs(M),l_max+1)
ct = 0
mat = []
for L in l_range:
    tensorLabels = chains.generateChains(spins,L)
    ct += len(tensorLabels)
    print("There are {0} states in this sector, and {1} different tensors with L = {2}.".format(len(basis), len(tensorLabels),L))
    print("The different tensors are labeled as follows: \n",tensorLabels)
    # print(basis)


    compCoef = lambda chain, state : chains.chainToCoef(chain,state,M,clebschDict)
    print("Computing matrix: ")
    coefTable = [[compCoef(chain,state) for state in basis] for chain in tensorLabels]
    mat += coefTable
    print("Done.")
    # print(len(coefTable), coefTable[0])
    #coefMat = Matrix(coefTable)
    #coefMatTranspose = coefMat.transpose()
    #print("Check that the tensors are orthonormal: ",coefMat*coefMatTranspose == Matrix.eye(len(tensorLabels)))
    # print("Computed CG coefficients: \n",clebschDict)
print("Is the number of tensors OK? ",ct == len(basis))
mat = Matrix(mat)
print("Is this a unitary change of basis? ",mat*mat.transpose() == Matrix.eye(len(basis)))
