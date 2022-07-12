#######################################
## A module for constructing bases.
##
## A basis is always initialized with a tuple of spins (l1,...,ln).
## States are of the form ((l1,m1),...,(ln,mn)). We distinguish ordinary and symmetrized states.
##
#######################################

#from itertools import permutations
from group_theory import testSpin
from sympy import ImmutableSparseMatrix, sqrt, diag, eye, factorial

def genAllData(spinTup,M=0):
    '''
    Return all information about a basis, both symmetrized and non-symmetrized.
    '''

    print("raw")
    rawBasis, rawLookup = spinBasis(spinTup,M,False)
    #rawDec = decomposeBasis(spinTup,False)
    print("symm")
    symBasis, symLookup = spinBasis(spinTup,M,True)
    symDec = decomposeBasis(spinTup,True)
    print("norms")
    symNorms = tuple(map(normOfState2,symBasis))
    sqrtOfGram = ImmutableSparseMatrix(diag(*tuple(map(sqrt,symNorms))))
    print("proj")
    mat = symmProjector(rawLookup, symLookup)
    
    return [rawBasis, [symBasis,symLookup,symDec,sqrtOfGram],mat.transpose()]
    
def spinBasis(spinTup,M=0,symm=False):
    '''
    Return a list of states with J_z = M as well as a lookup table (in the form of a dictionary).
    The format to be used is for example
    
    basis, lookupTable = spinBasis((0,1,1,2),-2,True)
    '''
    Jz = lambda state : sum([x[1] for x in state])
    
    out = []
    for state in fullBasis(spinTup,symm):
        if Jz(state) == M: out.append(state)
    out = tuple(out)

    lookUp = {}
    for n, state in enumerate(out):
        lookUp[state] = n

    return out, lookUp

# def normOfState(state):
#     '''Compute the norm of a state.'''
#     perms, ct = permutations(state), 0
#     for x in perms:
#         if x == state: ct += 1
#     return ct

def normOfState(state):
    '''Compute the norm of a state.'''
    nd = {}
    for x in state:
        nd[x] = nd.get(x,0) + 1
    nd = [factorial(n) for n in nd.values()]
    
    out = 1
    while nd: out *= nd.pop()
    return out

        
def fullBasis(spinTup,symm=False):
    '''
    Given a tuple of spins (l1,...,ln), return all states in the tensor product.
    Returns as an iterator.
    '''
    for l in spinTup: testSpin(l)

    if len(spinTup) == 0:
        return [()]
    else:
        rep = lambda l: [(l,m) for m in range(-l,l+1)]
        spinList = list(spinTup)
        
        l = spinList.pop(0)
        out = [ [x] for x in rep(l) ]
        
        while spinList:
            l = spinList.pop(0)
            outNw = [ state + [x] for state in out for x in rep(l) ]
            if symm:
                # sort and remove duplicates: 
                outNw = { tuple(sorted(state)) for state in outNw }
                outNw = [ list(state) for state in outNw ]
            out = outNw

        out = sorted([tuple(state) for state in out])
        return iter(out)

def symmProjector(rawLookup,symmLookup):
    '''
    Given lookup tables of both unsymmetrized and symmetrized states, returns a (sympy) matrix
    that projects from unsymmetrized to symmetrized states.
    '''
    
    nRaw, nSymm = len(rawLookup), len(symmLookup)
    dd = {}
    for state in rawLookup.keys():
        i = rawLookup[state]
        symmState = tuple(sorted(state))
        j = symmLookup[symmState]
        dd[(j,i)] = 1
        
    return ImmutableSparseMatrix(nSymm,nRaw,dd) 
        
def decomposeBasis(spinTup,symm=False):
    '''
    Given a tuple of spins (l1,l2,...), compute the multiplicities of all irreps appearing in the tensor product, 
    using explicit bases. Can use both unsymmetrized and unsymmetrized tensor products.
    Output is a list [n0,n1,n2,...] where n_L is the multiplicity of spin L.
    '''
    Jz = lambda state : sum([x[1] for x in state])
    
    basis = list(fullBasis(spinTup,symm))

    # count the number of states with 
    l_max = sum(spinTup) 
    mVec = [0 for _ in range(0,l_max+1)]
    for state in basis:
        m = Jz(state)
        if m >= 0: mVec[m] += 1
        
    multVec = []
    while mVec:
        n = mVec.pop()
        if n < 0:
            raise ValueError("Not a complete basis; one of the multiplicities became negative.")
        multVec = [n] + multVec
        mVec = [j-n for j in mVec]
        
    return { l : multVec[l] for l in range(0,l_max+1) if multVec[l] != 0 } 

def casimirMatrix(spinTup,M=0,symm=False):
    '''
    For a given sector (symmetrized or not), return the action of the Casimir.
    The eigenvalues should be of the form j*(j+1) with integer j.
    '''
    
    l_max = sum(spinTup)
    if abs(M) > l_max: return ImmutableSparseMatrix(0,0,{})
        
    normBas, normLookup = spinBasis(spinTup,M,symm)
    lowerBas, lowerLookup = spinBasis(spinTup,M-1,symm)

    gammaPlus = lambda l,m : sqrt((l-m)*(l+m+1))
    gammaMin = lambda l,m : sqrt((l+m)*(l-m+1))
    
    numNorm, numLower = len(normBas), len(lowerBas)

    diagPart = M*(M-1)*ImmutableSparseMatrix(eye(numNorm))
    if M == -l_max: return diagPart
        
    lowerDict = {}
    for state in normBas:
        i = normLookup[state]
        for k in range(0,len(state)):
            tp = list(state)
            l,m = tp[k]
            if m > -l:
                tp[k] = (l,m-1)
                if symm: tp = tuple(sorted(tp))
                else: tp = tuple(tp)
                j = lowerLookup[tp]
                lowerDict[(j,i)] = lowerDict.get((j,i),0) + gammaMin(l,m)
    lowerMat = ImmutableSparseMatrix(numLower,numNorm,lowerDict)

    upperDict = {}
    for state in lowerBas:
        i = lowerLookup[state]
        for k in range(0,len(state)):
            tp = list(state)
            l,m = tp[k]
            if m < l:
                tp[k] = (l,m+1)
                if symm: tp = tuple(sorted(tp))
                else: tp = tuple(tp)
                j = normLookup[tp]
                upperDict[(j,i)] = upperDict.get((j,i),0) + gammaPlus(l,m)
    upperMat = ImmutableSparseMatrix(numNorm,numLower,upperDict)

    #print(upperMat.transpose() == lowerMat)
    
    return upperMat*lowerMat + diagPart
    
def bruteForce(casMatrix,L=0):
    '''Given a Casimir matrix, find all states with eigenvalue L(L+1).''' 
    nr = casMatrix.rows
    mat = ImmutableSparseMatrix(casMatrix - L*(L+1)*eye(nr))
    return mat.nullspace()

#################################################################################################################################
#################################################################################################################################

# def symmProjector(rawBasis,rawLookup,symmLookup):
#     '''
#     Given a basis of unsymmetrized states and the two lookup tables, return the dictionary
#     that maps an unsymmetrized state to a symmetrized one.
#     '''
    
#     nRaw, nSymm = len(rawBasis), len(symmLookup)
#     if nRaw != len(rawLookup):
#         raise ValueError("Mismatch between the number of basis elements and the lookup table.")

#     out = {}
#     for state in rawBasis:
#         i = rawLookup[state]
#         symmState = tuple(sorted(state))
#         j = symmLookup[symmState]
#         out[i] = j

#     return out
