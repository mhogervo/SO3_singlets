#######################################
## A module for constructing bases.
##
## A basis is always initialized with a tuple of spins (l1,...,ln).
## States are of the form ((l1,m1),...,(ln,mn)). We distinguish ordinary and symmetrized states.
##
#######################################

from group_theory import testSpin
from sympy import ImmutableSparseMatrix

def genAllData(spinTup,M=0):
    '''
    Return all information about a basis, both symmetrized and non-symmetrized.
    '''
    
    rawBasis, rawLookup = spinBasis(spinTup,M,False)
    rawDec = decomposeBasis(spinTup,False)
    
    symBasis, symLookup = spinBasis(spinTup,M,True)
    symDec = decomposeBasis(spinTup,True)
    
    mat = symmProjector(rawLookup, symLookup)
    
    return [[rawBasis, rawLookup, rawDec],[symBasis,symLookup,symDec],mat]
    
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
        
        while len(spinList) > 0:
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
    while len(mVec) > 0:
        n = mVec.pop()
        if n < 0:
            raise ValueError("Not a complete basis; one of the multiplicities became negative.")
        multVec = [n] + multVec
        mVec = [j-n for j in mVec]
        
    return { l : multVec[l] for l in range(0,l_max+1) if multVec[l] != 0 } 


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
