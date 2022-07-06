###############################
## Constructing bases.
## A state is of the form ((l1,m1),...,(ln,mn)). We distinguish ordinary and symmetrized states.
###############################

import group_theory

def Jz(state):
    return sum([s[1] for s in state])

def fullBasis(spinTup,symm=False):
    '''
    Given a tuple of spins (l1,...,ln), return all states in the tensor product.
    '''
    for l in spinTup: group_theory.testSpin(l)

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

def spinBasis(spinTup,L=0,symm=False):
    '''
    Return a list of states with J_z = L as well as a lookup table.
    '''
    
    out = []
    for state in fullBasis(spinTup,symm):
        if Jz(state) == L: out.append(state)
    out = tuple(out)

    lookUp = {}
    for n, state in enumerate(out):
        lookUp[state] = n

    return out, lookUp
    
def symmProjector(rawBasis,rawLookup,symmLookup):
    '''
    Given a basis of unsymmetrized states and the two lookup tables, return the dictionary
    that maps an unsymmetrized state to a symmetrized one.
    '''
    
    nRaw, nSymm = len(rawBasis), len(symmLookup)
    if nRaw != len(rawLookup):
        raise ValueError("Mismatch between the number of basis elements and the lookup table.")

    out = {}
    for state in rawBasis:
        i = rawLookup[state]
        symmState = tuple(sorted(state))
        j = symmLookup[symmState]
        out[i] = j

    return out
    
    
def decomposeBasis(spinTup,symm=False):
    '''
    Given a tuple of spins (l1,l2,...), compute the multiplicities of spin-L
    irreps in the tensor product. 
    Output is a list [n0,n1,n2,...] where n_j is the multiplicity of spin j.
    '''

    basis = list(fullBasis(spinTup,symm))

    if len(basis) == 0:
        return []
    else:
        
        l_max = sum([x[0] for x in basis[0]])     
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
            
        return { l : multVec[l] for l in range(0,l_max+1) }
