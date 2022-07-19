from bases import fullBasis, normOfState
from sympy import sqrt, diag, eye, ImmutableSparseMatrix

def listToLookup(l):
    '''
    Given a list (or tuple, or even a set) l, return a lookup table in the form of a dict.
    '''
    lookup = {}
    for n, x in enumerate(l):
        lookup[x] = n
    return lookup
        
def computeOnes(N):
    '''
    This program computes irreducible representations living in the symmetrized tensor product
    of N copies of the 1 representation.
    '''
    spinTup = tuple([1 for _ in range(N)])
    allStates = tuple(fullBasis(spinTup,True))
    mRange = range(-N,N+1)

    spinStates = {}
    sqrtGram = {}

    # Split the full Hilbert space into states with magnetic quantum number Jz = M,
    # with M running from -N, ..., N. For every M, compute the
    # list of states, a lookup table and the Gram matrix.
    Jz = lambda state : sum([x[1] for x in state])
    for M in mRange:
        projState = lambda state : Jz(state) == M
        spinning = tuple(filter(projState,allStates))

        spinStates[M] = spinning, listToLookup(spinning) 

        norms = [sqrt(normOfState(st)) for st in spinning]
        sqrtGram[M] = ImmutableSparseMatrix(diag(*norms))

    # Compute the matrices J+ and J-. They are _rectangular_, not square.
    Jmin, Jplus = {}, {}
    
    Jmin[-N] = ImmutableSparseMatrix(0,1,{})
    for M in range(-N+1,N+1):            
        normBas, normLookup = spinStates[M]
        _, lowerLookup = spinStates[M-1]
        gammaMin = lambda l,m : sqrt((l+m)*(l-m+1))
        numNorm, numLower = len(normBas), len(lowerLookup)
        lowerDict = {}
        for state in normBas:
            i = normLookup[state]
            for k in range(0,len(state)):
                tp = list(state)
                l,m = tp[k]
                if m > -l:
                    tp[k] = (l,m-1)
                    tp = tuple(sorted(tp))
                    j = lowerLookup[tp]
                    lowerDict[(j,i)] = lowerDict.get((j,i),0) + gammaMin(l,m)
        Jmin[M] = ImmutableSparseMatrix(numLower,numNorm,lowerDict)
        # appropriate normalization
        sqg0, sqgmn = sqrtGram[M], sqrtGram[M-1]
        Jmin[M] = sqgmn * Jmin[M] * sqg0**(-1)

    #Jmax[N] = ImmutableSparseMatrix(0,1,{})
    for M in range(-N,N):            
        normBas, normLookup = spinStates[M]
        _, upperLookup = spinStates[M+1]
        gammaPlus = lambda l,m : sqrt((l-m)*(l+m+1))
        numNorm, numUpper = len(normBas), len(upperLookup)
        upperDict = {}
        for state in normBas:
            i = normLookup[state]
            for k in range(0,len(state)):
                tp = list(state)
                l,m = tp[k]
                if m < l:
                    tp[k] = (l,m+1)
                    tp = tuple(sorted(tp))
                    j = upperLookup[tp]
                    upperDict[(j,i)] = upperDict.get((j,i),0) + gammaPlus(l,m)
        Jplus[M] = ImmutableSparseMatrix(numUpper,numNorm,upperDict)
        # appropriate normalization
        sqg0, sqgpl = sqrtGram[M], sqrtGram[M+1]
        Jplus[M] = sqgpl * Jplus[M] * sqg0**(-1)

    # the spins that appear in the tensor product:
    irrepRange = tuple(range(N,-1,-2).__reversed__())
    #print(irrepRange)
    irr = {}
    for L in irrepRange:
        # generate all irreps starting from the lowest-weight state, moving up using J+
        M = -L
        # there will be a unique state with Jz = -L annihilated by J-,
        # and it's the bottom rung of the spin-L multiplet.
        nulls = Jmin[M].nullspace()
        if len(nulls) != 1:
            print("Kernel of J- has the wrong dimension.")
        tp = nulls[0]
        tp /= tp.norm()
        irr[(L,M)] = tp
        # generate the other states in the multiplet
        while M < L:
            raising = Jplus[M]
            M += 1
            tp = raising * tp
            tp /= tp.norm()
            irr[(L,M)] = tp

    return spinStates, irrepRange, irr


############## main file

# max = 20
# outAr = {}
# for f in range(0,max+1):
#     print("f = {}".format(f))
#     outAr[f] = computeOnes(f)




