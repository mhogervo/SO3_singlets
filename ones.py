from bases import fullBasis, normOfState
from sympy import sqrt, diag, eye, ImmutableSparseMatrix

def computeOnes(factors):
    spinTup = tuple([1 for _ in range(factors)])
    allStates = tuple(fullBasis(spinTup,True))
    mRange = range(-factors,factors+1)

    spinStates = {}
    sqrtGram = {}
    
    Jz = lambda state : sum([x[1] for x in state])
    for M in mRange:
        projState = lambda state : Jz(state) == M
        spinning = tuple(filter(projState,allStates))
        lookup = {}
        for n, st in enumerate(spinning):
            lookup[st] = n
        spinStates[M] = spinning, lookup
        norms = [sqrt(normOfState(st)) for st in spinning]
        sqrtGram[M] = ImmutableSparseMatrix(diag(*norms))


    Jmin, Jplus = {}, {}
    
    nmin = len(spinStates[-factors][0])
    Jmin[-factors] = ImmutableSparseMatrix(0,nmin,{})
    for M in range(-factors+1,factors+1):            
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
        
    # don't need Jplus for the top state:
    for M in range(-factors,factors):            
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
    
    irrepRange = range(factors,-1,-2)
    #print(irrepRange)
    irr = {}
    for L in irrepRange:
        #print("Working on spin {}".format(L))
        M = -L
        #print(Jmin[-L])
        tp = Jmin[M].nullspace()[0]
        tp /= tp.norm()
        irr[(L,M)] = tp
        # generate the other states in the multiplet
        while M < L:
            raising = Jplus[M]
            M += 1
            tp = raising * tp
            tp /= tp.norm()
            irr[(L,M)] = tp

    return spinStates, list(irrepRange), irr


############## main file

max = 20
outAr = {}
for f in range(0,max+1):
    print("f = {}".format(f))
    outAr[f] = computeOnes(f)







        # print("Casimir?")
        # norms = []
        # for M in range(-L,L+1):
        #     ns = len(spinStates[M][0])
        #     shiftedCasimir = Cas[M] - L*(L+1)*eye(ns)
        #     tp = shiftedCasimir*goodStates[(L,M)]
        #     norms.append(tp.norm())
        # print(norms)
    
        
#print(goodStates)

#print([(M,Cas[M].eigenvals()) for M in mRange])


# print("Checking commutation relations:")
# for M in mRange:
#     ns = len(spinStates[M][0])
    
#     if M > -factors:
#         ft =  Jplus[M-1] * Jmin[M]
#     else:
#         ft = ImmutableSparseMatrix(ns,ns,{})
#     if M < factors:
#         st =  Jmin[M+1] * Jplus[M]
#     else:
#         st = ImmutableSparseMatrix(ns,ns,{})
#     diff = ft - st - 2*M*eye(ns)
#     if diff.norm() != 0:
#         print("Commutation relations are violated.")
# print("Done.")        

  
    # Jplus[factors] = ImmutableSparseMatrix(0,len(spinStates[factors][0]),{})   

    #change basis to make the states orthonormal
    #for M in mRange:
        
   #     if M < factors:
            
            
        # if M > -factors:
        #     
        #     

    # for M in mRange:
    #     #print(M)
    #     ns = len(spinStates[M][0])
    #     if M == -factors:
    #         cm = M*(M-1)*eye(ns)
    #     else:
    #         cm = Jplus[M-1]*Jmin[M] + M*(M-1)*eye(ns)
    #     Cas[M] = cm
