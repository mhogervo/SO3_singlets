from sympy.physics.quantum.cg import CG as symCG

def testSpinList(spinAr):
    for l in spinAr:
        if not isinstance(l,int) or l < 0:
            raise ValueError("{} is not a positive integer.".format(l))

def testSpin(m,l):
    if not isinstance(l,int) or l < 0:
        raise ValueError("{} is not a positive integer.".format(l))
    if not isinstance(m,int) or abs(m) > l:
        raise ValueError("m = {0} is not an integer in the range [-{1},{1}]".format(m,l))

def Jz(state):
    '''Given a state [(l1,m1),...,(ln,mn)], return the value of J3 = m1 + ... + mn.'''
    return sum([s[1] for s in state])

###############################
## Generating explicit bases
###############################

def fullBasisList(spinTup,symm=False):
    '''Given a list of spins (l1,...,ln), return all states in the tensor product.'''
    testSpinList(spinTup)
    
    if len(spinTup) == 0:
        return [[]]
    else:
        rep = lambda l: [(l,m) for m in range(-l,l+1)]

        spinList = list(spinTup)
        l = spinList.pop(0)
        out = [ [s] for s in rep(l) ]
        
        while len(spinList) > 0:
            l = spinList.pop(0)
            outNw = [ state + [x] for state in out for x in rep(l) ]
            out = outNw

        out = [ tuple(s) for s in out ]
        
        if symm:
            out = { tuple(sorted(s)) for s in out }
            
        return iter(out)

def spinBasisList(spinTup,L=0,symm=False):
    '''Iterator that returns basis elements of spin J3  = L.'''
    mag = lambda s : Jz(s) == L
    for s in fullBasisList(spinTup,symm):
        if mag(s): yield s
      
# def symmSpinBasisList(spinTup,L=0):
#     '''Iterator that returns _symmetrized_ basis elements of spin J3 = L.'''
#     out = { sorted(s) for s in spinBasisList(spinTup,L) }
#     for s in out: yield s

###############################
## Group-theoretical functions
###############################

def triangleInEq(spinAr):
    '''Tests if three spins [l1, l2, l3] obey the triangle inequality.'''
    if len(spinAr) != 3:
        raise ValueError("Input {} is not an array of three positive integers.".format(spinAr))
    else:
        testSpinList(spinAr)           
        l1,l2,l3 = spinAr
        return abs(l1-l2) <= l3 <= l1+l2
    
def decomposeTensorProduct(spinTup,symm=False):
    '''
    Given a tuple of spins (l1,l2,...), compute the multiplicities of spin-L
    irreps in the tensor product. 
    Output is a list [n0,n1,n2,...] where n_j is the multiplicity of spin j.
    '''
    l_max = sum(spinTup)
    mVec = [0 for l in range(0,l_max+1)]
    for state in fullBasisList(spinTup,symm):
        m = Jz(state)
        if m >= 0: mVec[m] += 1

    multVec = []
    while len(mVec) > 0:
        n = mVec.pop()
        multVec = [n] + multVec
        mVec = [j-n for j in mVec]
    return multVec

def checkMultiplicities(spinTup,multVec,symm=False):
    '''
    Given a tuple of spins (l1,...,ln) and a number of multiplicities
    (e.g. computed using decomposeTensorProduct), check that the
    dimensions agree. This only works for non-symmetrized bases.
    '''
    dim_rep = lambda l : 2*l+1
    
    if symm:
        dim_1 = sum([1 for s in fullBasisList(spinTup,symm)])
    else:
        dim_1 = 1
        for l in spinTup: dim_1 *= dim_rep(l)
    
    dim_2 = 0
    for l in range(0,len(multVec)):
        dim_2 += multVec[l]*dim_rep(l)
        
    return dim_1 == dim_2

def tensorProduct(l1,l2):
    '''Return the spins that appear in the tensor product [l1] x [l2].'''
    min = abs(l1-l2)
    max = l1+l2
    return range(min,max+1)

###############################

def spinListToTuples(spinTup,L=0):
    '''Given a list of spins (l1,...,ln), find all intermediate spins {j1,...,j_{n-2}}.'''
    spinList, N = list(spinTup), len(spinTup)
    testSpinList(spinList + [L])

    if N == 0:
        if L==0: return [[]]
        else: return []
    elif N == 1:
        l = spinList[0]
        if l == L: return [[l]]
        else: return []
    elif N == 2:
        l1,l2 = spinList
        if triangleInEq([l1,l2,L]): return [[l1,l2,L]]
        else: return []
    else: #N >= 3
        
        # leftmost CG symbol:
        l1,l2 = spinList[:2]
        out = [[l1,l2,l] for l in tensorProduct(l1,l2)]
        spinList = spinList[2:]
        
        # intermediate CG symbols:
        while len(spinList) >= 2:
            print(spinList)
            outNew = []
            m = spinList.pop(0)
            for tup in out:
                l = tup[-1]
                for r in tensorProduct(l,m): outNew.append(tup + [m,r])
            out = outNew

        # last CG symbol:
        outNew = []
        m = spinList.pop()
        for tup in out:
            l = tup[-1]
            if triangleInEq([l,m,L]): outNew.append(tup + [m,L])
        out = outNew
    
        return out

def clebsch(l1,m1,l2,m2,l3,m3):
    testSpin(m1,l1), testSpin(m2,l2), testSpin(m3,l3)
    if m1+m2 != m3 or not triangleInEq([l1,l2,l3]):
        return 0
    else:
        return symCG(l1,m1,l2,m2,l3,m3).doit()
    
def chainToCoef(spinChain,mList,M=0):
    '''
    Given a chain of spins, a set of m_i and a total magnetic quantum number M,
    compute the numerical value of the tensor T^J.
    '''
    tup = list(spinChain)
    L = tup[-1]
    testSpin(M,L)
    N = (len(tup)+1)/2
    if len(mList) != M:
        raise ValueError("Numer of m_i = {} does not match the length of the tuple, which implies N = {}.".format(len(mList),N))
    if sum(mList) != M: return 0

    else:
        # compute the leftmost CG coefficient
        l1, l2, j_r = tup[:3]
        del tup[:3]
        m1, m2 = mList.pop(0), mList.pop(0)
        mu_r = m1+m2
        testSpin(m1,l1), testSpin(m2,l2), testSpin(mu_r,j_r)
        #print([l1,m1,l2,m2,j_r,mu_r])
        tp = clebsch(l1,m1,l2,m2,j_r,mu_r)
        #print("first: ",tp)
        
        while len(tup) >= 4:
            # fill in all the middle CG coefficients
            j_l, mu_l = j_r, mu_r
            l, j_r = tup.pop(0), tup.pop(0)
            m = mList.pop(0)
            mu_r = mu_l + m
            testSpin(m,l), testSpin(mu_r,j_r)
            tp *= clebsch(j_l,mu_l,l,m,j_r,mu_r)
            #print("intermediate: ",tp)
            # the rightmost CG coefficient
        if len(tup) == 2:
            j_l, mu_l = j_r, mu_r
            l, j_r = tup[-2:] # here j_r = L
            m = mList.pop(0) # should be empty now
            testSpin(m,l)
            tp *= clebsch(j_l,mu_l,l,m,L,M)
            #print("final:",tp)
            
    #print("remaining: ",mList,tup)
    return tp
