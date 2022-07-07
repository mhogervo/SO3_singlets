import group_theory
from sympy.physics.quantum.cg import CG as symCG

# def testSpinList(spinAr):
#     for l in spinAr:
#         if not isinstance(l,int) or l < 0:
#             raise ValueError("{} is not a positive integer.".format(l))

# def testSpin(m,l):
#     if not isinstance(l,int) or l < 0:
#         raise ValueError("{} is not a positive integer.".format(l))
#     if not isinstance(m,int) or abs(m) > l:
#         raise ValueError("m = {0} is not an integer in the range [-{1},{1}]".format(m,l))

#def Jz(state):
#    '''Given a state [(l1,m1),...,(ln,mn)], return the value of J3 = m1 + ... + mn.'''
#    return sum([s[1] for s in state])

###############################
## Group-theoretical functions
###############################

# def triangleInEq(spinAr):
#     '''Tests if three spins [l1, l2, l3] obey the triangle inequality.'''
#     if len(spinAr) != 3:
#         raise ValueError("Input {} is not an array of three positive integers.".format(spinAr))
#     else:
#         testSpinList(spinAr)           
#         l1,l2,l3 = spinAr
#         return abs(l1-l2) <= l3 <= l1+l2
    
# def decomposeTensorProduct(spinTup,symm=False):
#     '''
#     Given a tuple of spins (l1,l2,...), compute the multiplicities of spin-L
#     irreps in the tensor product. 
#     Output is a list [n0,n1,n2,...] where n_j is the multiplicity of spin j.
#     '''
#     l_max = sum(spinTup)
#     mVec = [0 for l in range(0,l_max+1)]
#     for state in fullBasisList(spinTup,symm):
#         m = Jz(state)
#         if m >= 0: mVec[m] += 1

#     multVec = []
#     while len(mVec) > 0:
#         n = mVec.pop()
#         multVec = [n] + multVec
#         mVec = [j-n for j in mVec]
#     return multVec

# def checkMultiplicities(spinTup,multVec,symm=False):
#     '''
#     Given a tuple of spins (l1,...,ln) and a number of multiplicities
#     (e.g. computed using decomposeTensorProduct), check that the
#     dimensions agree. This only works for non-symmetrized bases.
#     '''
#     dim_rep = lambda l : 2*l+1
    
#     if symm:
#         dim_1 = sum([1 for s in fullBasisList(spinTup,symm)])
#     else:
#         dim_1 = 1
#         for l in spinTup: dim_1 *= dim_rep(l)
    
#     dim_2 = 0
#     for l in range(0,len(multVec)):
#         dim_2 += multVec[l]*dim_rep(l)
        
#     return dim_1 == dim_2

# def tensorProduct(l1,l2):
#     '''Return the spins that appear in the tensor product [l1] x [l2].'''
#     min = abs(l1-l2)
#     max = l1+l2
#     return range(min,max+1)

###############################

def generateChains(spinTup,L=0):
    '''
    Given a list of spins (l1,...,ln), find all intermediate spins {j1,...,j_{n-2}}.
    '''
    spinList, N = list(spinTup), len(spinTup)
    
    for l in spinTup: group_theory.testSpin(l)
    group_theory.testSpin(L)

    if N == 0:
        if L==0: return [()]
        else: return []
    elif N == 1:
        l = spinList[0]
        if l == L: return [(l,)]
        else: return []
    elif N == 2:
        l1,l2 = spinList
        if group_theory.triangleInEq(l1,l2,L): return [(l1,l2,L)]
        else: return []
    else: #N >= 3
        
        # leftmost CG symbol:
        l1,l2 = spinList[:2]
        out = [[l1,l2,l] for l in group_theory.tensorProduct(l1,l2)]
        spinList = spinList[2:]
        
        # intermediate CG symbols:
        while len(spinList) >= 2:
            #print(spinList)
            outNew = []
            mid = spinList.pop(0)
            for tup in out:
                l = tup[-1]
                for r in group_theory.tensorProduct(l,mid): outNew.append(tup + [mid,r])
            out = outNew

        # last CG symbol:
        outNew = []
        mid = spinList.pop()
        for tup in out:
            l = tup[-1]
            if group_theory.triangleInEq(l,mid,L): outNew.append(tup + [mid,L])
        out = outNew

        out = [tuple(chain) for chain in out]
        return out

# def clebsch(s1,s2,s3):
#     '''
#     Return the Clebsch-Gordan symbol
#     < l_1 m_1; l_2 m_2 | l_3, m_3 >.
#     '''
#     l1,m1 = s1
#     l2,m2 = s2
#     l3,m3 = s3
#     # sanity check the input
#     for l in [l1,l2,l3]: group_theory.testSpin(l)
#     for m in [m1,m2,m3]:
#         if not isinstance(m,int): raise TypeError("{} is not an integer.".format(m))

#     # check if |m| > l for any of the columns
#     for x in [s1,s2,s3]:
#         if abs(x[1]) > x[0]: return 0
    
#     if m1+m2 == m3 and group_theory.triangleInEq(l1,l2,l3):
#         return symCG(l1,m1,l2,m2,l3,m3).doit()
#     else:
#         # invalid CG symbol
#         return 0
        
    
def chainToCoef(spinChain,state,M=0,cd = {}):
    '''
    Given a chain of spins, a state ((l1,m1),...) i and a total magnetic quantum number M,
    compute the numerical value of the tensor T^J evaluated at this state.
    '''

    if spinChain == () and M == 0 and state == ():
        # vacuum state
        return 1
    elif len(spinChain) == 1 and len(state) == 1:
        # single-particle state | l,m > where l = L
        if state[0][1] == M:
            return 1
    elif len(spinChain) == 3 and len(state) == 2:
        # two-particle state
        L = spinChain[-1]
        if state[0][1] + state[1][1] == M:
            return group_theory.clebsch(state[0],state[1],(L,M),cd)
    
    elif len(spinChain) >= 5:
        
        mList = [x[1] for x in state]
        if sum(mList) != M: return 0
        
        L = spinChain[-1]
        group_theory.testSpinLabel(L,M)

        chainList = list(spinChain)[2:]
        fullChain = []
        fullChain.append(state[0])
        fullChain.append(state[1])
        j_r, m_r = chainList.pop(0), state[0][1] + state[1][1]
        fullChain.append((j_r,m_r))
        del mList[:2]
        
        while len(chainList) >= 4:
            j_l, m_l = j_r, m_r
            l, j_r = chainList.pop(0), chainList.pop(0)
            m = mList.pop(0)
            m_r = m_l + m
            fullChain.append((l,m))
            fullChain.append((j_r,m_r))

        if len(chainList) == 2:
            l, m = chainList.pop(0), mList.pop(0)
            fullChain.append((l,m))
            fullChain.append((L,M))
            
        # chain is complete; let's see if any m's are out of bounds
        for x in fullChain:
            if abs(x[1]) > x[0]: return 0

        # now actually compute the value:
        x_1, x_2, x_r = fullChain[:3]
        del fullChain[:3]
        tp = group_theory.clebsch(x_1,x_2,x_r,cd)
        while len(fullChain) >= 4:
            x_l = x_r
            x_m, x_r = fullChain.pop(0), fullChain.pop(0)
            tp *= group_theory.clebsch(x_l,x_m,x_r,cd)
        if len(fullChain) == 2:
            x_l = x_r
            x_m, x_r = fullChain.pop(0), fullChain.pop(0)
            tp *= group_theory.clebsch(x_l,x_m,x_r,cd)
        return tp
    
    else:
        return 0
        
        
        
        # N = (len(tup)+1)/2
        # if len(state) != N:
        #     raise ValueError("The size of the state = {} does not match the length of the tuple, which implies N = {}.".format(len(mList),N))
    
    
        # # compute the leftmost CG coefficient
        # l1, l2, j_r = tup[:3]
        # del tup[:3]
        # m1, m2 = mList.pop(0), mList.pop(0)
        # mu_r = m1+m2
        # tp = clebsch((l1,m1),(l2,m2),(j_r,mu_r))
     
        # while len(tup) >= 4:
        #     # fill in all the middle CG coefficients
        #     j_l, mu_l = j_r, mu_r
        #     l, j_r = tup.pop(0), tup.pop(0)
        #     m = mList.pop(0)
        #     mu_r = mu_l + m
        #     tp *= clebsch((j_l,mu_l),(l,m),(j_r,mu_r))
            
        # # the rightmost CG coefficient
        # if len(tup) == 2:
        #     j_l, mu_l = j_r, mu_r
        #     l, j_r = tup[-2:] # here j_r = L
        #     m = mList.pop(0) # should be empty now
            
        #     tp *= clebsch((j_l,mu_l),(l,m),(L,M))
        #     #print("final:",tp)
            
        #     #print("remaining: ",mList,tup)
        #     return tp
