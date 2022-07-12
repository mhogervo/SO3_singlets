import group_theory

def generateChains(spinTup,L=0):
    '''
    Given a list of spins (l1,...,ln), find all intermediate spins {j1,...,j_{n-2}}.
    '''
    for l in spinTup: group_theory.testSpin(l)
    group_theory.testSpin(L)
    
    spinList, N = list(spinTup), len(spinTup)
     
    if N == 0:
        # vacuum 
        if L==0: return [()]
        else: return []
    elif N == 1:
        # single-particle state
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
            outNew = []
            mid = spinList.pop(0)
            for tup in out:
                left = tup[-1]
                for right in group_theory.tensorProduct(left,mid): outNew.append(tup + [mid,right])
            out = outNew

        # last CG symbol:
        outNew = []
        mid = spinList.pop()
        for tup in out:
            left = tup[-1]
            if group_theory.triangleInEq(left,mid,L): outNew.append(tup + [mid,L])
        out = outNew

        out = [tuple(chain) for chain in out]
        return out


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
        l,m = state[0]
        L = spinChain[0]
        if l == L and m == M and abs(M) <= L:
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

        #l1, l2 = state[0][0], state[1][0]
        #if [l1,l2] != spinChain[:2]:
            # chain does not match with the desired state
            #return 0
        chainList = list(spinChain)[2:]
        fullChain = []
        fullChain.append(state[0])
        fullChain.append(state[1])
        j_r, m_r = chainList.pop(0), state[0][1] + state[1][1]
        if abs(m_r) > j_r: return 0
        fullChain.append((j_r,m_r))
        del mList[:2]
        
        while len(chainList) >= 4:
            j_l, m_l = j_r, m_r
            l, j_r = chainList.pop(0), chainList.pop(0)
            m = mList.pop(0)
            m_r = m_l + m
            if abs(m_r) > j_r: return 0
            fullChain.append((l,m))
            fullChain.append((j_r,m_r))

        if len(chainList) == 2:
            l, m = chainList.pop(0), mList.pop(0)
            fullChain.append((l,m))
            fullChain.append((L,M))

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
        
