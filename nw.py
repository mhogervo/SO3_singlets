from itertools import product

def fullBasisList(spinAr):
    for l in spinAr:
        if not isinstance(l,int) or l < 0: raise ValueError("{} is not a positive integer.".format(l))

    if spinAr == []:
        return [[]]
    else:
        rep = lambda l: [(l,m) for m in range(-l,l+1)]
        l = spinAr.pop(0)
        out = [ [s] for s in rep(l) ]
        #print(out)
        while len(spinAr) > 0:
            l = spinAr.pop(0)
            outNw = [ state + [x] for state in out for x in rep(l) ]
            out = outNw
        return iter(out)

def magneticSpin(state):
    return sum([s[1] for s in state])

def getBasisList(spinAr):
    mag = lambda s : magneticSpin(s) == 0
    for s in fullBasisList(spinAr):
        if mag(s): yield s

def symmBasisList(spinAr):
    out = []
    for s in getBasisList(spinAr): out += sorted(s)
    out = set(out)
    for s in out: yield s
            
def triangleInEq(spinAr):
    '''Tests if three spins [l1, l2, l3] obey the triangle inequality.'''
    if len(spinAr) != 3:
        raise ValueError("Input needs to be an array of three positive integers.")
    else:
        for l in spinAr:
            if not isinstance(l,int) or l < 0: raise ValueError("{} is not a positive integer.".format(l))
                
    l1,l2,l3 = spinAr
    if abs(l1-l2) <= l3 <= l1+l2:
        return True
    else:
        return False
    
def tensorProduct(l1,l2):
    '''Return the spins that appear in the tensor product [l1] x [l2].'''
    min = abs(l1-l2)
    max = l1+l2
    return range(min,max+1)

def genFullList(spinList):

    rep = lambda l: range(-l,l+1)

    if len(spinList) == 0:
        return [ [] ]
    elif len(spinList) == 1:
        l = spinList[0]
        if l ==0: return [ [0] ]
        else: return []
    elif len(spinList) == 2:
        l1,l2 = spinList
        if l1 == l2: return [ [l1,l1] ]
        else: return []
    elif len(spinList) == 3:
        if triangleInEq(spinAr): return [ spinList ]
        else: return []            
    else:
        l1,l2 = spinList[:2]
        out = [[l1,l2,l] for l in tensorProduct(l1,l2)]
        spinList = spinList[2:]
        print(out)
        
        while len(spinList) > 2:
            print("Remaining spins: ",spinList)
            curSpin = spinList.pop(0)
            
            outNw = []
            for tup in out:
                l1 = tup[-1]
                for l in tensorProduct(l1,curSpin):
                    outNw.append(tup + [l])
            out = outNw
            print(out)
            
            
        # finally, add the last two spins:
        out = [tup + spinList[-2:] for tup in out]

        testTriIneq = lambda ar : triangleInEq(ar[-3:])
        out = filter(testTriIneq,out)
        return out
        
