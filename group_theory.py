#########################################################################################################
### This module contains purely SO(3) group-theoretical functions.
#########################################################################################################

from sympy.physics.quantum.cg import CG as symCG

def triangleInEq(l1,l2,l3):
    '''Test if {l1,l2,l3} obey the triangle inequality. '''
    
    testSpin(l1,l2,l3)
    return abs(l1-l2) <= l3 <= l1+l2

def testSpin(*l_list):
    '''Given an input {l_0,l_1,...}, test if all l_i are positive integers.'''

    for l in l_list:
        if not isinstance(l,int):
            raise TypeError("{} is not an integer.".format(l))
        elif l < 0:
            raise ValueError("{} is a negative integer.".format(l))
        else:
            pass

def testSpinLabel(l,m):
    '''Given quantum numbers (l,m), see if they are appropriate SO(3) labels.'''
    
    testSpin(l)
    if not isinstance(m,int):
        raise TypeError("{} is not an integer.".format(m))
    elif abs(m) > l:
        raise ValueError("m = {0} is not an integer in the range [-{1},{1}]".format(m,l))
    else:
        pass
 
def tensorProduct(l1,l2,symm=False):
    '''
    Returns a list of spins that appear in the tensor product [l_1] x [l_2].
    The flag symm=True returns the symmetrized tensor product.
    '''
    testSpin(l1,l2)
    if symm and l1==l2:
        return range(0,2*l1+1,2)
    else:
        return range(abs(l1-l2),l1+l2+1)

def decomposeTensorProduct(spinTup):
    '''
    Input: a tuple of spins (l1,...,ln).
    Computes the tensor product decomposition [l_1] x [l_2] x .... as a sum [n_0, n_1, ...] where n_j is the number of copies of [j].
    Only works for non-symmetrized tensor products.
    '''
    for l in spinTup: testSpin(l)

    maxSpin = sum(spinTup)
    empty = [0 for _ in range(0,maxSpin+1)]
    
    tp = empty.copy()
    tp[0] = 1

    l_list = list(spinTup)
        
    while len(l_list) > 0:
        l_new = l_list.pop(0)
        nw = empty.copy()
        for l_old in range(0,maxSpin+1):
            n_old = tp[l_old]
            if n_old > 0:
                for j in tensorProduct(l_old,l_new): nw[j] += n_old
        tp = nw
    return tp
    
def countDimensionSum(n_List):
    '''
    Input: a vector [n_0, n_1, ...]. Let V be the representation with n_0 copies of [0], n_1 copies of [1] etc. Returns the dimension dim(V).
    '''
    
    for n in n_List:
        if not isinstance(n,int) or n < 0:
            raise ValueError("{} is not a positive integer, thus not a multiplicty.".format(n))

    dimRep = lambda l : 2*l + 1
    l_max = len(n_List) - 1
    
    return sum([n_List[l]*dimRep(l) for l in range(0,l_max+1)])

def countDimensionProduct(spinTup):
    '''
    Given a tensor product V = [l_0] x [l_1] x ..., return dim(V).
    ''' 

    for l in spinTup: testSpin(l)

    dimRep = lambda l : 2*l + 1
    
    tp = 1
    for l in spinTup: tp *= dimRep(l)
    return tp

def checkMultiplicities(spinTup):
    '''
    Given a tensor product V = [l_0] x [l_1] x ..., compute dim(V) in two ways.
    Either directly, or by decomposing into irreps and taking the sum of 
    dimensions.
    Works only for the unsymmetrized tensor product.
    '''

    dec = decomposeTensorProduct(spinTup)
    dim1 = countDimensionSum(dec)
    
    dim2 = countDimensionProduct(spinTup)
    
    if dim1 == dim2:
        print("The tensor product has dim(V) = {}.".format(dim1))
        return True
    else:
        print("Mismatch!")
        return False

    
def clebsch(s1,s2,s3, memoDict = {}):
    '''
    Return the Clebsch-Gordan symbol
    < l_1 m_1; l_2 m_2 | l_3, m_3 >.
    The entries are two-tuples: si = (li,mi).
    Can use a dictionary memoDict to store results, since this is costly.
    '''
    try:
        return memoDict[(s1,s2,s3)]
    except KeyError:
        l1,m1 = s1
        l2,m2 = s2
        l3,m3 = s3
        
        # sanity check the input
        testSpin(l1,l2,l3)
        for m in [m1,m2,m3]:
            if not isinstance(m,int): raise TypeError("{} is not an integer.".format(m))

        # check if |m| > l for any of the columns
        for x in [s1,s2,s3]:
            if abs(x[1]) > x[0]:  return 0
        
        if m1+m2 == m3 and triangleInEq(l1,l2,l3):
            # only non-trivial case
            tp = symCG(l1,m1,l2,m2,l3,m3).doit()
            memoDict[(s1,s2,s3)] = tp
            return tp
        else:
        # invalid CG symbol
            return 0
        
