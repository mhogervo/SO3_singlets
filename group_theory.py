#########################################################################################################
### This module contains purely SO(3) group-theoretical functions.
#########################################################################################################

def triangleInEq(l1,l2,l3):
    ''' Test if {l1,l2,l3} obey the triangle inequality. '''
    
    testSpin(l1,l2,l3)
    return abs(l1-l2) <= l3 <= l1+l2

def testSpin(*l_list):
    '''
    Given an input {l_0,l_1,...}, test if all l_i are positive integers.
    '''

    for l in l_list:
        if not isinstance(l,int):
            raise TypeError("{} is not an integer.".format(l))
        elif l < 0:
            raise ValueError("{} is a negative integer.".format(l))
        else:
            pass

def testSpinLabel(l,m):
    '''
    Given quantum numbers (l,m), see if they are appropriate SO(3) labels.
    '''
    
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

def decomposeTensorProduct(l_tup):
    '''
    Compute the tensor product decomposition [l_1] x [l_2] x .... as a sum [n_0, n_1, ...] where n_j is the number of copies of [j].
    '''
    
    for l in l_tup: testSpin(l)

    maxSpin = sum(l_tup)
    empty = [0 for _ in range(0,maxSpin+1)]
    
    tp = empty.copy()
    tp[0] = 1

    if isinstance(l_tup,list):
        l_list = l_tup.copy()
    elif isinstance(l_tup,tuple):
        l_list = list(l_tup)
    else:
        raise TypeError
    
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

def countDimensionProduct(l_List):
    '''
    Given a tensor product V = [l_0] x [l_1] x ..., return dim(V).
    ''' 

    for l in l_List: testSpin(l)

    dimRep = lambda l : 2*l + 1
    
    tp = 1
    for l in l_List: tp *= dimRep(l)
    return tp

def checkMultiplicities(l_List):
    '''
    Given a tensor product V = [l_0] x [l_1] x ..., compute dim(V) in two ways.
    Either directly, or by decomposing into irreps and taking the sum of 
    dimensions.
    '''

    dec = decomposeTensorProduct(l_List)
    dim1 = countDimensionSum(dec)
    dim2 = countDimensionProduct(l_List)
    
    if dim1 == dim2:
        print("The tensor product has dim(V) = {}.".format(dim1))
        return True
    else:
        return False
