from bases import decomposeBasis
from group_theory import triangleInEq

spinTup = (1,1,2,2,3)
L=2

oneTup, restTup = tuple(filter((1).__eq__,spinTup)), tuple(filter((1).__ne__,spinTup))
num_ones = len(oneTup)

oneIrreps, restIrreps = decomposeBasis(oneTup,True), decomposeBasis(restTup,True)
spinPoss = [(i,j) for i in oneIrreps.keys() for j in restIrreps.keys() if triangleInEq(i,j,L)]
print(spinPoss)
