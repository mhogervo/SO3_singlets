from bases import decomposeBasis
import pickle

cutoff = 20
fn = "tuples_{}.pickle".format(cutoff)

lRange = range(0,int(cutoff)) # modes with these spins can appear

tupEnergy = lambda tup : sum([l+1/2 for l in tup])
belowCutoff = lambda tup : tupEnergy(tup) <= cutoff
evenParity = lambda tup : sum(tup) % 2 == 0
oddParity = lambda tup : sum(tup) % 2 == 1

# generate a list of all states
allStates, newStates = [], [()]
while len(newStates) > 0:
    allStates.extend(newStates)
    prevStates = newStates
    newStates = { tuple(sorted(list(s) + [l])) for s in prevStates for l in lRange }
    newStates = list(filter(belowCutoff, newStates))
# split into even and odd parity
evenStates = sorted(list(filter(evenParity, allStates)),key=tupEnergy)
oddStates = sorted(list(filter(oddParity, allStates)),key=tupEnergy)


evenIrreps = [(tup,decomposeBasis(tup,True)) for tup in evenStates]
oddIrreps = [(tup,decomposeBasis(tup,True)) for tup in oddStates]
print(evenIrreps)
print(oddIrreps)
with open(fn,'wb') as f:
    pickle.dump((cutoff,evenIrreps,oddIrreps),f)
    
