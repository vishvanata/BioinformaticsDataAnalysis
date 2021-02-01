'''
Uses z-score for each protein at n and cutoff to retrieve only the proteins common to n genes
that are also significant, i.e, z-score > 1.645
'''
''' n = 5 and Cutoff = 0.5 for the initial analysis. '''
from idiosync import idiosync
import numericalComparison_zscore as ncz
import numericalComparison_relfreq as ncrf
n = 2
cutoff = 0.5
idiosyncDict = idiosync(n)
zscoreDict = ncz.zscoredict(n,cutoff)
keyErrorCount = 0
for key in idiosyncDict:
    originalList = idiosyncDict[key]
    newList = []
    for protein in originalList:
        try:
            if zscoreDict[protein] > 1.645:
                newList.append(protein)
        except KeyError:
            keyErrorCount += 1
            continue
    idiosyncDict[key] = newList

'''
5. Determine drivers (based on DA ratio) with highest frequency z-scores.
    Notes: Provides a way to rank drivers based on frequency.
    Observation: adding "and relFreqDict[protein] > 1.645" gives no proteins, so I will display each z-score next to protein.
'''
relFreqDict = ncrf.relFreqDict(n,cutoff)
relFreqDrivers = {}
for key in idiosyncDict:
    originalList = idiosyncDict[key]
    newList = []
    newListRF = []
    for protein in originalList:
        try:
            if zscoreDict[protein] > 1.645:
                newList.append(protein)
                newListRF.append((protein, relFreqDict[protein]))
        except KeyError:
            keyErrorCount += 1
            continue
    relFreqDrivers[key] = newListRF
print(idiosyncDict)
print(relFreqDrivers)
