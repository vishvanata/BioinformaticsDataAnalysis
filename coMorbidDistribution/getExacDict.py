
'''
For ExAC:
1) For each protein, calculate the percent it is DA in ExAC
   represented as: {proteinName : (%, # DA, # normal), ...}
'''
from collections import Counter

def getExacDict(cutoffParam):
    cutoff = cutoffParam

    filenameEx = "textfiles/final_entpxentp_out.txt"
    finEx = open(filenameEx)
    finEx.readline()
    readfileEx = finEx.readlines()
    for i in range(len(readfileEx)):
        readfileEx[i] = readfileEx[i].split()
    exacProteins = []
    for list in readfileEx:
        exacProteins.append(list[0])
    exacProteinCount = Counter(exacProteins)

    exacDict = {}
    for list in readfileEx:
        if float(list[2]) > cutoff:
            if list[0] not in exacDict:
                exacDict[list[0]] = 1
            else:
                exacDict[list[0]] += 1

    for protein in exacDict:
        exacDict[protein] = (round(((exacDict[protein]/exacProteinCount[protein]) * 100),1), exacDict[protein], exacProteinCount[protein])

    return exacDict
