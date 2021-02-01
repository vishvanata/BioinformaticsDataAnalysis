
'''
For ExAC:
1) Get list of proteins with at least 1 variant of DA score > 0.5
2) For each protein, calculate the percent it is DA in all TCGA files
   represented as: {proteinName : (%, # DA, # normal), ...}
'''
from collections import Counter
import distrPathProteins as dpp

def numCompareExac(nin,value):
    n = nin
    cutoff = value

    filename = "textfiles/final_entpxentp_out.txt"
    fin = open(filename)
    fin.readline()
    readfile = fin.readlines()
    for i in range(len(readfile)):
        readfile[i] = readfile[i].split()

    exacProteins = []
    for list in readfile:
        exacProteins.append(list[0])
    exacProteinCount = Counter(exacProteins)
    # print(exacProteinCount)

    exacDict = {}
    for list in readfile:
        if float(list[2]) > cutoff:
            if list[0] not in exacDict:
                exacDict[list[0]] = 1
            else:
                exacDict[list[0]] += 1
    for protein in exacDict:
        exacDict[protein] = (round(((exacDict[protein]/exacProteinCount[protein]) * 100),1), exacDict[protein], exacProteinCount[protein])
    # print(returnDict)
    ''' Get only those proteins from fromCancerList, i.e, those that occur in only n diseases'''
    returnDict = {}
    fromCancerList = dpp.fromCancerList(n,cutoff)
    keyErrorCount = 0
    for protein in fromCancerList:
        try:
            returnDict[protein] = exacDict[protein]
        except KeyError:
            keyErrorCount += 1
    return returnDict
# print(returnDict)
# print(keyErrorCount)
