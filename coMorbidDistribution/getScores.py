'''
Given all proteins in {}_freq.txt goes through {} disease file in TCGA and
creates dict as follows:
    calculate the percent it is DA in TCGA
    represented as: {proteinName : (%, # DA, # normal), ...}
'''

from collections import Counter

def getScores(number, disease):

    fileNameTCGA = "epex_outputs/epexchosen_{}.txt".format(disease)
    finTCGA = open(fileNameTCGA)
    readfileTCGA = finTCGA.readlines()
    tcgaProteins = []
    for i in range(len(readfileTCGA)):
        readfileTCGA[i] = readfileTCGA[i].split()
        if readfileTCGA[i][0] not in tcgaProteins:
            tcgaProteins.append(readfileTCGA[i][0])
    tcgaProteinCount = Counter(tcgaProteins)

    tcgaDict = {}
    for list in readfileTCGA:
        if float(list[2]) > cutoff:
            if list[0] not in tcgaDict:
                tcgaDict[list[0]] = 1
            else:
                tcgaDict[list[0]] += 1

    returnDict = {}
    for protein in tcgaDict:
        tcgaDict[protein] = (round(((tcgaDict[protein]/tcgaProteinCount[protein]) * 100),1), tcgaDict[protein], tcgaProteinCount[protein])
        if protein in proteinList:
            returnDict[protein] = tcgaDict[protein]

    return returnDict, comorbidityScoreDict
