'''
D-score z-score for all proteins regardless of cancer.
'''

from collections import Counter

def getScoresAll(cutoffInput):
    cutoff = cutoffInput
    ''' Create master protein : dzscore dict '''
    nameTonumber = {"BRCA": 2022, "SKCM": 3297, "LAML": 3340, "HNSC": 2058, "COAD": 2168, "ACC": 2640, "THCA": 2646, "ESCA": 1073, "MESO": 2292}
    proteinList = []
    tcgaDict = {}
    for cancer in nameTonumber:
        fin = open("C:/Users/vishv/Desktop/Skolnick Lab/Comorbid_proteins/epex_outputs/epexchosen_{}.txt".format(cancer))
        readfile = fin.readlines()
        addToList = []
        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()
            addToList.append(readfile[i][0])
            if float(readfile[i][2]) > cutoff:
                if readfile[i][0] not in tcgaDict:
                    tcgaDict[readfile[i][0]] = 1
                else:
                    tcgaDict[readfile[i][0]] += 1

        proteinList += addToList

    proteinCount = Counter(proteinList)
    returnDict = {}
    for protein in tcgaDict:
        returnDict[protein] = (round(((tcgaDict[protein]/proteinCount[protein]) * 100),1), tcgaDict[protein], proteinCount[protein])

    return returnDict
# print(getScoresAll())

'''Function to calculate list containing driver score z-scores for each protein '''
def getdzscoreDict():
    import getExacDict as ged
    exacDict = ged.getExacDict(0.5)
    tcgaDict = getScoresAll(0.5)
    returnDict = {}
    dscoreList = []
    dscoreDict = {}
    keyErrorCount = 0
    for protein in tcgaDict:
        tcgaProp = round(float(tcgaDict[protein][0])/100,2)
        try:
            dscore = (tcgaProp - (round(float(exacDict[0][protein][0])/100),2))*float(tcgaDict[protein][2])
            if dscore < 0:
                dscore = 0
            dscoreDict[protein] = dscore
        except KeyError:
            keyErrorCount += 1
            dscore = tcgaProp*float(tcgaDict[protein][2])
            if dscore < 0:
                dscore = 0
            dscoreDict[protein] = dscore
            continue
    ''' Calculate list of z-scores '''

    dscoreList = []
    for protein in dscoreDict:
        dscoreList.append(float(dscoreDict[protein]))
    ''' Compute mean '''
    def Average(lst):
        return sum(lst) / len(lst)
    mean = Average(dscoreList)

    ''' Compute standard deviation '''
    import statistics
    standardDev = statistics.stdev(dscoreList)

    ''' Compute z-score for each protein. Create dictionary like the following: {protein: z-score, protein : z-score ...} '''
    dzscoreDict = {}
    for protein in dscoreDict:
        zscore = round((dscoreDict[protein] - mean)/standardDev, 2)
        dzscoreDict[protein] = zscore
    return dzscoreDict
