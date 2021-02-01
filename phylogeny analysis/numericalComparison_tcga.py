
'''
For TCGA:
1) Get list of proteins with at least 1 variant of DA score > 0.5
2) For each protein, calculate the percent it is DA in all TCGA files
   represented as: {proteinName : (%, # DA, # normal), ...}
'''

'''
1) Get list of proteins with at least 1 variant of DA score > 0.5
    Notes: Start with default n = 5.
'''
import numlist
from collections import Counter
def numCompareTcga(nin,value):
    n = nin
    cutoff = value

    casedict = {'ACC':92, 'BRCA':1044, 'ESCA':184, 'HNSC': 510 , 'LAML': 149 , 'MESO': 83 , 'SKCM': 470 , 'THCA': 496 , 'COAD': 433}
    diseaseList = casedict.keys()
    listCandidateProteins = [] #list of all proteins in all files with DA > cutoff (not filtered for idiosync)
    for disease in diseaseList:
        filename = "epex_outputs\epexchosen_{}.txt".format(disease)
        fin = open(filename)
        readfile = fin.readlines()
        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()
            if float(readfile[i][2]) > cutoff and readfile[i][0] not in listCandidateProteins:
                # scorelist.append(round(float(readfile[i][2]),5))
                # listCandidateProteins.append((readfile[i][0],readfile[i][2]))  #[(protein, DAscore), (protein, DAscore),(protein, DAscore)...]
                listCandidateProteins.append(readfile[i][0])

    fromCancerList = []
    nProteins = numlist.cancerlist(n)
    for protein in listCandidateProteins:
        if protein in nProteins[1]:
            fromCancerList.append(protein)

    '''
    2) For each protein, calculate the percent it is DA in all TCGA files
       represented as: (%, # DA, # normal)
    '''
    returnDict = {}
    tcgaProteins = []
    for disease in diseaseList:
        filename = "epex_outputs\epexchosen_{}.txt".format(disease)
        fin = open(filename)
        readfile = fin.readlines()
        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()
            if readfile[i][0] in fromCancerList and float(readfile[i][2]) > cutoff:
                if readfile[i][0] not in returnDict:
                    returnDict[readfile[i][0]] = 1
                else:
                    returnDict[readfile[i][0]] += 1
            tcgaProteins.append(readfile[i][0])
    tcgaProteinCount = Counter(tcgaProteins)
    for protein in returnDict:
        returnDict[protein] = (round(((returnDict[protein]/tcgaProteinCount[protein]) * 100),1), returnDict[protein], tcgaProteinCount[protein])
    # print(returnDict)
    return returnDict
