

'''
Returns a dict of the form {protein : freq, protein : freq ... }
'''

import numericalComparison_tcga as nct
import numlist
from collections import Counter

def tcgaFreqDict(num, value):
    n = num
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
                listCandidateProteins.append(readfile[i][0])

    fromCancerList = []
    nProteins = numlist.cancerlist(n)
    for protein in listCandidateProteins:
        if protein in nProteins[1]:
            fromCancerList.append(protein)

    freq_dict = {}
    for disease in diseaseList:
        filename = "epex_outputs\epexchosen_{}.txt".format(disease)
        fin = open(filename)
        readfile = fin.readlines()
        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()
            if readfile[i][0] in fromCancerList:
                if readfile[i][0] not in freq_dict:
                    freq_dict[readfile[i][0]] = 1
                else:
                    freq_dict[readfile[i][0]] += 1
    return freq_dict
# print(tcgaFreqDict(5,0.5))
