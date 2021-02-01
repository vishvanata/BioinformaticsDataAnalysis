'''

1)  Goal is to obtain DA distribution for proteins with at least 1 DA variant
    in the score distribution of n diseases.
2)  Plot DA score distribution of just those proteins from the ExAC population.

Ideal Result:
    TCGA will be left-skewed while ExAC will be right-skewed for same proteins.
Implication:
    To what extent can filtered proteins be candidate drivers.

'''

'''
1)  Goal is to obtain DA distribution for proteins with at least 1 DA variant
    in the score distribution of n diseases.
'''
import numlist
def fromCancerList(n, value):
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
    return fromCancerList
