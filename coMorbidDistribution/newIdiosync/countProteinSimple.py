'''
Idiosyncratic proteins based on MEDICASCY comorbid ranked
Output: file 9_coProt.txt
'''
def getIdioMed():
    from collections import Counter

    finalProteinList = []
    cancerProteinDict = {}
    nameTonumber = {"BRCA": 2022, "SKCM": 3297, "LAML": 3340, "HNSC": 2058, "COAD": 2168, "ACC": 2640, "THCA": 2646, "ESCA": 1073, "MESO": 2292}
    for disease in nameTonumber:
        addingList = []
        fin = open("LemeDiscoFiles/{}_freq.txt".format(nameTonumber[disease]))
        readfile = fin.readlines();
        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()
            if readfile[i][5] not in addingList:
                addingList.append(readfile[i][5])
        cancerProteinDict[disease] = addingList
        finalProteinList += addingList

    sharedProteinCount = Counter(finalProteinList);
    # print(cancerProteinDict)
    # print(sharedProteinCount)

    ''' Create tuples '''
    returnDict = {}
    for protein in sharedProteinCount:
        cancersTup = ()
        for cancer in cancerProteinDict:
            if protein in cancerProteinDict[cancer]:
                cancersTup += (cancer,)
        if cancersTup not in returnDict:
            returnDict[cancersTup] = [protein]
        else:
            returnDict[cancersTup].append(protein)
    return returnDict


# ''' Write all this info to a file '''
# fOut = open("coMorbidDistribution/newIdiosync/9_coProt.txt", "w")
# for tupl in returnDict:
#     tupString = ""
#     for item in tupl:
#         tupString += (item + ",")
#     tupString = tupString[:-1] #remove end commas
#     for protein in returnDict[tupl]:
#         fOut.write("{} {}\n".format(protein, tupString))

''' Calculate driver scores for each protein in each tuple '''
