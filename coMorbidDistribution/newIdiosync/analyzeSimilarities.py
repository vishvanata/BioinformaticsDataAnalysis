''' Follows creation of files from analyzeDifferences.py '''

finM1 = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/9_sigSharedMed.txt")
finM2 = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/9_sigUniqueMed.txt")

readfileM1 = finM1.readlines()
readfileM2 = finM2.readlines()
m1Dict = {}
m2Dict = {}
for i in range(len(readfileM1)):
    readfileM1[i] = readfileM1[i].split()
    if readfileM1[i][1] not in m1Dict:
        m1Dict[readfileM1[i][1]] = [readfileM1[i][0]]
    else:
        m1Dict[readfileM1[i][1]].append(readfileM1[i][0])
for i in range(len(readfileM2)):
    readfileM2[i] = readfileM2[i].split()
    if readfileM2[i][1] not in m2Dict:
        m2Dict[readfileM2[i][1]] = [readfileM2[i][0]]
    else:
        m2Dict[readfileM2[i][1]].append(readfileM2[i][0])

m1Dict.update(m2Dict)
medProteinsList = []
for key in m1Dict:
    for protein in m1Dict[key]:
        medProteinsList.append(protein)

finT = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/9_sigUniqueTcga.txt")
readfileT = finT.readlines()
tDict = {}
for i in range(len(readfileT)):
    readfileT[i] = readfileT[i].split()
    if readfileT[i][1] not in tDict:
        tDict[readfileT[i][1]] = [readfileT[i][0]]
    else:
        tDict[readfileT[i][1]].append(readfileT[i][0])
# print(tDict)
tcgaProteinsList = []
for key in tDict:
    for protein in tDict[key]:
        tcgaProteinsList.append(protein)
''' Get similar proteins between two lists '''
similarCount = 0
differentCount = 0
sharedList = []
for protein in tcgaProteinsList:
    if protein in medProteinsList:
        sharedList.append(protein)
        similarCount+=1
    else:
        differentCount+=1
print("Similar: {}".format(similarCount))
print("Different: {}".format(differentCount))
print(sharedList)

proteinGroupDict = {}
for protein in sharedList:
    for key in m1Dict:
        for value in m1Dict[key]:
            if value == protein:
                medGroup = key
    for key in tDict:
        for value in tDict[key]:
            if value == protein:
                tcgaGroup = key

    proteinGroupDict[protein] = (medGroup, tcgaGroup)
# print(proteinGroupDict)
fOut = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/analyzeSimilaritiesOut.txt", "w")
for key in proteinGroupDict:
    fOut.write("{}   {}   {}\n".format(key, proteinGroupDict[key][0], proteinGroupDict[key][1] ))
fOut.close()
