

hDict = {}
fIN = open("textfiles/listgene.txt")
header = fIN.readline()
readfileIn = fIN.readlines()
for i in range(len(readfileIn)):
    readfileIn[i] = readfileIn[i].split()
    hDict[readfileIn[i][1]] = readfileIn[i][3]


nameTonumber = {"BRCA": 2022, "SKCM": 3297, "LAML": 3340, "HNSC": 2058, "COAD": 2168, "ACC": 2640, "THCA": 2646, "ESCA": 1073, "MESO": 2292}
retDict = {}
for disease in nameTonumber:
    retDict[disease] = {}
    fOpen = open("textfiles/medicascy_diseasefiles/{}_ef.txt".format(nameTonumber[disease]))
    readfile = fOpen.readlines()
    proteinDict = {}
    for i in range(len(readfile)):
        readfile[i] = readfile[i].split()
        try:
            protein = hDict[readfile[i][0]]
        except:
            continue
        enrichmentFactor = readfile[i][5]
        proteinDict[protein] = enrichmentFactor
    retDict[disease] = proteinDict

fSim = open("coMorbidDistribution/newIdiosync/analyzeSimilaritiesOut.txt")
readfileSim = fSim.readlines()
for i in range(len(readfileSim)):
    readfileSim[i] = readfileSim[i].split()
    # readfileSim[i][1] = readfileSim[i][1].split(",")
print(readfileSim)

fout = open("coMorbidDistribution/newIdiosync/analyze_Similar_coUpdated.txt","w")
for list in readfileSim:
    efString = "";
    for item in list[1].split(","):
        efString = efString + retDict[item][list[0]] + ","
    efString = efString[0:-1]
    fout.write("{} | {} | {} | {}\n".format(list[0],efString,list[1],list[2]))
# print(retDict["ESCA"])
