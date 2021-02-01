'''
Analyze differences between 9_coProt.txt and 9_coProtTcga.txt
'''

''' Turn info in each file into dict. '''
finMed = open("C:/Users/vishv/Desktop/Skolnick Lab/Comorbid_proteins/coMorbidDistribution/newIdiosync/9_coProt.txt")
finTcga = open("C:/Users/vishv/Desktop/Skolnick Lab/Comorbid_proteins/coMorbidDistribution/newIdiosync/9_coProtTcga.txt")

readfileMed = finMed.readlines()
readfileTcga = finTcga.readlines()
for i in range(len(readfileMed)):
    readfileMed[i] = readfileMed[i].split()
for i in range(len(readfileTcga)):
    readfileTcga[i] = readfileTcga[i].split()

medDict = {}
tcgaDict = {}

for list in readfileMed:
    if list[1] not in medDict:
        medDict[list[1]] = [list[0]]
    else:
        medDict[list[1]].append(list[0])
for list in readfileTcga:
    if list[1] not in tcgaDict:
        tcgaDict[list[1]] = [list[0]]
    else:
        tcgaDict[list[1]].append(list[0])

''' How many shared keys? '''

sharedKeys = 0
for key in tcgaDict:
    if key in medDict:
        sharedKeys+=1
print("Total keys for MED: {}".format(len(medDict.keys())))
print("Total keys in TCGA: {}".format(len(tcgaDict.keys())))
print("Shared keys: {}".format(sharedKeys))

''' Create dictionary of shared keys & dictionary for unique keys '''

sharedMedDict = {}
sharedTcgaDict = {}
uniqueMedDict = {}
uniqueTcgaDict = {}
for key in medDict:
    if key in tcgaDict:
        sharedMedDict[key] = medDict[key]
        sharedTcgaDict[key] = tcgaDict[key]
    else:
        uniqueMedDict[key] = medDict[key]
for key in tcgaDict:
    if key not in medDict:
        uniqueTcgaDict[key] = tcgaDict[key]


# print(sharedMedDict)
# print(sharedTcgaDict)

''' Get number of similarity and difference between the two dicts for each key '''
similarityCount = 0
differenceCount = 0
for key in sharedMedDict:
    for value in sharedMedDict[key]:
        if value in sharedTcgaDict[key]:
            similarityCount+=1
    differenceCount = len(sharedMedDict[key]) - similarityCount
    print(key, "Similar: {}".format(similarityCount), "Different: {}".format(differenceCount), "Med Length: {}".format(len(sharedMedDict[key])), "TCGA length: {}".format(len(sharedTcgaDict[key])))

''' Repeat the previous analysis after filtration by driver z-score (Shared) '''
import getScoresAll as gs
dzscoreDict = gs.getdzscoreDict()
# f prefix stands for filtered
fSharedMedDict = {}
fSharedTcgaDict = {}
keyErrorListMed = []
keyErrorListTcga = []
# reasons for KeyErrors: (Empirically) not in TCGA cancers, not in ExAC
for key in sharedMedDict:
    for protein in sharedMedDict[key]:
        try:
            if dzscoreDict[protein] > 2:
                if key not in fSharedMedDict:
                    fSharedMedDict[key] = [protein]
                else:
                    fSharedMedDict[key].append(protein)
        except KeyError:
            keyErrorListMed.append(protein);
            continue
    for protein in sharedTcgaDict[key]:
        try:
            if dzscoreDict[protein] > 2:
                if key not in fSharedTcgaDict:
                    fSharedTcgaDict[key] = [protein]
                else:
                    fSharedTcgaDict[key].append(protein)
        except KeyError:
            keyErrorListTcga.append(protein)
            continue
''' Write the data to file (>2 driver z-score, shared, Med) '''
# fOutMed = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/9_sigSharedMed.txt", "w")
# for key in fSharedMedDict:
#     for protein in fSharedMedDict[key]:
#         fOutMed.write("{} {}\n".format(protein, key))
# fOutMed.close()
print(fSharedMedDict)
print(fSharedTcgaDict)

''' Repeat the previous analysis after filtration by driver z-score (Unique) '''
import getScoresAll as gs
dzscoreDict = gs.getdzscoreDict()
# f prefix stands for filtered
fUniqueMedDict = {}
fUniqueTcgaDict = {}
keyErrorListUniqueMed = []
keyErrorListUniqueTcga = []
# reasons for KeyErrors: (Empirically) not in TCGA cancers, not in ExAC
for key in uniqueMedDict:
    for protein in uniqueMedDict[key]:
        try:
            if dzscoreDict[protein] > 2:
                if key not in fSharedMedDict:
                    fUniqueMedDict[key] = [protein]
                else:
                    fUniqueMedDict[key].append(protein)
        except KeyError:
            keyErrorListUniqueMed.append(protein);
            continue
for key in uniqueTcgaDict:
    for protein in uniqueTcgaDict[key]:
        try:
            if dzscoreDict[protein] > 2:
                if key not in fUniqueTcgaDict:
                    fUniqueTcgaDict[key] = [protein]
                else:
                    fUniqueTcgaDict[key].append(protein)
        except KeyError:
            keyErrorListUniqueTcga.append(protein)
            continue

print(fUniqueMedDict)
print(fUniqueTcgaDict)

''' Write the data to file (>2 driver z-score, unique, Med) '''
# fOutUniqueMed = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/9_sigUniqueMed.txt", "w")
# for key in fUniqueMedDict:
#     for protein in fUniqueMedDict[key]:
#         fOutUniqueMed.write("{} {}\n".format(protein, key))
# fOutUniqueMed.close()

''' Write the data to file (>2 driver z-score, unique, Tcga) '''
# fOutUniqueTcga = open("Comorbid_proteins/coMorbidDistribution/newIdiosync/9_sigUniqueTcga.txt", "w")
print("Number of TCGA unique: {}".format(len(fUniqueTcgaDict.keys())) )
# for key in fUniqueTcgaDict:
#     for protein in fUniqueTcgaDict[key]:
#         fOutUniqueTcga.write("{} {}\n".format(protein, key))
# fOutUniqueTcga.close()
