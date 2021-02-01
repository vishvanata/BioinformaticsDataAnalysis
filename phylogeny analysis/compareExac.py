'''

From ExAC:
{Protein: (# Above Cutoff, # Total, Percent)}
{'NP_940978.2': (48, 643, 7.5), 'NP_003318.1': (3, 89, 3.4), 'NP_060341.2': (68, 220, 30.9),

'''
from collections import Counter

cutoff = 0.5

filename = "C:/Users/vishv\Desktop\Skolnick Lab\Project_2_clean/textfiles/final_entpxentp_out.txt"
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

returnDict = {}
for list in readfile:
    if float(list[2]) > cutoff:
        if list[0] not in returnDict:
            returnDict[list[0]] = 1
        else:
            returnDict[list[0]] += 1
for protein in returnDict:
    returnDict[protein] = (returnDict[protein], exacProteinCount[protein], round(((returnDict[protein]/exacProteinCount[protein]) * 100),1) )
# print(returnDict)
