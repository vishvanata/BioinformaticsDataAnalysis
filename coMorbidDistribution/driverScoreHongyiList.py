

'''
Calculates d-score dict. ExAC and one disease file.
Inputs:
    diseaseName
Output:
    driver score dict: {protein: d-score, protein: d-score ... }
'''
import getExacDict
import getScores as gs

def calculateDriverScores(diseaseOI):
    diseaseName = diseaseOI
    nameTonumber = {"BRCA": 2022, "SKCM": 3297, "LAML": 3340, "HNSC": 2058, "COAD": 2168, "ACC": 2640, "THCA": 2646, "ESCA": 1073, "MESO": 2292}

    exacDict = getExacDict.getExacDict(0.5)
    getScores = gs.getScores(nameTonumber[diseaseName], diseaseName)

    driverScoreDict = {}
    for protein in getScores[0]:
        tcgaFrequency = getScores[0][protein][1] + getScores[0][protein][2]
        '''d-score = (tcga_ratio - exac_ratio) * tcga_frequency'''
        try:
            driverScore = (getScores[0][protein][0] - exacDict[protein][0])/100 * tcgaFrequency
        except KeyError:
            driverScore = (getScores[0][protein][0])/100 * tcgaFrequency
        if driverScore >= 0:
            driverScoreDict[protein] = round(driverScore, 2)
        else:
            driverScoreDict[protein] = 0

    return driverScoreDict, getScores[1]
print(calculateDriverScores("LAML")
lamlDict = calculateDriverScores("LAML")
for protein in lamlDict[0]:
    if lamlDict[0][protein] > 5:
        print(protein, lamlDict[0][protein])
