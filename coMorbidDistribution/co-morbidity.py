'''
Goal:
    Correlation: co-morbidity ranking and d-score?
    Write to file and put in folder: (comorbidityVsDscore_perDisease)
    Protein | co-morbidity ranking | D-score
'''
import driverScoreHongyiList as dshl
import matplotlib.pyplot as plt

def plotCovsDscore(disease):
    diseaseName = disease
    driverScoreTuple = dshl.calculateDriverScores(diseaseName)
    proteinDscoreDict = driverScoreTuple[0]
    proteinCoPercentDict = driverScoreTuple[1]

    x = []
    y = []

    fOut = open("coMorbidDistribution/comorbidityVsDscore_perDisease/{}_cVsDs".format(diseaseName), "w")
    for protein in proteinDscoreDict:
        fOut.write("{} {} {}\n".format(protein, proteinCoPercentDict[protein], proteinDscoreDict[protein]))
        x.append(round(float(proteinCoPercentDict[protein]),2))
        y.append(proteinDscoreDict[protein])
    fOut.close()

    ''' Plot Correlation: co-morbidity ranking and d-score '''


    plt.scatter(x, y, c="orange")
    plt.title("% Comorbid diseases with protein vs. Driver-score: {}".format(diseaseName))
    plt.xlabel("Fraction of comorbid diseases")
    plt.ylabel("Driver-score")
    plt.savefig("coMorbidDistribution/comorbidityVsDscore_perDisease/{}_cVsDsScatter".format(diseaseName))
    return True
nameTonumber = {"BRCA": 2022, "SKCM": 3297, "LAML": 3340, "HNSC": 2058, "COAD": 2168, "ACC": 2640, "THCA": 2646, "ESCA": 1073, "MESO": 2292}
for key in nameTonumber.keys():
    plotCovsDscore(key)
