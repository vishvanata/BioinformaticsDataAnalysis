import statistics
nameTonumber = {"BRCA": 2022, "SKCM": 3297, "LAML": 3340, "HNSC": 2058, "COAD": 2169, "ACC": 2640, "THCA": 2646, "ESCA": 1073, "MESO": 2292}
for key in nameTonumber.keys():
    diseaseName = key
    fin = open("coMorbidDistribution\comorbidityVsDscore_perDisease\{}_cVsDs".format(diseaseName))
    readfile = fin.readlines()

    def Average(lst):
        return sum(lst) / len(lst)

    dscore_list = []
    for i in range(len(readfile)):
        readfile[i] = readfile[i].split()
        dscore_list.append(round(float(readfile[i][2]),3))
    # print(dscore_list)
    ''' Compute mean '''
    mean = Average(dscore_list)
    ''' StatisticsError is raised if only 1 z-score because at least 2 are required for variance.  '''
    exceptCount = 0
    ''' Compute standard deviation '''
    standardDev = statistics.stdev(dscore_list)
    ''' Compute z-score for each protein. Create dictionary like the following: {protein: z-score, protein : z-score ...} '''
    dzscore_dict = {}
    fOut = open("coMorbidDistribution/comorbidityVsDscore_perDisease/{}_cVsDs".format(diseaseName), "w")
    for row in readfile:
        zscore = round((round(float(row[2]),3) - mean)/standardDev, 2)
        """ Write to file """
        fOut.write("{} {} {} {}\n".format(row[0], row[1], row[2], zscore))
        """ Make it returnable by putting in dictionary """
        dzscore_dict[row[0]] = zscore

# print(dzscore_dict)
