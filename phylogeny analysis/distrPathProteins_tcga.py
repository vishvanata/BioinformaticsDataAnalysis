import distrPathProteins
import matplotlib.pyplot as plt


''' For TCGA'''
for i in range(1,10):
    n = i
    for x in range(5,10):
        x = round(x/10,2)
        cutoff = x
        fromCancerList = distrPathProteins.fromCancerList(n,x)

        casedict = {'ACC':92, 'BRCA':1044, 'ESCA':184, 'HNSC': 510 , 'LAML': 149 , 'MESO': 83 , 'SKCM': 470 , 'THCA': 496 , 'COAD': 433}
        diseaseList = casedict.keys()
        for disease in diseaseList:
            filename = "epex_outputs\epexchosen_{}.txt".format(disease)
            fin = open(filename)
            fin.readline()
            readfile = fin.readlines()
            for i in range(len(readfile)):
                readfile[i] = readfile[i].split()


        ''' create dictionary: {protein: [score, score, score], protein2: [score...]}'''
        # print(readfile)
        mydict = {}
        for list in readfile:
            if list[0] not in mydict:
                mydict[list[0]] = [round(float(list[2]),2)]
            else:
                mydict[list[0]].append(round(float(list[2]),2))

        ''' use dictionary '''
        tcgaScorelist = []
        keyErrorCount = 0

        for protein in fromCancerList:
            try:
                tcgaScorelist.extend(mydict[protein])
            except KeyError:
                keyErrorCount += 1
                continue

        ''' Plot TCGA score list '''
        scoreList = tcgaScorelist
        bins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
        plt.hist(scoreList, bins, histtype='stepfilled', rwidth=0.8, color='red')
        plt.xlabel('Pathogenicity')
        plt.ylabel('Frequency')
        # plt.title('Genes common to {} cancers'.format(input))
        plt.title('TCGA distribution | n={} | DA score > {}'.format(n, x))
        plt.savefig('compare_ExAC_distributions/{}/tcga_n{}.png'.format(x,n))
        # plt.show()
