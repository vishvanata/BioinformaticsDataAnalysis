import distrPathProteins
import matplotlib.pyplot as plt


''' For ExAC'''
for i in range(1,10):
    n=i
    for x in range(5,10):
        x = round(x/10,2)
        cutoff = x
        fromCancerList = distrPathProteins.fromCancerList(n,x)

        templist = []
        filename = "C:/Users/vishv\Desktop\Skolnick Lab\Project_2_clean/textfiles/final_entpxentp_out.txt"
        fin = open(filename)
        fin.readline()
        readfile = fin.readlines()
        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()

        ''' create dictionary: {protein: [score, score, score], protein2: [score...]}'''
        mydict = {}
        for list in readfile:
            if list[0] not in mydict:
                mydict[list[0]] = [round(float(list[2]),2)]
            else:
                mydict[list[0]].append(round(float(list[2]),2))
        # print(mydict)
        ''' use dictionary '''
        exacScorelist = []
        keyErrorCount = 0

        for protein in fromCancerList:
            try:
                exacScorelist.extend(mydict[protein])
            except KeyError:
                keyErrorCount += 1
                continue

        ''' Plot ExAC score list '''
        scoreList = exacScorelist
        bins = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
        plt.hist(scoreList, bins, histtype='stepfilled', rwidth=0.8, color='dodgerblue')
        plt.xlabel('Pathogenicity')
        plt.ylabel('Frequency')
        # plt.title('Genes common to {} cancers'.format(input))
        plt.title('ExAC distribution | n={} | DA score > {}'.format(n,x))
        plt.savefig('compare_ExAC_distributions/{}/exac_n{}.png'.format(x,n))
        # plt.show()
