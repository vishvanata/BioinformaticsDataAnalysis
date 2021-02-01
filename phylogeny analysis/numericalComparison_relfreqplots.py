

'''
Comparisons between TCGA and ExAC distributions for genes implicated in n diseases.
'''
'''
Reasoning:
    Which proteins are drivers in an individual case can vary based on epistatic context.

1. Get dictionary like the following: {protein : frequency, protein : frequency ...}
    Notes: Make another file with weighted frequencies. Q: But why would we need to weight it?
    Notes: This is only using all the TCGA files, not ExAC.
2. Make list of all of the frequencies. Compute mean, std, and z-score for each protein.
3. Create dictionary like the following: {protein: z-score, protein : z-score ...}
4. Plot the z-scores, and go from there to determine proteins with highest z-scores.
    Notes: How is DA ratio correlated with frequency z-score?
4. Plot the z-scores, and go from there to determine proteins with highest z-scores.
5. Determine drivers (based on DA ratio) with highest frequency z-scores.
    Notes: Provides a way to rank drivers based on frequency.
'''
'''
1. Get dictionary like the following: {protein : frequency, protein : frequency ...}
'''
''' Get list with cancerlist for given n and cutoff. n = 5 and Cutoff = 0.5 will be used for initial analysis. '''
import numericalComparison_tcga as nct
import numlist
from collections import Counter
for i in range(1,10):
    n = i
    for x in range(5,10):
        mynum = x/10
        cutoff = mynum

        casedict = {'ACC':92, 'BRCA':1044, 'ESCA':184, 'HNSC': 510 , 'LAML': 149 , 'MESO': 83 , 'SKCM': 470 , 'THCA': 496 , 'COAD': 433}
        diseaseList = casedict.keys()
        listCandidateProteins = [] #list of all proteins in all files with DA > cutoff (not filtered for idiosync)
        for disease in diseaseList:
            filename = "epex_outputs\epexchosen_{}.txt".format(disease)
            fin = open(filename)
            readfile = fin.readlines()
            for i in range(len(readfile)):
                readfile[i] = readfile[i].split()
                if float(readfile[i][2]) > cutoff and readfile[i][0] not in listCandidateProteins:
                    # scorelist.append(round(float(readfile[i][2]),5))
                    # listCandidateProteins.append((readfile[i][0],readfile[i][2]))  #[(protein, DAscore), (protein, DAscore),(protein, DAscore)...]
                    listCandidateProteins.append(readfile[i][0])

        fromCancerList = []
        nProteins = numlist.cancerlist(n)
        for protein in listCandidateProteins:
            if protein in nProteins[1]:
                fromCancerList.append(protein)

        freq_dict = {}
        for disease in diseaseList:
            filename = "epex_outputs\epexchosen_{}.txt".format(disease)
            fin = open(filename)
            readfile = fin.readlines()
            for i in range(len(readfile)):
                readfile[i] = readfile[i].split()
                if readfile[i][0] in fromCancerList:
                    if readfile[i][0] not in freq_dict:
                        freq_dict[readfile[i][0]] = 1
                    else:
                        freq_dict[readfile[i][0]] += 1
        # print(freq_dict)
        '''
        2. Make list of all of the frequencies. Compute mean, std, and z-score for each protein.
        '''
        ''' Make list '''
        freq_list = []
        for protein in freq_dict:
            freq_list.append(freq_dict[protein])

        ''' Compute mean '''
        def Average(lst):
            return sum(lst) / len(lst)
        mean = Average(freq_list)

        ''' try/except: StatisticsError is raised if only 1 z-score because at least 2 are required for variance.  '''
        exceptCount = 0
        try:
            ''' Compute standard deviation '''
            import statistics
            standardDev = statistics.stdev(freq_list)

            ''' Compute z-score for each protein. Create dictionary like the following: {protein: z-score, protein : z-score ...} '''
            zscore_dict = {}
            for protein in freq_dict:
                zscore = round((freq_dict[protein] - mean)/standardDev, 2)
                zscore_dict[protein] = zscore

            ''' 4. Plot the z-scores, and go from there to determine proteins with highest z-scores. '''
            zscore_list = []
            for protein in zscore_dict:
                zscore_list.append(zscore_dict[protein])

            ''' Actually make plot '''
            import matplotlib.pyplot as plt

            scoreList = zscore_list
            plt.hist(scoreList, histtype='stepfilled', rwidth=0.8, color='orange')
            plt.xlabel('z-score')
            plt.ylabel('Frequency')
            plt.title('Z-score distribution for protein frequency: n = {} and Cutoff = {}'.format(n, cutoff))
            plt.savefig('compare_ExAC_distributions/relfreq/{}/relfreq_n{}_c{}.png'.format(n,n,cutoff))
            # plt.show()
            print('compare_ExAC_distributions/relfreq/{}/relfreq_n{}_c{}.png'.format(n,n,cutoff))
        except:
            exceptCount += 1
            continue
