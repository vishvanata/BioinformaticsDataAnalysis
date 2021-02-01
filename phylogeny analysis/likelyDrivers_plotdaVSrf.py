'''
Uses z-score for each protein at n and cutoff to retrieve only the proteins common to n genes
that are also significant, i.e, z-score > 1.645
'''
''' n = 5 and Cutoff = 0.5 for the initial analysis. '''
from idiosync import idiosync
import numericalComparison_zscore as ncz
import numericalComparison_relfreq as ncrf
for i in range(1,10):
    n = i
    for x in range(5,10):
        mynum = x/10
        cutoff = mynum
        try:
            idiosyncDict = idiosync(n)
            zscoreDict = ncz.zscoredict(n,cutoff)
            keyErrorCount = 0
            for key in idiosyncDict:
                originalList = idiosyncDict[key]
                newList = []
                for protein in originalList:
                    try:
                        if zscoreDict[protein] > 1.645:
                            newList.append(protein)
                    except KeyError:
                        keyErrorCount += 1
                        continue
                idiosyncDict[key] = newList
            # for key in idiosyncDict:
            #     for protein in idiosyncDict[key]:
            #         print(zscoreDict[protein])

            '''
            5. Determine drivers (based on DA ratio) with highest frequency z-scores.
                Notes: Provides a way to rank drivers based on frequency.
                Observation: adding "and relFreqDict[protein] > 1.645" gives no proteins, so I will display each z-score next to protein.
            '''
            relFreqDict = ncrf.relFreqDict(n,cutoff)
            relFreqDrivers = {}
            for key in idiosyncDict:
                originalList = idiosyncDict[key]
                newList = []
                newListRF = []
                for protein in originalList:
                    try:
                        if zscoreDict[protein] > 1.645:
                            newList.append(protein)
                            newListRF.append((protein, relFreqDict[protein]))
                    except KeyError:
                        keyErrorCount += 1
                        continue
                relFreqDrivers[key] = newListRF
            # print(idiosyncDict)
            # print(relFreqDrivers)
            '''
            Plot DA z-score as a function of relative frequency z-score.
            '''
            import matplotlib.pyplot as plt
            zscoreplot_list = []
            rfplot_list = []
            for key in relFreqDrivers:
                for protein in relFreqDrivers[key]:
                    rfplot_list.append(protein[1])
                    zscoreplot_list.append(zscoreDict[protein[0]])
            plt.scatter(rfplot_list, zscoreplot_list, c ="orange");
            plt.title('DA Z-score as a function of RF Z-score: n = {} and Cutoff = {}'.format(n,cutoff))
            plt.xlabel('Relative frequency Z-score')
            plt.ylabel('Disease-associated Z-score')
            plt.savefig('compare_ExAC_distributions/da_vs_rf/{}/davsrf_n{}_c{}.png'.format(n,n,cutoff))
            # plt.show()
        except:
            continue
