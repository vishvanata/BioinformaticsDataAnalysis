
'''
Comparisons between TCGA and ExAC distributions for genes implicated in n diseases.
'''
'''
Iterative plots of d-score z-scores.
Driver-score is defined as:
    d-score = (tcga_ratio - exac_ratio) * tcga_frequency.
    Notes: This function takes into account frequency in determining driver status.
Steps in implementation of this function are identical to "numericalComparison_zscore.py".
'''
''' First get both dictionaries. n = 5 and cutoff = 0.5 used for initial analysis. '''
import numericalComparison_tcga as nct
import numericalComparison_exac as nce
import tcgaFreqDict as tfd
import numericalComparison_dZscore as ndzs
import matplotlib.pyplot as plt

for i in range(1,10):
    n = i
    for x in range(5,10):
        try:
            mynum = x/10
            cutoff = mynum
            dzscore_dictionary = ndzs.dzscoredict(n,cutoff)
            ''' Get score list for histogram. '''
            dzscore_list = []
            for protein in dzscore_dictionary:
                dzscore_list.append(dzscore_dictionary[protein])
            ''' Plot histogram '''
            scoreList = dzscore_list
            plt.hist(scoreList, histtype='stepfilled', rwidth=0.8, color='green')
            plt.xlabel('Z-score of Driver Status')
            plt.ylabel('Frequency')
            plt.title('Driver Status Z-score distribution: n = {} and Cutoff = {}'.format(n, cutoff))
            plt.savefig('compare_ExAC_distributions/dzscore/{}/dzscore_n{}_c{}.png'.format(n,n,cutoff))
            # plt.show()
            print('compare_ExAC_distributions/dzscore/{}/dzscore_n{}_c{}.png'.format(n,n,cutoff))
        except:
            continue
