
'''
Comparisons between TCGA and ExAC distributions for genes implicated in n diseases.
'''
'''
First get both dictionaries. n = 5 and cutoff = 0.5 used for initial analysis.
'''
import numericalComparison_tcga as nct
import numericalComparison_exac as nce
import matplotlib.pyplot as plt

for i in range(1,10):
    n = i
    for x in range(5,10):
        cutoff = x/10

        tcga_dict = nct.numCompareTcga(n,cutoff)
        exac_dict = nce.numCompareExac(n,cutoff)
        ''' Note: Takes 33 seconds to fetch above 2 dictionaries '''

        '''
        Plot the distribution of ratios at the same time and side-by-side.
        '''
        ''' First get list of scores. '''

        tcga_ratios = []
        exac_ratios = []
        for protein in tcga_dict:
            tcga_ratios.append(tcga_dict[protein][0])
        for protein in exac_dict:
            exac_ratios.append(exac_dict[protein][0])

        ''' Now actually plot them. '''

        fig, ax = plt.subplots(1,2)
        bins = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        ax[0].hist(tcga_ratios, bins, histtype='stepfilled', color = 'r')
        ax[1].hist(exac_ratios, bins, histtype='stepfilled', color = 'g')
        ax[0].set_title("TCGA")
        ax[1].set_title("ExAC")
        fig.suptitle("% Disease Associated: n = {} & Cutoff = {}".format(n,cutoff))
        # plt.show()
        plt.savefig('compare_ExAC_distributions/ratios/n_{}_c_{}.png'.format(n,cutoff))
        print('compare_ExAC_distributions/ratios/n_{}_c_{}.png'.format(n,cutoff))
