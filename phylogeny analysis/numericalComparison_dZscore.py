
'''
Comparisons between TCGA and ExAC distributions for genes implicated in n diseases.
'''
'''
Driver-score is defined as:
    d-score = (tcga_ratio - exac_ratio) * tcga_frequency.
    Notes: This function takes into account frequency in determining driver status.
Steps in implementation of this function are identical to "numericalComparison_zscore.py".
'''
''' First get both dictionaries. n = 5 and cutoff = 0.5 used for initial analysis. '''
import numericalComparison_tcga as nct
import numericalComparison_exac as nce
import tcgaFreqDict as tfd

def dzscoredict(num,value):
    n = num
    cutoff = value

    tcga_dict = nct.numCompareTcga(n,cutoff)
    exac_dict = nce.numCompareExac(n,cutoff)
    tcgafreq_dict = tfd.tcgaFreqDict(n, cutoff)
    ''' Note: Takes 33 seconds to fetch above 2 dictionaries '''

    ''' 1. Get dictionary like the following: {protein : ratio, protein : ratio ...} '''
    dscore_dict = {}
    keyErrored_proteins = []
    keyErrorCount = 0
    for protein in tcga_dict:
        try:
            value = round(((tcga_dict[protein][0] - exac_dict[protein][0])/100)*tcgafreq_dict[protein],1)
            if value < 0:
                value = 0
        except KeyError:
            value = tcga_dict[protein][0]*tcgafreq_dict[protein]
            keyErrored_proteins.append((protein,value))
            keyErrorCount += 1
        dscore_dict[protein] = value
    # print(dscore_dict)
    '''2. Make list of all of the d-scores (see header comment for what d-score is). Compute mean, std, and z-score for each protein.'''
    dscore_list = []
    for protein in dscore_dict:
        dscore_list.append(dscore_dict[protein])
    # print(dscore_list)
    ''' Compute mean '''
    def Average(lst):
        return sum(lst) / len(lst)
    mean = Average(dscore_list)
    ''' StatisticsError is raised if only 1 z-score because at least 2 are required for variance.  '''
    exceptCount = 0

    ''' Compute standard deviation '''
    import statistics
    standardDev = statistics.stdev(dscore_list)

    ''' Compute z-score for each protein. Create dictionary like the following: {protein: z-score, protein : z-score ...} '''
    dzscore_dict = {}
    for protein in dscore_dict:
        zscore = round((dscore_dict[protein] - mean)/standardDev, 2)
        dzscore_dict[protein] = zscore
    return dzscore_dict
    # print(zscore_dict['NP_003904.3'])
    # print(zscore_dict['NP_005284.2'])
    # print(zscore_dict['NP_009145.1'])
print(dzscoredict(5,0.5))
