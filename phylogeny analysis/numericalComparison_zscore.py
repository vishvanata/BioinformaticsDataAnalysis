
'''
Comparisons between TCGA and ExAC distributions for genes implicated in n diseases.
'''
'''
1. Get dictionary like the following: {protein : ratio, protein : ratio ...}
    Notes: ratio = tcga_ratio - exac_ratio. If ratio < 0, set equal to 0.
2. Make list of all of the ratios. Compute mean, std, and z-score for each protein.
3. Create dictionary like the following: {protein: z-score, protein : z-score ...}
4. Plot the z-scores, and go from there to determine proteins with highest z-scores.
'''
''' First get both dictionaries. n = 5 and cutoff = 0.5 used for initial analysis. '''
import numericalComparison_tcga as nct
import numericalComparison_exac as nce

def zscoredict(num,value):
    n = num
    cutoff = value

    tcga_dict = nct.numCompareTcga(n,cutoff)
    exac_dict = nce.numCompareExac(n,cutoff)

    ''' Note: Takes 33 seconds to fetch above 2 dictionaries '''

    ''' 1. Get dictionary like the following: {protein : ratio, protein : ratio ...} '''
    ratio_dict = {}
    keyErrored_proteins = []
    keyErrorCount = 0
    for protein in tcga_dict:
        try:
            value = round(tcga_dict[protein][0] - exac_dict[protein][0],1)
            if value < 0:
                value = 0
        except KeyError:
            value = tcga_dict[protein][0]
            keyErrored_proteins.append((protein,value))
            keyErrorCount += 1
        ratio_dict[protein] = value

    '''2. Make list of all of the ratios. Compute mean, std, and z-score for each protein.'''
    ratio_list = []
    for protein in ratio_dict:
        ratio_list.append(ratio_dict[protein])

    ''' Compute mean '''
    def Average(lst):
        return sum(lst) / len(lst)
    mean = Average(ratio_list)
    ''' StatisticsError is raised if only 1 z-score because at least 2 are required for variance.  '''
    exceptCount = 0

    ''' Compute standard deviation '''
    import statistics
    standardDev = statistics.stdev(ratio_list)

    ''' Compute z-score for each protein. Create dictionary like the following: {protein: z-score, protein : z-score ...} '''
    zscore_dict = {}
    for protein in ratio_dict:
        zscore = round((ratio_dict[protein] - mean)/standardDev, 2)
        zscore_dict[protein] = zscore
    return zscore_dict
    # print(zscore_dict['NP_003904.3'])
    # print(zscore_dict['NP_005284.2'])
    # print(zscore_dict['NP_009145.1'])
# zscoredict(5,0.5)
