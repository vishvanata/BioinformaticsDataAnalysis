
'''

$ print(cancerlist(6))
$ (6, ['NP_060602.2', 'NP_065915.1', 'NP_037418.3', 'NP_001128360.1', 'NP_001888.2', 'NP_000446.1',

'''

from collections import Counter
from countproteins import countproteins



def cancerlist(num):

    cancerdictlist = countproteins()

    # Build list of all genes once.
    allgenes = []
    for dict in cancerdictlist:
        for key in dict:
            for gene in dict[key]:
                    allgenes.append(gene)

    # Count occurence of each gene.
    genedict = Counter(allgenes)

    returnlist = []
    for key in genedict.keys():
        if genedict[key] == num and key not in returnlist:
            returnlist.append(key)
    return num,returnlist
