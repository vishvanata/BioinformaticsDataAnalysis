'''
Returns the diseases sharing common genes as key, followed by what those genes are.
$ print(idiosync(2))
$ {('ACC', 'BRCA'): ['GTF2A1L', 'RP11-438E8.2', 'RP11-410N8.4'], ('ACC', 'COAD'): ['MYADML2', 'RP11-38L1

'''

from countproteins import countproteins
from numlist import cancerlist

def idiosync(num):
    cancerdictlist = countproteins()

    genelistnum = cancerlist(num)

    returndict = {}
    for genel in genelistnum[1]:
        tupl = ()
        for dict in cancerdictlist:
            for key in dict.keys():
                for gene in dict[key].keys():
                    if genel == gene:
                        tupl += (key,)
        if tupl not in returndict.keys():
            returndict[tupl] = [genel]
        else:
            returndict[tupl].append(genel)
    # print(returndict)

    return returndict
