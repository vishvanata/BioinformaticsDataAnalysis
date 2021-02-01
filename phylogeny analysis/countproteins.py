
# from makedicts import npgene_dict, hidnp_dict, hidgene_dict
'''
returns
{'ACC': {'NP_060602.2': 0.2392, 'NP_000098.1': 0.2392, 'NP_004425.2': 0.2392,
relative frequency of each NPid for a cancer
'''
from collections import Counter

#only whats avaialble in exac is in npgene dict, but we dont want to exclude those
def countproteins():
    returnlist = []
    casedict = {'ACC':92, 'BRCA':1044, 'ESCA':184, 'HNSC': 510 , 'LAML': 149 , 'MESO': 83 , 'SKCM': 470 , 'THCA': 496 , 'COAD': 433}
    casetotal = 0
    for key in casedict.keys():
        casetotal += casedict[key]
    casetotal_average = casetotal/(len(casedict))

    # fout = open("output_files/countproteins/proteincount_all.txt", "w")
    for keya in casedict.keys():
        filename = "epex_outputs\epexchosen_{}.txt".format(keya)
        fin = open(filename)
        cancername = keya
        readfile = fin.readlines()

        for i in range(len(readfile)):
            readfile[i] = readfile[i].split()[0]
        # print(readfile)
        protein_freq_dict = dict(Counter(readfile))

        # [ (# of people who have cancer A) / (average number of people that have any cancer)] * protein frequency

        for gene in protein_freq_dict.keys():
            protein_freq_dict[gene] = round(protein_freq_dict[gene] * (casedict[cancername]/casetotal_average),4)
        returnlist.append({cancername : protein_freq_dict})
        # print({cancername : protein_freq_dict})

    return returnlist



# print(countproteins())
