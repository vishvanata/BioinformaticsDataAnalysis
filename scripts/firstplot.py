
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import Counter
"""
Thought: this can be done across all cancers, for an exploratory analysis
"""

#plot number above pathogenicity cutoff  as a function of  gene variation frequency
file = "ep_outputs/tcga_TCGA.MESO.txt_pred.out"
fin = open(file,"r")
cancer = file[-17:-13]

readfile = fin.readlines()
cutoff = 0.5
#long block is to fix spacing problem of entprise output
for i in range(len(readfile)):
    readfile[i] = readfile[i].split()
for i in range(len(readfile)):
    if i < (len(readfile)-1):
        for x in range(len(readfile[i+1])):
            readfile[i].append(readfile[i+1][x])
for i in range(len(readfile)):
    try:
        if type(float(readfile[i][0])):
            readfile.remove(readfile[i])
    except:
        continue

"""find how many times each gene appears, and then the amount of times that frequency appears."""
#frequency dict made here
genes = []
for item in readfile:
    genes.append(item[0])
frequencies_dict = Counter(genes)

#unique gene list
unique_genes = []
for item in readfile:
    if item not in unique_genes:
        unique_genes.append(item[0])
frequencies_dict = Counter(genes)

#above cutoff dict made here
abovecutoff_dict = {}
for gene in unique_genes:
    abovecutoff_dict[gene] = 0
    for item in readfile:
        if item[0] == gene:
            if float(item[4]) > cutoff:
                abovecutoff_dict[gene] += 1

#list of tuples
tuplist = []
for gene in frequencies_dict.keys():
    (x,y) = (frequencies_dict[gene], abovecutoff_dict[gene])
    tuplist.append((x,y))
#make scatterplot
zip(*tuplist)
plt.scatter(*zip(*tuplist), color = "green")
plt.title("{}: Gene Mutation frequency vs frequency >{} ENTPRISE cutoff".format(cancer,cutoff))
plt.xlabel("Gene Mutation frequency ")
plt.ylabel("Frequency >{} ENTPRISE cutoff".format(cutoff))
plt.show()
