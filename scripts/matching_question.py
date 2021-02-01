

#different functions for different statistics
from collections import Counter

fin = open("output_files/tcga_TCGA.SKCM.txt")
read_lines = fin.readlines()
for i in range(len(read_lines)):
    read_lines[i] = read_lines[i].split()

"""
common genes and frequencies of variations for each
"""

def gene_frequencies(indvs):
    genes = []
    genes_in_all = {}


    for i in range(len(read_lines)):
        genes.append(read_lines[i][0])
    genedict = Counter(genes)
    for gene in genedict:
        if genedict[gene] >= indvs:
            genes_in_all[gene] = genedict[gene]
    print(genes_in_all)


"""
Not sure why i made this - matching genes
"""
# def matching_genes():
#
#     gene_match_count = 0
#     gene_match_list = []
#     tuplist = []   #to prevent repeats
#     #next i want to store each common gene and the frequency of matches as a tuple.
#
#     for i in range(len(read_lines)):
#         try:
#             for x in range(len(read_lines)):
#                 if read_lines[i][0] == read_lines[x][0] and read_lines.index(read_lines[x]) != i:
#                     if (read_lines.index(read_lines[x]),i) not in tuplist and (i,read_lines.index(read_lines[x])) not in tuplist:
#                         tuplist.append((read_lines.index(read_lines[x]),i))
#                         gene_match_count += 1
#                         if read_lines[i][0] not in gene_match_list:
#                             gene_match_list.append(read_lines[i][0])
#         except:
#             continue
#     print(gene_match_count)
#     print(len(gene_match_list))





"""
in a common gene, how many variations are at the same position?
"""

def samegene_sameposition():

    tuplist = []
    same_gene_pos_list = []
    same_gene_pos_count = 0

    for i in range(len(read_lines)):
        try:
            for x in range(len(read_lines)):
                if read_lines[i][0] == read_lines[x][0] and read_lines.index(read_lines[x]) != i:
                    if (read_lines.index(read_lines[x]),i) not in tuplist and (i,read_lines.index(read_lines[x])) not in tuplist:
                        tuplist.append((read_lines.index(read_lines[x]),i))
                        if int(read_lines[i][1]) == int(read_lines[x][1]) and read_lines.index(read_lines[x]) != i:
                            same_gene_pos_count += 1
                            same_gene_pos_list.append((read_lines[i],read_lines[x]))
        except:
            continue

    print(same_gene_pos_count)
    print(same_gene_pos_list)

"""
for how many genes is there a number of shared positions = individuals?
"""


# matching_genes()
samegene_sameposition()
# gene_frequencies(4)
