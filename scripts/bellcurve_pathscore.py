
# for each disease, how many variants at each pathogenicity cutoff?


import matplotlib.pyplot as plt
import seaborn as sns
from makedicts import npgene_dict, hidnp_dict, hidgene_dict
sns.set_style("white")

medcutoff = 15

filename = "ep_outputs/tcga_TCGA.SKCM.txt_pred.out"
fin = open(filename) #only proteins associated with BRCA
medfile = open("medicascy_diseasefiles/2022_ef.txt")
readmed = medfile.readlines()
cancer = filename[-17:-13]
#medicascy list created
proteinlist = []
for i in range(len(readmed)):
    readmed[i] = readmed[i].split()
    if float(readmed[i][5]) >= medcutoff:
        try:
            proteinlist.append(hidnp_dict[readmed[i][0]])
        except:
            pass
# print(proteinlist)
read_tcga = fin.readlines()
scorelist = []

for i in range(len(read_tcga)):
    read_tcga[i] = read_tcga[i].split()
# print(read_tcga)
for i in range(1,len(read_tcga),2):
    if read_tcga[i-1][0] in proteinlist:
        scorelist.append(round(float(read_tcga[i][0]),3))
scorelist_subset = scorelist[0:1000]
# print(len(scorelist))
bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
plt.hist(scorelist, bins, histtype='stepfilled', color='pink')
plt.xlabel('Missense Mutation Pathogenicity')
plt.ylabel('Frequency')
plt.title('Distribution of Pathogenicity Scores: {}'.format(cancer))
# plt.savefig('bellcurves/epbellcurve_{}.png'.format(cancer))
plt.show()
