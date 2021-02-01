
# for each disease, how many variants at each pathogenicity cutoff?


import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
filename = "textfiles/entpoutu.txt"
fin = open(filename) #only proteins associated with BRCA
cancer = filename[-17:-13]

read_tcga = fin.readlines()
scorelist = []
scorelistt = [0.45,0.56,0.57,0.58,0.65,0.78]
for i in range(len(read_tcga)):
    read_tcga[i] = read_tcga[i].split()
    scorelist.append(round(float(read_tcga[i][4]),3))


print(len(scorelist))
print((scorelist[0:10]))
bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
plt.hist(scorelist, bins, histtype='stepfilled', rwidth=0.8, color='green')
plt.xlabel('Missense Mutation Pathogenicity')
plt.ylabel('Frequency')
plt.title('Distribution of Pathogenicity Scores: ExAC')
plt.savefig('bellcurves/epbellcurve_exac.png')
plt.show()
