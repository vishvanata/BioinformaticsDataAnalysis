
# for each disease, how many variants at each pathogenicity cutoff?


import matplotlib.pyplot as plt


filename = "tcga_TCGA.BRCA.txt_pred.out"
fin = open(filename) #only proteins associated with BRCA
cancer = filename[-17:-13]

read_tcga = fin.readlines()
scorelist = []
scorelistt = [0.45,0.56,0.57,0.58,0.65,0.78]
for i in range(len(read_tcga)):
    read_tcga[i] = read_tcga[i].split()

for i in range(1,len(read_tcga),2):
    scorelist.append(read_tcga[i][0])

print(len(scorelist))
bins=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
plt.hist(scorelistt, bins, histtype='stepfilled', rwidth=0.8, color='pink')
plt.xlabel('Missense-Mutation Pathogenicity')
plt.ylabel('Frequency')
plt.title('Distribution of Pathogenicity Scores: {}'.format(cancer))

plt.show()
