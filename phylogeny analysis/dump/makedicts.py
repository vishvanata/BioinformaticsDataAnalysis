

# id conversion dicts straight from listgene.txt provided by Hongyi.

npgene_dict = {}
hidgene_dict = {}
hidnp_dict = {}
genenp_dict = {}
fin = open("textfiles/listgene.txt")
fin.readline()
readfin = fin.readlines()

for i in range(len(readfin)):
    readfin[i] = readfin[i].split()
    genenp_dict[readfin[i][2]] = readfin[i][3]
    npgene_dict[readfin[i][3]] = readfin[i][2]
    hidgene_dict[readfin[i][1]] = readfin[i][2]
    hidnp_dict[readfin[i][1]] = readfin[i][3]

genelist = []
for i in range(len(readfin)):
    if readfin[i][2] not in genelist:
        genelist.append(readfin[i][2])

for gene in genelist:
    count = 0
    for i in range(len(readfin)):
        if readfin[i][2] == gene:
            count += 1
    if count == 1:
        genenp_dict[gene] = readfin[i][3]
    elif count > 1:
        nplist = []
        for i in range(len(readfin)):
            if readfin[i][2] == gene:
                nplist.append(readfin[i][3])
        genenp_dict[gene] = nplist
# print(readfin)

print(genenp_dict)
# print(npgene_dict)
# print(hidgene_dict)
# print(hidnp_dict)
