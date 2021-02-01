

# id conversion dicts straight from listgene.txt provided by Hongyi.

npgene_dict = {}
hidgene_dict = {}
hidnp_dict = {}
fin = open("textfiles/listgene.txt")
fin.readline()
readfin = fin.readlines()

for i in range(len(readfin)):
    readfin[i] = readfin[i].split()
    npgene_dict[readfin[i][3]] = readfin[i][2]
    hidgene_dict[readfin[i][1]] = readfin[i][2]
    hidnp_dict[readfin[i][1]] = readfin[i][3]
# print(npgene_dict)
# print(hidgene_dict)
# print(hidnp_dict)
