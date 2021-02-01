
# === export

npToGene = {}
fin2 = open("textfiles/listgene.txt")
fin2.readline()
readfile2 = fin2.readlines()

for i in range(len(readfile2)):
    readfile2[i] = readfile2[i].split()
    npToGene[readfile2[i][3]] = readfile2[i][2]
