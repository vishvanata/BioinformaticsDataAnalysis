

fin = open("epex_outputs\epexchosen_BRCA.txt")
fin2 = open("textfiles/listgene.txt")
readfile = fin.readlines()
fin2.readline()
readfile2 = fin2.readlines()


# print(readfile)

npToGene = {}
for i in range(len(readfile2)):
    readfile2[i] = readfile2[i].split()
    npToGene[readfile2[i][3]] = readfile2[i][2]
# print(readfile2)
print(npToGene)
count = 0
notFound = []
for i in range(len(readfile)):
    readfile[i] = readfile[i].split()
    try:
        readfile[i][0] = npToGene[readfile[i][0]]
    except KeyError:
        count += 1
        if readfile[i][0] not in notFound:
            notFound.append(readfile[i][0])
print(count, len(readfile))
print(notFound)
