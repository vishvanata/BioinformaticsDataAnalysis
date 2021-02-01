

fin = open("textfiles\TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf")
for i in range(6):
    fin.readline()
readfile = fin.readlines()

for i in range(len(readfile)):
    readfile[i] = readfile[i].split()


for i in range(len(readfile)-100000):
    print(readfile[i][14])
    # for x in range(len(readfile[i])):
    #     # print(str(i) + " " + str(readfile[0][i]))
    #     print(readfile[])
# print(len(readfile))
