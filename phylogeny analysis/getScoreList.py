
'''
input dictionary
add all pathogenicity score into same list for n
'''

from idiosync import idiosync

def getScoreList(input):

    count = 0
    filenamelist = []
    scorelist = []

    for tup in input.keys():
        for gene in tup:

            filename = "epex_outputs\epexchosen_{}.txt".format(gene)
            filenamelist.append(filename)
            fin = open(filename)
            readfile = fin.readlines()
            for i in range(len(readfile)):
                readfile[i] = readfile[i].split()
                for NPid in input[tup]:
                    if readfile[i][0] == NPid:
                        scorelist.append(round(float(readfile[i][2]),5))

    return scorelist  #, tup
#     # print(idiodict)

# print(getScoreList(idiosync(7)));
