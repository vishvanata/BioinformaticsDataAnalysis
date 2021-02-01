from makedicts import npgene_dict, hidnp_dict, hidgene_dict
from collections import Counter

ecutoff = 0.5
mcutoff = 1

exac = open("textfiles/entpoutu.txt")
med = open("medicascy_diseasefiles/2022_ef.txt")
tcga = open("ep_outputs/tcga_TCGA.BRCA.txt_pred.out")
read_exac = exac.readlines()
read_med = med.readlines()
read_tcga = tcga.readlines()
hiddict = {}

#MEDICASCY protein list established
proteinlist = []
for i in range(len(read_med)):
    read_med[i] = read_med[i].split()
    if float(read_med[i][5]) >= mcutoff:
        try:
            proteinlist.append(hidnp_dict[read_med[i][0]])
        except:
            pass
print("GOT HERE!")
if len(proteinlist) == 0:
    print("MEDICASCY cutoff too high!")

# rest of code is filtration based off of cutoffs
else:

    #tcga filter
    count_tcga = 0
    tcga_variantlist = []

    for i in range(len(read_tcga)):
        read_tcga[i] = read_tcga[i].split()
    for protein in proteinlist:
        for i in range(1,len(read_tcga)):
            if read_tcga[i-1][0] == protein and float(read_tcga[i][0]) > ecutoff:
                count_tcga += 1
                tcga_variantlist.append((read_tcga[i-1][0]))
    tcga_freqdict = Counter(tcga_variantlist)




    # exac filter
    count_exac = 0
    exac_variantlist = []

    for i in range(len(read_exac)):
        read_exac[i] = read_exac[i].split()
        if read_exac[i][0] in proteinlist and float(read_exac[i][4]) > ecutoff:
            count_exac += 1
            exac_variantlist.append((read_exac[i][0]))
    exac_freqdict = Counter(exac_variantlist)

    #compare frequency at which certain protein in tcga is pathogenic
    for key in tcga_freqdict:
        if tcga_freqdict[key] == 3:
            print(key,exac_freqdict[key])


    print("{} variants which meet cutoffs in ExAC".format(count_exac))
    print("{} variants which meet cutoffs in TCGA".format(count_tcga))
    print(tcga_freqdict)
    print(exac_freqdict)

# 0.000107  = 3/len(BRCA)
# 0.0000075 = 20/len()

    # what cutoffs give largest difference between internal and external variant match?
    # define variant as NPid + position
    # how many internal matches are found in the ExAC population? What is the frequency of that variant (NPid + pos) in the ExAC pop?
    # how many of each proteins are crumby above the cutoffs in TCGA? Max frequency?

    # <dump>

    #compare tcga variants to ExAC variants to look for any variants specific to tcga
    # ext_matchcount = 0
    # for variant_tuple in tcga_variantlist:
    #     if variant_tuple in exac_variantlist:
    #         ext_matchcount += 1
