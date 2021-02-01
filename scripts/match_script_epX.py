
#turns MAF file into ENTPRISE input

openfile_path = "textfiles/TCGA.SKCM.mutect.4b7a5729-b83e-4837-9b61-a6002dce1c0a.DR-10.0.somatic.maf"
fin = open("textfiles/exac_dict.txt")
fin2 = open(openfile_path)
filename = openfile_path[10:19]

read_exac = fin.readlines()
for i in range(6):
    fin2.readline()

read_tcga = fin2.readlines()


exdict = {}
for i in range(len(read_exac)):
    read_exac[i] = read_exac[i].split("\t")
    read_exac[i][1] = read_exac[i][1].strip()
    exdict[read_exac[i][0]] = read_exac[i][1]


fout = open("output_files/eX_tcga_{}.txt".format(filename), "w")
for i in range(len(read_tcga)):
    read_tcga[i] = read_tcga[i].split("\t")

    if read_tcga[i][8] in ["Nonsense_Mutation", "Frame_Shift_Del","Frame_Shift_Ins"]:
        # print(read_tcga[i][8], read_tcga[i][71], read_tcga[i][35])
        wild = read_tcga[i][35][2:5]
        pos_var = ""
        for char in read_tcga[i][35][5:]:
            if char in "0123456789":
                pos_var += char
        if "NM_" in read_tcga[i][71]:
            try:
                if ";" in read_tcga[i][71]:
                    NPid = exdict[read_tcga[i][71][:read_tcga[i][71].find(";")]]
                else:
                    NPid = exdict[read_tcga[i][71]]
            except:
                continue
        if len(pos_var) > 1:
            fout.write(NPid + " " + pos_var + " "+ wild.upper() + "\n")
fout.close()
