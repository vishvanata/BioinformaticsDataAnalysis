
#turns MAF file into ENTPRISE input

openfile_path = "textfiles/TCGA.MESO.varscan.f62c1539-7c81-4240-ad22-b125f80686ff.DR-10.0.somatic.maf"
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

fout = open("output_files/tcga_{}.txt".format(filename), "w")
for i in range(len(read_tcga)):
    read_tcga[i] = read_tcga[i].split("\t")
    if read_tcga[i][8] == "Missense_Mutation":
        try:
            if read_tcga[i][71][0] == "N" and ";" not in read_tcga[i][71]:
                (NPid,proteinid) = (exdict[read_tcga[i][71]], read_tcga[i][35])
                (NPid2, wild, pos, mutated) = (NPid, proteinid[2:5].upper(), proteinid[5:-3] , proteinid[-3:].upper())
                if len(wild)>=1 and len(pos) >=1 and len(mutated) >=1 :
                    fout.write(NPid + " " + pos + " " + wild + " " + mutated + "\n")


            else:
                if read_tcga[i][71][0] == "N":
                    (NPid,proteinid) = (exdict[read_tcga[i][71][:read_tcga[i][71].find(";")]], read_tcga[i][35])
                    (NPid2, wild, pos, mutated) = (NPid, proteinid[2:5].upper(), proteinid[5:-3] , proteinid[-3:].upper())
                    if len(wild)>=1 and len(pos) >=1 and len(mutated) >=1 :
                        fout.write(NPid + " " + pos + " " + wild + " " + mutated + "\n")
        except:
            continue

fout.close()
