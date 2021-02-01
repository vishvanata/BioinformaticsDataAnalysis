import pandas as pd
import var_above_cutoff as vac

"""
calls abovecutoff.py and converts protein IDs to intermediary HID's for use with KnowGene & MEDICASCY
return file format: HID, % above cutoff, frequency

"""
gv_cutoff = 0.9
def create_ac_rosetta(number):

    df = vac.above_cutoff(number).reset_index()
    df2 = pd.read_csv('text_files/listgene.txt', sep= '\t')

    df3 = pd.DataFrame(data=df2, columns=['protein', 'HID'])
    dictdf3 = df3.set_index('protein').T.to_dict('list')
    for keys in dictdf3:
        dictdf3[keys] = str(dictdf3[keys])
    df.replace({'protein': dictdf3}, inplace = True)
    df['protein'] = df['protein'].str.strip("['']")

    df4 = pd.DataFrame(data=df2, columns=['protein', 'exon'])
    dictdf4 = df4.set_index('protein').T.to_dict('list')
    for keys in dictdf4:
        dictdf4[keys] = str(dictdf4[keys])
    df.replace({'exon': dictdf4}, inplace = True)
    df['exon'] = df['exon'].str.strip("['']")

    df.to_csv("vars_above_{}.txt".format(number), sep='\t', index=False)

create_ac_rosetta(gv_cutoff)
