import pandas as pd
import epex_distr_above_cutoff as da

"""
Calls distr_above_cutoff.py and converts protein IDs to intermediary HID's and exon ID's for use with KnowGene & MEDICASCY
Creates file with columns: HID, exon, % above cutoff, frequency

"""
"""
for i in range (6,10,1) and i=i/10: [Finished in 667.701s]
"""

def create_ac_rosetta(number):

    df = da.above_cutoff(number).reset_index()
    df.insert(1, "exon", df['protein'])

    df2 = pd.read_csv('text_files/listgene.txt', sep= '\t')

    df3 = pd.DataFrame(data=df2, columns=['protein', 'HID'])
    dictdf3 = df3.set_index('protein').T.to_dict('list')
    for keys in dictdf3:
        dictdf3[keys] = str(dictdf3[keys])
    df.replace({'protein': dictdf3}, inplace = True)
    df['protein'] = df['protein'].str.strip("['']")
    df.rename(columns={"protein": "HID"})

    df4 = pd.DataFrame(data=df2, columns=['protein', 'exon'])
    dictdf4 = df4.set_index('protein').T.to_dict('list')
    for keys in dictdf4:
        dictdf4[keys] = str(dictdf4[keys])
    df.replace({'exon': dictdf4}, inplace = True)
    df['exon'] = df['exon'].str.strip("['']")

    dfnice = pd.read_csv('text_files/new_translated_allinputfreq.txt', sep= '\t')
    df.insert(4, "ExAC_freq", df['exon'])
    dfnice2 = pd.DataFrame(data=dfnice, columns=['exon', 'frequency'])
    df.rename(columns={"frequency": "ExAC_freq"})


    dictdfnice2 = dfnice2.set_index('exon').T.to_dict('list')
    for keys in dictdfnice2:
        dictdfnice2[keys] = str(dictdfnice2[keys])

    df.replace({'ExAC_freq': dictdfnice2}, inplace = True)
    df['ExAC_freq'] = df['ExAC_freq'].str.strip("['']")
    df['ExAC_freq'] = df['ExAC_freq'].astype(int)

    df.to_csv("entpentpx_distrs_files/entpxentp_distr_above_{}.txt".format(number), sep='\t', index=False)

# list = [0.4]
#
# for i in list:
#     gv_cutoff = i
#     create_ac_rosetta(gv_cutoff)
