import pandas as pd

'''

Takes file and returns the percent of variants above the given threshold (e.g > 0.5) and the frequency of the variant
returns df with columns: NP_XXX, pathogenicity

'''


def above_cutoff(cutoff):
    df = pd.read_csv("text_files/final_entpxentp_out.txt", sep='\t')
    #calculate percent above cutoff for each variant
    df2 = pd.DataFrame(df.groupby('protein')['pathogenicity'].apply(lambda x: (x>cutoff).sum()/len(x)))
    df2.columns = ['% above cutoff']
    df2['frequency'] = df['protein'].value_counts()
    #df2 = df2[df2['% above cutoff']>0] #we dont want those '% above cutoff' of 0
    return df2
