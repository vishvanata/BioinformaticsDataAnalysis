import pandas as pd

'''

takes file and returns the percent of variants above the given threshold (e.g > 0.5)
return df format: NP_XXX, pathogenicity

'''

def above_cutoff(cutoff):
    df = pd.read_csv("text_files/exout.txt", sep='\t')
    df.insert (1, "exon", df['protein'])
    df = df[df['pathogenicity']>cutoff]
    return df
