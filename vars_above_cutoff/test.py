import pandas as pd
import ac_rosetta as ac

df = pd.read_csv('var_above_{}.txt'.format(ac.gv_cutoff), sep='\t')
print(df)
