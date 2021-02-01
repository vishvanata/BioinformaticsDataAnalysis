import pandas as pd
import numpy as np
import all_med as am


# can be module file for iterator
# Predicted prevalence using MEDICASCY files
# creates new df mergedStuff with common exons and sums (% above cutoff * freq) in ExAC over all exons
# run time: ~20 sec
# if has to create new entp/entpx distr file, takes ~3 min per disease


def make_prediction(diseasenum, pathcutoff, medcutoff):

    df = pd.read_csv("entpentpx_distrs_files/entpxentp_distr_above_{}.txt".format(pathcutoff), sep='\t') #pathogenicity cutoff is chosen
    df.sort_index(inplace=True)
    dfdisease = am.get_medicascy(diseasenum, medcutoff)
    df1 = pd.DataFrame()
    mergedStuff = pd.merge(df, dfdisease, on=['protein'], how='inner')

    num_variants_in_ExAC = df['ExAC_freq'].sum()
    mergedStuff['product'] = mergedStuff['% above cutoff'] * mergedStuff['ExAC_freq']
    pred = round(((mergedStuff['product'].sum())/num_variants_in_ExAC)*100, 5)

    return (diseasenum,pred,pathcutoff,medcutoff)


    # return "Disease ID: {} | Pred Prevalence: {} % | ENTP/ENTPX: {} | Know-GENE: {} | Know-Gene 0or1: {}".format(diseasenum,pred,pathcutoff,kgcutoff,args)

# make_pred(69,0.4,0)
# print(make_prediction(2227,0.5, 0.065))
